#include "allen.h" 
#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm> 
#include <iomanip> 
#include <string>

// Global parameter: Interface thickness
double epsilon = 0.01;     

// Initial condition
double _u_(const double* p)
{
    double x = p[0];
    double y = p[1];
    
    double term1 = (x * x) / 0.04 + (y * y) / 0.36 - 1.0;
    
    double term2 = (x * x) / 0.36 + (y * y) / 0.04 - 1.0;

    double val = 5.0 * term1 * term2;

    return std::tanh(val);
}

// AVF nonlinear stabilization term corresponding to Eq. (2.10)
inline double AVF_nonlinear(double u_k, double u_n)
{
    return 0.25 * (u_k + u_n) * (u_k * u_k + u_n * u_n - 2.0);
}

// Calculate local element matrix: M / dt + 0.5 * K
void Matrix::getElementMatrix(const Element<double, DIM>& ele0,
                              const Element<double, DIM>& ele1,
                              const ActiveElementPairIterator<DIM>::State state)
{
    // Calculate the volume of the reference element
    double volume = ele0.templateElement().volume();

    // Retrieve the quadrature rule of the reference element.
    const QuadratureInfo<DIM>& quad_info = ele0.findQuadratureInfo(algebricAccuracy());
    
    // Calculate the Jacobian of the coordinate transformation at the quadrature points
    std::vector<double> jacobian = ele0.local_to_global_jacobian(quad_info.quadraturePoint());
    
    // Number of quadrature points
    int n_quad_p = quad_info.n_quadraturePoint();

    // Transform quadrature points from reference to physical coordinates
    std::vector<AFEPack::Point<DIM>> q_point = ele0.local_to_global(quad_info.quadraturePoint());
    
    // Evaluate basis function values and gradients at quadrature points
    std::vector<std::vector<double>> bas_val = ele0.basis_function_value(q_point);
    std::vector<std::vector<std::vector<double>>> bas_grad = ele0.basis_function_gradient(q_point);
    
    // Number of degrees of freedom per element
    int n_ele_dof = ele0.dof().size();

    for (int l = 0; l < n_quad_p; l++) {
        double Jxw = quad_info.weight(l) * jacobian[l] * volume;
        for (int j = 0; j < n_ele_dof; j++) {
            for (int k = 0; k < n_ele_dof; k++) {
                elementMatrix(j, k) += Jxw * (
                    // M / dt
                    bas_val[j][l] * bas_val[k][l] / _dt + 
                    // 0.5 * K
                    0.5 * innerProduct(bas_grad[j][l], bas_grad[k][l]) 
                );
            }
        }
    }
}

void allen::initialize()
{

    //Read and process the mesh
    h_tree.readEasyMesh("D");
    ir_mesh = new IrregularMesh<DIM>(h_tree);
    ir_mesh->globalRefine(2);
    ir_mesh->semiregularize();
    ir_mesh->regularize(false);
    
    //Configure reference elements
    template_geometry.readData("triangle.tmp_geo");
    coord_transform.readData("triangle.crd_trs");
    template_dof.reinit(template_geometry);
    template_dof.readData("triangle.1.tmp_dof");
    basis_function.reinit(template_dof);
    basis_function.readData("triangle.1.bas_fun");
    
    template_geometry1.readData("twin_triangle.tmp_geo");
    coord_transform1.readData("twin_triangle.crd_trs");
    template_dof1.reinit(template_geometry1);
    template_dof1.readData("twin_triangle.1.tmp_dof");
    basis_function1.reinit(template_dof1);
    basis_function1.readData("twin_triangle.1.bas_fun");
    
    template_element.resize(2);
    template_element[0].reinit(template_geometry, template_dof, coord_transform, basis_function);
    template_element[1].reinit(template_geometry1, template_dof1, coord_transform1, basis_function1);
    
    buildFEMSpace();
    
    // Initialize time parameters
    t = 0.0; 
    dt = 1.0e-05; 
}

//Build the finite element space
void allen::buildFEMSpace()
{
    RegularMesh<DIM>& mesh = ir_mesh->regularMesh();
    u_int n_ele = mesh.n_geometry(DIM);
    fem_space = new FEMSpace<double, DIM>(mesh, template_element);
    fem_space->element().resize(n_ele);
    for (u_int i = 0; i < n_ele; ++i) {
        u_int n_vtx = mesh.geometry(DIM, i).n_vertex();
        if (n_vtx == 3) {
            fem_space->element(i).reinit(*fem_space, i, 0);
        } else {
            fem_space->element(i).reinit(*fem_space, i, 1);
        }
    }
    fem_space->buildElement();
    fem_space->buildDof();
    fem_space->buildDofBoundaryMark();

    // Initialize the FEM function 
    u_h = new FEMFunction<double, DIM>(*fem_space);
    Operator::L2Interpolate(_u_, *u_h);
}

// Compute the a posteriori error indicator
void allen::getIndicator(Indicator<DIM>& ind)
{

    RegularMesh<DIM>& mesh = ir_mesh->regularMesh();
    u_int n_face = mesh.n_geometry(DIM - 1);
    std::vector<bool> flag(n_face, false);
    std::vector<double> jump(n_face);

    // Calculate the jump of normal derivatives across each interior face
    FEMSpace<double, DIM>::ElementIterator
        the_ele = fem_space->beginElement(),
        end_ele = fem_space->endElement();
    for (u_int i = 0; the_ele != end_ele; ++the_ele, ++i) {
        AFEPack::Point<DIM> p;
        std::vector<double> u_h_grad = u_h->gradient(p, *the_ele);
        GeometryBM& geo = the_ele->geometry();
        u_int n_bnd = geo.n_boundary();
        for (u_int j = 0; j < n_bnd; ++j) {
            u_int sid_idx = geo.boundary(j);
            GeometryBM& bnd = mesh.geometry(DIM - 1, sid_idx);
            AFEPack::Point<DIM>& p0 = mesh.point(bnd.vertex(0));
            AFEPack::Point<DIM>& p1 = mesh.point(bnd.vertex(1));
            double grad_normal = (u_h_grad[0] * (p1[1] - p0[1]) -
                                  u_h_grad[1] * (p1[0] - p0[0]));
            if (flag[sid_idx] == false) {
                jump[sid_idx] = grad_normal;
                flag[sid_idx] = true;
            } else {
                jump[sid_idx] -= grad_normal;
                flag[sid_idx] = false;
            }
        }
    }

     // Aggregate jumps to compute the error indicator for each element
    the_ele = fem_space->beginElement();
    end_ele = fem_space->endElement();
    for (u_int i = 0; the_ele != end_ele; ++the_ele, ++i) {
        GeometryBM& geo = the_ele->geometry();
        u_int n_bnd = geo.n_boundary();
        ind[i] = 0.0;
        for (u_int j = 0; j < n_bnd; ++j) {
            u_int sid_idx = geo.boundary(j);
            if (flag[sid_idx]) continue;
            ind[i] += jump[sid_idx] * jump[sid_idx];
        }
    }
}

void allen::adaptMesh()
{
    Indicator<DIM> ind(ir_mesh->regularMesh());
    getIndicator(ind);
    
    IrregularMesh<DIM>* old_ir_mesh = ir_mesh;
    FEMSpace<double, DIM>* old_fem_space = fem_space;
    FEMFunction<double, DIM>* old_u_h = u_h;
    
    u_h = nullptr; 
    fem_space = nullptr; 

    ir_mesh = new IrregularMesh<DIM>(*old_ir_mesh);
    MeshAdaptor<DIM> mesh_adaptor(*old_ir_mesh, *ir_mesh);
    mesh_adaptor.setIndicator(ind);


    // Adaptivity control parameters
    mesh_adaptor.tolerence() = 5.0e-3;
    mesh_adaptor.convergenceOrder() = 1;
    mesh_adaptor.refineStep() = 1; 
    mesh_adaptor.adapt();
    
    ir_mesh->semiregularize();
    ir_mesh->regularize(false);

    buildFEMSpace();
    Operator::L2Interpolate(*old_u_h, *u_h);
    
    delete old_u_h;
    delete old_fem_space;
    delete old_ir_mesh;
}

void allen::stepForward()
{
    FEMFunction<double, DIM> last_u_h(*u_h); 
    FEMFunction<double, DIM> u_k(*u_h);      

    // Prepare the linear system matrix
    Matrix stiff_matrix(*fem_space, dt);
    stiff_matrix.algebricAccuracy() = 3;
    stiff_matrix.build(); 

    int nonlinear_iter = 0;
    double error = 0.0; 
    int max_iter = 100;

    // Start nonlinear iteration loop
    for (; nonlinear_iter < max_iter; ++nonlinear_iter) { 
        
        // Prepare the right-hand side vector
        Vector<double> rhs;
        rhs.reinit(fem_space->n_dof());
        rhs = 0.0;
        
        FEMSpace<double, DIM>::ElementIterator the_ele = fem_space->beginElement();
        FEMSpace<double, DIM>::ElementIterator end_ele = fem_space->endElement();

        for (; the_ele != end_ele; ++the_ele) {

            double vol = the_ele->templateElement().volume();
            const QuadratureInfo<DIM>& quad_info = the_ele->findQuadratureInfo(3);
            int n_qp = quad_info.n_quadraturePoint();

            std::vector<double> jac = the_ele->local_to_global_jacobian(quad_info.quadraturePoint());

            std::vector<AFEPack::Point<DIM>> q_p = the_ele->local_to_global(quad_info.quadraturePoint());

            std::vector<std::vector<double>> bas_val = the_ele->basis_function_value(q_p);

            std::vector<std::vector<std::vector<double>>> bas_grad = the_ele->basis_function_gradient(q_p);

            // The solution and gradient from the previous step
            std::vector<double> u_n_val = last_u_h.value(q_p, *the_ele);
            
            std::vector<std::vector<double>> u_n_grad = last_u_h.gradient(q_p, *the_ele);

            // The solution at the current step
            std::vector<double> u_k_val = u_k.value(q_p, *the_ele);

            const std::vector<int>& ele_dof = the_ele->dof();

            for (int l = 0; l < n_qp; ++l) {

                double Jxw = vol * jac[l] * quad_info.weight(l);

                double nl_term = AVF_nonlinear(u_k_val[l], u_n_val[l]);

                for (u_int j = 0; j < ele_dof.size(); ++j) {

                    double grad_term = 0.0;
                    for (int d = 0; d < DIM; ++d)
                        grad_term += u_n_grad[l][d] * bas_grad[j][l][d];
                        
                        rhs(ele_dof[j]) += Jxw * ((u_n_val[l] / dt) * bas_val[j][l] - 0.5 * grad_term - (nl_term / (epsilon * epsilon)) * bas_val[j][l]);
                }
            }
        }
        
         // Solve the linear system using an Algebraic Multigrid (AMG) solver
        AMGSolver solver;
        solver.lazyReinit(stiff_matrix);
        solver.solve(*u_h, rhs); 
        
        // Check convergence by calculating the L2 norm of the difference between iterations
        FEMFunction<double, DIM> diff(*u_h);
        diff.add(-1.0, u_k);
        error = Functional::L2Norm(diff, 3);

        u_k = *u_h; 
        
        // Check if convergence criterion is met
        if (error < 1e-8) {
            std::cout << "Converged. Iterations: "
                      << nonlinear_iter + 1
                      << ", Error: "
                      << error << std::endl;
            break;
        }
    } 

    if (nonlinear_iter == max_iter) {
        std::cout << "WARNING: Nonlinear iteration "
                  << "failed to converge! "
                  << "Final error: "
                  << error << std::endl;
    }
}

// Compute energy
double  allen::calculateEnergy()
{
    double energy = 0.0;
    
    FEMSpace<double, DIM>::ElementIterator the_ele = fem_space->beginElement();
    FEMSpace<double, DIM>::ElementIterator end_ele = fem_space->endElement();

    for (; the_ele != end_ele; ++the_ele) {
        double vol = the_ele->templateElement().volume();
        const QuadratureInfo<DIM>& quad_info = the_ele->findQuadratureInfo(3);
        int n_qp = quad_info.n_quadraturePoint();
        
        std::vector<double> jac = the_ele->local_to_global_jacobian(quad_info.quadraturePoint());
        std::vector<AFEPack::Point<DIM>> q_p = the_ele->local_to_global(quad_info.quadraturePoint());
        
        std::vector<double> u_val = u_h->value(q_p, *the_ele);
        std::vector<std::vector<double>> u_grad = u_h->gradient(q_p, *the_ele);

        for (int l = 0; l < n_qp; ++l) {
            double Jxw = vol * jac[l] * quad_info.weight(l);

            double u = u_val[l];

            double grad_sq = 0.0;

            // Calculate squared magnitude of the gradient: |\nabla u|^2
            for(int d=0; d<DIM; ++d) {

                grad_sq += u_grad[l][d] * u_grad[l][d];

            }

             // Gradient energy density: 0.5 * |\nabla u|^2
            double energy_grad = 0.5 * grad_sq;

            // Potential energy density: (u^2 - 1)^2 / (4 * epsilon^2)
            double term = (u * u - 1.0);
            double energy_pot = 0.25 * term * term / (epsilon * epsilon);

            energy += (energy_grad + energy_pot) * Jxw;
        }
    }

    return energy;
}

// Compute the L2 difference norm between two FEM functions
double calc_L2_diff(const FEMFunction<double, DIM>& u1, 
                    const FEMFunction<double, DIM>& u2)
{
    FEMFunction<double, DIM> diff(u1);
    diff.add(-1.0, u2); 
    return Functional::L2Norm(diff, 3);
}

void allen::adaptTimeStep()
{

    // Set tolerances and parameters for adaptive time stepping
    double tol_time_max = 1.0e-2; 
    double tol_time_min = 1.0e-3; 
    double delta1 = 0.5;          
    double delta2 = 2.0;          
    double dt_min = 1.0e-8;
    double dt_max = 1.0e-5; 
    
    //u^{n-1} 
    FEMFunction<double, DIM> u_prev(*u_h); 

    bool step_accepted = false;

    while (!step_accepted ) {

        stepForward(); 

        double eta = calc_L2_diff(*u_h, u_prev);

        // Case 1: Error is too large, reduce dt and recompute
        if (eta > tol_time_max) {
            dt *= delta1;

            // If dt hits minimum limit, accept the step to avoid infinite loop
            if (dt < dt_min) {
                std::cout << "Warning: dt reached min limit (" << dt_min << "). Accepting step." << std::endl;
                step_accepted = true;
            } else {
                *u_h = u_prev;
            }
        }

        // Case 2: Error is very small, increase dt to improve efficiency
        else if (eta < tol_time_min) {

            // If dt reaches maximum limit, accept the step
            if (dt >= dt_max) {

                step_accepted = true; 

            } else {

                dt *= delta2;

                if (dt > dt_max) dt = dt_max; 

                *u_h = u_prev;
            }
        }
        else {
            step_accepted = true;
        }
    }
    
    if (!step_accepted) {
        std::cout << "Warning: Adaptive loop limit reached. Step accepted forcedly." << std::endl;
    }
}
void allen::run()
{
    initialize();
    Operator::L2Interpolate(&_u_, *u_h);

    // Initial mesh adaptation
    std::cout << "Adapting initial mesh..." << std::endl;
    for (int i = 0; i <10 ; ++i) {
        adaptMesh(); 
        Operator::L2Interpolate(_u_, *u_h); 
        std::cout << "Initial Adaptation Step " << i + 1 
                  << ", DOFs: " << fem_space->n_dof() << std::endl;
    }

    u_h->writeOpenDXData("eg32_0.dx"); 

    // Set up the file to log energy data
    std::ofstream eng_file("energy.dat");
    eng_file << "# Time \t Energy \t dt \t DOFs \t Elements \t EnergyDiff" << std::endl;
    eng_file << std::scientific << std::setprecision(8);

    double last_energy = calculateEnergy(); 
    double current_energy = last_energy;
    
    int n_dof = fem_space->n_dof();
    int n_ele = fem_space->element().size();

    // Log initial energy state (Step 0)
    eng_file << t << "\t" << last_energy << "\t" << dt << "\t" << n_dof << "\t" << n_ele << "\t" << 0.0 << std::endl; 

    int step = 0;
    
    double energy_tol = 1.0e-06; 
    
    double energy_diff = 1.0; 
    
    double T_MAX = 1.0; 

    std::cout << "Starting Simulation: Run until Energy Difference <= " << energy_tol << std::endl;
    std::cout << "Step \t Time \t\t dt \t\t Energy \t E_Diff" << std::endl;

    // Main evolution loop: continue until energy converges or time reaches T_MAX
    while (energy_diff > energy_tol && t < T_MAX) {
        step++;

        adaptTimeStep();

        t += dt; 
        
        current_energy = calculateEnergy();
        
        energy_diff = std::abs(current_energy - last_energy);

        n_dof = fem_space->n_dof();
        n_ele = fem_space->element().size();

        eng_file << t << "\t" 
                 << current_energy << "\t" 
                 << dt << "\t" 
                 << n_dof << "\t" 
                 << n_ele << "\t" 
                 << energy_diff << std::endl;

        if (step % 1 == 0) { 
            adaptMesh();
        }

        if (step % 10 == 0) {
            std::cout << step << "\t " 
                      << t << "\t " 
                      << dt << "\t " 
                      << current_energy << "\t " 
                      << energy_diff 
                      << std::endl;

            std::string filename = "eg32_" + std::to_string(t) + ".dx";
            u_h->writeOpenDXData(filename.c_str());
        }

        last_energy = current_energy;
    }
    
    eng_file.close();
}
