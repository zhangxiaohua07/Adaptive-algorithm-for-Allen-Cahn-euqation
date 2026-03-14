#ifndef _allen_h  
#define _allen_h

#include <AFEPack/AMGSolver.h>               
#include <AFEPack/TemplateElement.h>         
#include <AFEPack/FEMSpace.h>                
#include <AFEPack/Operator.h>                
#include <AFEPack/BilinearOperator.h>        
#include <AFEPack/Functional.h>              
#include <AFEPack/EasyMesh.h>               

#define DIM 2
#define PI (4.0*atan(1.0))

class allen {
private:

    TemplateGeometry<DIM> template_geometry;    
    CoordTransform<DIM, DIM> coord_transform;   
    TemplateDOF<DIM> template_dof;             
    BasisFunctionAdmin<double, DIM, DIM> basis_function; 

    TemplateGeometry<DIM> template_geometry1;   
    CoordTransform<DIM, DIM> coord_transform1;  
    TemplateDOF<DIM> template_dof1;             
    BasisFunctionAdmin<double, DIM, DIM> basis_function1; 

    std::vector<TemplateElement<double, DIM> > template_element; 

    HGeometryTree<DIM> h_tree;          
    IrregularMesh<DIM>* ir_mesh;        

   
    FEMSpace<double, DIM>* fem_space;   
    FEMFunction<double, DIM>* u_h;   

    double t;  // current time
    double dt; // time step size

public:
    void initialize();   // initialization
    void run();          // main loop
    void stepForward();  // single time step evolution
    void buildFEMSpace();// construct the finite element space
    void getIndicator(Indicator<DIM>& ind); // compute the adaptive indicator
    void adaptMesh();    // mesh adaptivity
    void adaptTimeStep(); // adaptive time stepping
    double calculateEnergy(); // compute energy
};


class Matrix : public L2InnerProduct<DIM, double>                                  
{
    private:
    double _dt;

    public:
    Matrix(FEMSpace<double, DIM> & sp,double dt) : L2InnerProduct<DIM, double>(sp,sp), _dt(dt) {}; 
    virtual ~Matrix() {};   
    public:
    virtual void getElementMatrix(const Element<double, DIM> & ele0,
                                  const Element<double, DIM> & e1e1,
                                  const ActiveElementPairIterator<DIM>::State state);
};

#endif