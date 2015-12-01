/** tutorial/Ex1
 * This example shows how to:
 * initialize a femus application;
 * define the multilevel-mesh object mlMsh;
 * read from the file ./input/square.neu the coarse-level mesh and associate it to mlMsh;
 * add in mlMsh uniform refined level-meshes;
 * define the multilevel-solution object mlSol associated to mlMsh;
 * add in mlSol different types of finite element solution variables;
 * initialize the solution varables;
 * define vtk and gmv writer objects associated to mlSol;
 * print vtk and gmv binary-format files in ./output directory.
 **/

#include "FemusInit.hpp"
#include "MultiLevelProblem.hpp"
#include "VTKWriter.hpp"
#include "GMVWriter.hpp"
#include "LinearImplicitSystem.hpp"
#include "NumericVector.hpp"
#include "Files.hpp"

using namespace femus;

bool SetBoundaryCondition(const std::vector < double >& x, const char solName[], double& value, const int faceName, const double time) {
  bool dirichlet = true; //dirichlet
  value = 7.*x[0];

//   if (faceName == 1)
//     dirichlet = false;

  return dirichlet;
}
// double InitalValueU(const std::vector < double >& x) {
//   return x[0] + x[1];
// }


double InitialValueU(const std::vector < double >& x) {
  return 7.*x[0];                           //boundary condition should be same as initialvalue
}

void AssemblePoissonProblem(MultiLevelProblem& ml_prob);

double GetExactSolutionLaplace(const std::vector < double >& x) {
  double pi = acos(-1.);
  return -pi * pi * cos(pi * x[0]) * cos(pi * x[1]) - pi * pi * cos(pi * x[0]) * cos(pi * x[1]);
};



int main(int argc, char** args) {

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  
  // ======= Files ========================
  Files files; 
        files.CheckIODirectories();
        files.RedirectCout();

  
  
  //=======================MESH BEGIN===================================================
  // define multilevel mesh
  MultiLevelMesh mlMsh;
  double scalingFactor = 1.;
  // read coarse level mesh and generate finers level meshes
//   mlMsh.ReadCoarseMesh("./input/square.neu", "seventh", scalingFactor);
  /* "seventh" is the order of accuracy that is used in the gauss integration scheme
      probably in the furure it is not going to be an argument of this function   */
//  const unsigned int nx, const unsigned int ny,const unsigned int nz,onst double xmin, const double xmax,st double ymin, const double ymax,onst double zmin, const double zmax, const ElemType type,
//                                const char GaussOrder[]


  mlMsh.GenerateCoarseBoxMesh( 32,32,0,-0.5,0.5,-0.5,0.5,0.,0.,QUAD9,"seventh");
  
  unsigned numberOfUniformLevels = 1;
  unsigned numberOfSelectiveLevels = 0;
  mlMsh.RefineMesh(numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL);
  mlMsh.PrintInfo();

  // define the multilevel solution and attach the mlMsh object to it
  MultiLevelSolution mlSol(&mlMsh);

  // add variables to mlSol
  mlSol.AddSolution("U", LAGRANGE, SECOND);
  

  mlSol.Initialize("U",InitialValueU);    
  // initialize all varaibles to one or any initial value with  "U",InitialValueU
  // initialize all varaibles to zero with "All"
  
 

      // attach the boundary condition function and generate boundary data
      mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
      mlSol.GenerateBdc("U");

      // define the multilevel problem attach the mlSol object to it
      MultiLevelProblem mlProb(&mlSol);

      // add system Poisson in mlProb as a Linear Implicit System
      LinearImplicitSystem& system = mlProb.add_system < LinearImplicitSystem > ("Poisson");

      // add solution "u" to system
      system.AddSolutionToSystemPDE("U");

      // attach the assembling function to system
      system.SetAssembleFunction(AssemblePoissonProblem);

      // initilaize and solve the system
      system.init();
      system.solve();


  //mlSol.Initialize("U", InitalValueU);
//    
  // print solutions
  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back("U");

  VTKWriter vtkIO(&mlSol);
  vtkIO.write(files.GetOutputPath(), "biquadratic", variablesToBePrinted);

  GMVWriter gmvIO(&mlSol);
  variablesToBePrinted.push_back("all");
  gmvIO.SetDebugOutput(false);
  gmvIO.write(files.GetOutputPath(), "biquadratic", variablesToBePrinted);

  return 0;
}

void AssemblePoissonProblem(MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data

  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  //  extract pointers to the several objects that we are going to use

  LinearImplicitSystem* mlPdeSys  = &ml_prob.get_system<LinearImplicitSystem> ("Poisson");   // pointer to the linear implicit system named "Poisson"
  const unsigned level = mlPdeSys->GetLevelToAssemble();
  const unsigned levelMax = mlPdeSys->GetLevelMax();
  const bool assembleMatrix = mlPdeSys->GetAssembleMatrix();

  Mesh*                    msh = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem*                     el = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution*    mlSol = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*                sol = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object

  LinearEquationSolver* pdeSys = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object
  SparseMatrix*             KK = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector*           RES = pdeSys->_RES; // pointer to the global residual vector object in pdeSys (level)

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  unsigned dim2 = (3 * (dim - 1) + !(dim - 1));        // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

  //solution variable
  unsigned soluIndex;
  soluIndex = mlSol->GetIndex("U");    // get the position of "u" in the ml_sol object
  unsigned soluType = mlSol->GetSolutionType(soluIndex);    // get the finite element type for "u"

  unsigned soluPdeIndex;
  soluPdeIndex = mlPdeSys->GetSolPdeIndex("U");    // get the position of "u" in the pdeSys object

  vector < double >  solu; // local solution
  solu.reserve(maxSize);

  vector < vector < double > > x(dim);    // local coordinates
  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  for (unsigned i = 0; i < dim; i++) {
    x[i].reserve(maxSize);
  }

  vector <double> phi;  // local test function
  vector <double> phi_x; // local test function first order partial derivatives
  vector <double> phi_xx; // local test function second order partial derivatives
  
  phi.reserve(maxSize);
  phi_x.reserve(maxSize * dim);
  phi_xx.reserve(maxSize * dim2);

  vector< int > l2GMap; // local to global mapping
  l2GMap.reserve(maxSize);
  
  
  //========================EQUATION=============================== 
  double weight; // gauss point weight

  vector< double > Res; // local redidual vector
  Res.reserve(maxSize);

  
  
  vector < double > Jac;
  Jac.reserve(maxSize * maxSize);

  if (assembleMatrix) {  KK->zero(); }// Set to zero all the entries of the Global Matrix
  
  
  RES->zero();// Set to zero all the entries of the Global RHS
    
  //========================EQUATION===============================  
    
  int counter = 0;
  
  
  // element loop: each process loops only on the elements that owns
  for (int iel = msh->IS_Mts2Gmt_elem_offset[iproc]; iel < msh->IS_Mts2Gmt_elem_offset[iproc + 1]; iel++) {

    unsigned kel = msh->IS_Mts2Gmt_elem[iel]; // mapping between paralell dof and mesh dof
    
    
    //====================GEOMETRY========================================
    short unsigned kelGeom = el->GetElementType(kel);    // element geometry type
    unsigned nDofx = el->GetElementDofNumber(kel, xType);    // number of coordinate element dofs
   
    for (int i = 0; i < dim; i++) {  //RESISE
      x[i].resize(nDofx);
    }
    
    for (unsigned i = 0; i < nDofx; i++) {   //FILL
      unsigned iNode = el->GetMeshDof(kel, i, xType);    // local to global coordinates node
      unsigned xDof  = msh->GetMetisDof(iNode, xType);    // global to global mapping between coordinates node and coordinate dof

      for (unsigned jdim = 0; jdim < dim; jdim++) {
        x[jdim][i] = (*msh->_coordinate->_Sol[jdim])(xDof);      // global extraction and local storage for the element coordinates
      }
    }
    
    
    //====================GEOMETRY========================================
    
    
    //====================SOLUTION========================================
    unsigned nDofu  = el->GetElementDofNumber(kel, soluType);    // number of solution element dofs
   
    // resize local arrays             //RESIZE
    l2GMap.resize(nDofu);
    solu.resize(nDofu);

    
     // local storage of global mapping and solution
    for (unsigned i = 0; i < nDofu; i++) {  //FILL
      unsigned iNode = el->GetMeshDof(kel, i, soluType);    // local to global solution node
      unsigned solDof = msh->GetMetisDof(iNode, soluType);    // global to global mapping between solution node and solution dof
      solu[i] = (*sol->_Sol[soluIndex])(solDof);      // global extraction and local storage for the solution
      l2GMap[i] = pdeSys->GetKKDof(soluIndex, soluPdeIndex, iNode);    // global to global mapping between solution node and pdeSys dof
    }
    
    //====================SOLUTION========================================
    
    
    //====================EQUATION(RHS)========================================
    
    
    Res.resize(nDofu);    //resize
    std::fill(Res.begin(), Res.end(), 0);    //set Res to zero

    Jac.resize(nDofu * nDofu);    //resize
    std::fill(Jac.begin(), Jac.end(), 0);    //set Jac to zero

     
    //====================EQUATION(RHS)=======================================
   
    

    // local storage of coordinates
    

    if (level == levelMax || !el->GetRefinedElementIndex(kel)) {      // do not care about this if now (it is used for the AMR)

      // *** Gauss point loop ***
      for (unsigned ig = 0; ig < msh->_finiteElement[kelGeom][soluType]->GetGaussPointNumber(); ig++) {
        // *** get gauss point weight, test function and test function partial derivatives ***
        msh->_finiteElement[kelGeom][soluType]->Jacobian(x, ig, weight, phi, phi_x, phi_xx);

        // evaluate the solution, the solution derivatives and the coordinates in the gauss point
        double solu_gss = 0;
        vector < double > gradSolu_gss(dim, 0.);
        vector < double > x_gss(dim, 0.);

        for (unsigned i = 0; i < nDofu; i++) {
          solu_gss += phi[i] * solu[i];

          for (unsigned jdim = 0; jdim < dim; jdim++) {
            gradSolu_gss[jdim] += phi_x[i * dim + jdim] * solu[i];
            x_gss[jdim] += x[jdim][i] * phi[i];
          }
        }

        // *** phi_i loop ***
        for (unsigned i = 0; i < nDofu; i++) {

          double laplace_mat = 0.;			//laplace is not necessary in this example
          //it is used to explain Ax=b where b=Ax0 and to solve with residue (error) method 
          //with exact solution

          for (unsigned jdim = 0; jdim < dim; jdim++) {
            laplace_mat   +=  phi_x[i * dim + jdim] * gradSolu_gss[jdim];
          }

          double srcTerm = - GetExactSolutionLaplace(x_gss);
          Res[i] += (srcTerm * phi[i]) * weight;

          if (assembleMatrix) {
            // *** phi_j loop ***
            for (unsigned j = 0; j < nDofu; j++) {
                            double laplace_mat = 0.;
	      counter++;
	      
	      
              for (unsigned kdim = 0; kdim < dim; kdim++) {
                laplace_mat += (phi_x[i * dim + kdim] * phi_x[j * dim + kdim]) * weight;
		

              }

              Jac[i * nDofu + j] += laplace_mat;
            } // end phi_j loop
          } // endif assemble_matrix

        } // end phi_i loop
      } // end gauss point loop
    } // endif single element not refined or fine grid loop
    
      
    
    
    //--------------------------------------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector

    //copy the value of the adept::adoube aRes in double Res and store
    RES->add_vector_blocked(Res, l2GMap);

    if (assembleMatrix) {
      //store K in the global matrix KK
      KK->add_matrix_blocked(Jac, l2GMap, l2GMap);
    }
  } //end element loop for each process
  
  
  std::cout << "counter value ................  " << counter << std::endl;

  RES->close();

  if (assembleMatrix) KK->close();

  // ***************** END ASSEMBLY *******************
}



