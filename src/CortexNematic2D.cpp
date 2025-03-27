
/// C++ headers
#include <iostream>
#include <mpi.h>
#include <fstream>

/// Trilinos headers
#include <Teuchos_RCP.hpp>

/// hiperlife headers
#include "hl_TypeDefs.h"
#include "hl_Geometry.h"
#include "hl_StructMeshGenerator.h"
#include "hl_DistributedMesh.h"
#include "hl_FillStructure.h"
#include "hl_DOFsHandler.h"
#include "hl_HiPerProblem.h"
#include "hl_LinearSolver_Direct_Amesos.h"
#include "hl_Tensor.h"
#include "hl_MeshLoader.h"
#include <Teuchos_RCP.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include "hl_ConfigFile.h"
#include "hl_ConsistencyCheck.h"
#include "hl_ConsistencyCheck.h"
#include "AuxCortexNematic2D.h"
#include "hl_LinearSolver_Direct_Amesos.h"
#include "hl_LinearSolver_Iterative_AztecOO.h"
#include "hl_NonlinearSolver_NewtonRaphson.h"
#include "hl_LinearSolver_Direct_Amesos2.h"
#include "hl_LinearSolver_Direct_MUMPS.h"

int main ( int argc, char *argv[]  )
{
    using namespace std;
    using namespace hiperlife;
    using Teuchos::rcp;
    using Teuchos::RCP;


    //Initialize MPI
    MPI_Init(nullptr, nullptr);

    bool printVtk{true};
    bool printFile{false};
    //bool diagnosticFlag{false};

    const char *config_filename = argv[1];
    ConfigFile config(config_filename);

    double restart{};
    config.readInto(restart, "restart");

    // path to store the .csv file for output
    string csvPathDiss = "tensionOnBoundary.csv" ;
    config.readInto(csvPathDiss, "csvPathDiss");


    //Number of Gauss points
    int gPts{};
    config.readInto(gPts, "gPts");

    // Geoemtric parameters for structured mesh
    double r_in = 3.0;
    config.readInto(r_in, "r_in");

    double r_out = 10.0;
    config.readInto(r_out, "r_out");

    double phi = 2.*M_PI;
    config.readInto(phi, "phi");

    int nx=100;
    config.readInto(nx, "nx");

    int ny=20;
    config.readInto(ny, "ny");

    //Model parameters
    double visc      = 20.0;
    config.readInto(visc, "visc");

    double cvisc     = 3.0;
    config.readInto(cvisc, "cvisc");

    double rvisc = 1.0;
    config.readInto(rvisc, "rvisc");

    double frank   = 1.0;
    config.readInto(frank, "frank");

    double lambda_iso   = 0.0;//10.0;  // Contractile case
    config.readInto(lambda_iso, "lambda_iso");

    double lambda_aniso = 0.0;//10.0; //10.0;
    config.readInto(lambda_aniso, "lambda_aniso");

    double lambda_rot = 0.0;//10.0;
    config.readInto(lambda_rot, "lambda_rot");

    double kp    = 0.2*0.1;
    config.readInto(kp, "kp");

    double kd    = 0.1;
    config.readInto(kd, "kd");

    double fric    = 0.0;
    config.readInto(fric, "fric");

    double sus20    = 200.0;
    config.readInto(sus20, "sus20");

    double sus40    = 800.0;
    config.readInto(sus40, "sus40");

    double hcrit    = 0.25;
    config.readInto(hcrit, "hcrit");

    //Time integrator
    double deltat  = 1E-2;
    config.readInto(deltat, "deltat");

    //Numerical parameters for linear and NR solver
    string linTol = "1e-12";
    config.readInto(linTol, "linTol");

    double resTol = 1e-12;
    config.readInto(resTol, "resTol");

    double solTol = 1e-12;
    config.readInto(solTol, "solTol");

    int maxIter = 10;
    config.readInto(maxIter, "maxIter");

    string linImax = "10000";
    config.readInto(linImax, "linImax");

    double maxDelt = 1.0;
    config.readInto(maxDelt, "maxDelt");
    //testCase=1 for the annulus
    //testCase=2 for full circle
    //testCase=3 for square domain
    int testCase = 3;
    config.readInto(testCase, "testCase");

    double kosm = 0.001;
    config.readInto(kosm, "kosm");

    int totalTimeSteps{1000000};
    config.readInto(totalTimeSteps, "totalTimeSteps");


    double adaptiveStepTime = 0.25;
    config.readInto(adaptiveStepTime, "adaptiveStepTime");


    int nPrint = 1;
    config.readInto(nPrint, "nPrint");

    double stab = 1.0;
    config.readInto(stab, "stab");

    double uPoly = -5.0/60.0;
    config.readInto(uPoly, "uPoly");

    double t_stall = -1.0;
    config.readInto(t_stall, "t_stall");

    //Consistency check
    bool cCheck{true};
    config.readInto(cCheck, "cCheck");

    double heqb = 0.20;
    config.readInto(heqb, "heqb");

    double noiseStrength = 0.001;
    config.readInto(noiseStrength, "noiseStrength");   


    double tol = .01; 
    config.readInto(tol, "tol");   


    double crossLinker = 20;
    config.readInto(crossLinker, "crossLinker");   


    double k_on = 0.02;
    config.readInto(k_on,"k_on");

    double k_off = 0.1;
    config.readInto(k_off,"k_off"); 


    double strainGradient = -3; //9;
    config.readInto(strainGradient,"strainGradient"); 

    int nIter = 2;
    
    ofstream fileNameDiss(csvPathDiss);
    fileNameDiss.close();

    //Input
    string fMesh{};
    config.readInto(fMesh, "fMesh");

    //Output
    string solname_v;
    string oname = "stressFiberCircle";
    config.readInto(oname, "oname");
    int timeStep=0;
    double time = 0.0;
    //Parameters of the mesh and geometry
    int bfOrder = 1;
    bool balanceMesh = true ;
    bool testCaseMotorsTurnOff = false;
    int testCaseStrainTest=0;
    //Add model parameters to userStr
    RCP<UserStructure> userStr = rcp(new UserStructure);
    RCP<UserStructure> userNum = rcp(new UserStructure);
    {
        userStr->dparam.resize(50);

        userStr->dparam[0] = deltat;
        userStr->dparam[1] = sus20/2;
        userStr->dparam[2] = sus40/8;
        userStr->dparam[3] = hcrit;

        userStr->dparam[4] = kp;
        userStr->dparam[5] = kd;
        userStr->dparam[6] = frank;

        userStr->dparam[7] = lambda_iso;
        userStr->dparam[8] = lambda_aniso;
        userStr->dparam[9] = lambda_rot;

        userStr->dparam[10] = visc;
        userStr->dparam[11] = rvisc;
        userStr->dparam[12] = cvisc;

        userStr->dparam[13] = fric;
        userStr->dparam[14] = stab;
        userStr->dparam[15] = phi;
        userStr->dparam[16] = r_in;
        userStr->dparam[17] = r_out;
        userStr->dparam[18] = uPoly;
        userStr->dparam[19] = t_stall;
        userStr->dparam[22] = heqb;
        userStr->dparam[23] = nx;
        userStr->dparam[24] = noiseStrength;
        userStr->dparam[25] = tol;
        userStr->dparam[26] = crossLinker;
        userStr->dparam[27] = k_on;
        userStr->dparam[28] = k_off;
        userStr->dparam[30] = testCase;
        userStr->dparam[31] = strainGradient;
        userStr->dparam[32] = testCaseStrainTest;
        //Numerical parameters
        userNum->dparam.resize(3);
        userNum->iparam.resize(3);
        userNum->dparam[0] = solTol;
        userNum->dparam[1] = resTol;
        userNum->iparam[0] = maxIter;
        userNum->iparam[1] = testCase;
    
    }
    RCP<StructMeshGenerator> structMesh = rcp (new StructMeshGenerator);
    RCP<MeshLoader> loadedMesh = rcp(new MeshLoader);


    RCP<DistributedMesh> disMesh;

    if (testCase==1)
    {
        structMesh->setBasisFuncType(BasisFuncType::Lagrangian);
        structMesh->setBasisFuncOrder(bfOrder);
        structMesh->setElemType(ElemType::Square);
        // Impose periodicity in the theeta direction
        vector<Axis> axis;
        axis.push_back(Axis::Yaxis);
        structMesh->setPeriodicBoundaryCondition(axis);
        //Generate an annulus
        structMesh->genAnnularSector(nx, ny, 6., 25., phi);
      //  structMesh->stretchElemsAsymmetric(2,1,1);
    }
    else if (testCase==2)
    {


        loadedMesh->setElemType(ElemType::Square);
        loadedMesh->setBasisFuncType(BasisFuncType::Lagrangian);
        loadedMesh->setBasisFuncOrder(bfOrder);
        loadedMesh->loadLegacyVtk(fMesh);
        loadedMesh->transformFree([](double x, double y) {
            x = 5 * x;
            y = 5 * y;
            return std::make_tuple(x, y);
        });
 
    }
    else if (testCase==3)
    {

        structMesh->setBasisFuncType(BasisFuncType::Lagrangian);
        structMesh->setBasisFuncOrder(bfOrder);
        structMesh->setElemType(ElemType::Square);
        // Impose periodicity in the theeta direction
        vector<Axis> axis;
        axis.push_back(Axis::Xaxis);
        if (testCaseStrainTest==0)
        {
            axis.push_back(Axis::Yaxis);
        }
        structMesh->setPeriodicBoundaryCondition(axis);
        //Generate an annulus
        structMesh->genRectangle(nx, ny, r_out,r_out);                            
    }
    try
    {
        disMesh = rcp(new DistributedMesh);

        if (testCase==1 or testCase==3)
            disMesh->setMesh(structMesh); 
        else if (testCase==2)
            disMesh->setMesh(loadedMesh);

        
        disMesh->setBalanceMesh(balanceMesh);
        disMesh->Update();
        disMesh->printFileLegacyVtk("DensityDependCortexNemacMesh");
        if (disMesh->myRank() == 0)
            cout << "Dismesh successfully created. " << endl;
    }
    catch (runtime_error)
    {
        throw runtime_error("Error in dismesh.");
        MPI_Finalize();
        return 1;
    }

    //DOFsHand
    RCP<DOFsHandler> dofHand;

    try
    {
        dofHand = rcp(new DOFsHandler(disMesh));
        dofHand->setNameTag("cortexHand");
        dofHand->setDOFs({"h", "q1", "q2", "vx", "vy"});
        dofHand->setNodeAuxF({"t","Tx", "Ty"});
        dofHand->Update();                                      

        if (dofHand->myRank() == 0)
            cout << "DOFsHandler successfully created. " << endl;
    }
    catch (runtime_error err)
    {
        cout << dofHand->myRank() << ":   DOFsHandler could not be created. " << err.what() << endl;
        MPI_Finalize();
        return 1;
    }


    // Initial and boundary conditions
    if (testCase==1)
    {
//        setCircleVelBC(dofHand, userStr);
//        setCircleNemBC(dofHand, userStr);
//        setCircleThiBC(dofHand, userStr);
    }
    else if (testCase==2)
    {

//        setCircleVelBC_uniform(dofHand, userStr);
//        setCircleNemBC_uniform(dofHand, userStr);
//        setCircleThiBC_uniform(dofHand, userStr);
    }
    else if (testCase==3)
    {
        setCircleNemBC_square(dofHand, userStr);
        setCircleThiBC_square(dofHand, userStr);
        setCircleVelBC_square(dofHand, userStr);
    }

    dofHand->nodeDOFs0->setValue(dofHand->nodeDOFs);
    dofHand->UpdateGhosts();
    //Hiperproblem
    RCP<HiPerProblem> hiperProbl;

    hiperProbl = rcp(new HiPerProblem);
    hiperProbl->setUserStruct(userStr);
    hiperProbl->setDOFsHandlers({dofHand});
    hiperProbl->setUserStruct(userStr);
    hiperProbl->setConsistencyCheckTolerance(1.E-7);
    hiperProbl->setConsistencyCheckDelta(1.E-7);
    // Only for the annulus sector
 //   if (testCase==1)
 //       setLinearConstraints(userStr->dparam, hiperProbl);

    //Set integration
    hiperProbl->setIntegration("Integ", {"cortexHand"});
    hiperProbl->setIntegration("BorderInteg",{"cortexHand"});
    hiperProbl->setCubatureBorderGauss("BorderInteg",2);
    hiperProbl->setElemFillings("BorderInteg",nullptr,nullptr,LS_CortexNematic_Border);
    hiperProbl->setGlobInteg({"biLinear", "sigma_xx","sigma_yy","sigma_partial_xx","sigma_partial_yy","gradV_xx","gradV_yy","fric_x","fric_y","Tx","Ty"});


    // annulus with a structured square mesh
    if (testCase==1 or testCase==3) 
        hiperProbl->setCubatureGauss("Integ", 4);
    else if (testCase==2)
        hiperProbl->setCubatureGauss("Integ",  4);

    //Element fillings
    if (cCheck == true)
    {
        hiperProbl->setElemFillings("Integ", ConsistencyCheck<LS_CortexNematic_2D>);
    }
    else
    {
        hiperProbl->setElemFillings("Integ", LS_CortexNematic_2D);
    }


    //Update
    hiperProbl->Update();                                                                                                                                                 

    if (hiperProbl->myRank() == 0)
        cout << "HiperProblem successfully updated." << endl;




    //Output initial condition
    string solName = oname + ".0";
    solname_v = oname + "." + to_string(0);
    if (printVtk)
        dofHand->printFileLegacyVtk(solname_v, true);
    if (printFile)
        dofHand->printFile("InitialMesh", OutputMode::Text, true, 0.0);


    RCP<AztecOOIterativeLinearSolver> linSolver = rcp ( new AztecOOIterativeLinearSolver());
    linSolver->setHiPerProblem(hiperProbl);
    linSolver->setTolerance(1e-10);
    linSolver->setMaxNumIterations(7000);
    linSolver->setSolver(AztecOOIterativeLinearSolver::Solver::Gmres);
    linSolver->setPrecond(AztecOOIterativeLinearSolver::Precond::DomainDecomp);
    linSolver->setSubdomainSolve(AztecOOIterativeLinearSolver::SubdomainSolve::Ilut);
    linSolver->setVerbosity(AztecOOIterativeLinearSolver::Verbosity::None);
    linSolver->setDefaultParameters();
    linSolver->Update();                                                                                                                                       


//    RCP<Amesos2DirectLinearSolver> linSolver = rcp(new Amesos2DirectLinearSolver());
//    linSolver->setHiPerProblem(hiperProbl);
//    linSolver->setSolver(Amesos2DirectLinearSolver::Solver::Mumps);
//    linSolver->setUseTranspose(false);
//    linSolver->setSymbolicFactorization(true);
//    linSolver->setNumericFactorization(true);
//    linSolver->setVerbosity(AztecOOIterativeLinearSolver::Verbosity::None);
//    linSolver->setScaleMethod(true);
//    linSolver->setRefactorize(false);
//    linSolver->Update();
                                                                                                       
 
//  RCP<MUMPSDirectLinearSolver> linSolver = rcp(new MUMPSDirectLinearSolver());
//  linSolver->setHiPerProblem(hiperProbl);
//  linSolver->setMatrixType(MUMPSDirectLinearSolver::MatrixType::General); //General, POD or Symmetric
//  linSolver->setAnalysisType(MUMPSDirectLinearSolver::AnalysisType::Parallel); //Sequential or parallel
//  linSolver->setOrderingLibrary(MUMPSDirectLinearSolver::OrderingLibrary::Auto); //Many options here
//  linSolver->setVerbosity(AztecOOIterativeLinearSolver::Verbosity::None);
//  linSolver->setDefaultParameters();
//  linSolver->Update(); 

                                                                                                  

    //Set non-lienar solver
    RCP<NewtonRaphsonNonlinearSolver> nonlinSolver =  rcp ( new NewtonRaphsonNonlinearSolver());
    nonlinSolver->setLinearSolver(linSolver);
    nonlinSolver->setMaxNumIterations(maxIter);
    nonlinSolver->setResTolerance(1E-8);
    nonlinSolver->setSolTolerance(1E-8);
    nonlinSolver->setLineSearch(false);
    nonlinSolver->setPrintIntermInfo(true);
    nonlinSolver->setConvRelTolerance(false);
    nonlinSolver->Update();                                                                                  
                                                                                                            
    std::ofstream myfile;
    myfile.open ("bilinear.csv");  //Loop over time steps
    while (time<100000)
    {
        
        //Write
        if (hiperProbl->myRank() == 0)
        {
           if (testCaseMotorsTurnOff)
           {
            cout << "    Time step " + to_string(timeStep) + " "
                 << "deltat= " << userStr->dparam[0] << " and maximum iteration " << maxIter << " and Rotational power "<< userStr->dparam[9] <<" Isotropic power "<< userStr->dparam[7] <<" Anisotropic power "<< userStr->dparam[26] <<endl;
           }
           else
           {

            cout << "    Time step " + to_string(timeStep) + " " << "deltat= " << userStr->dparam[0] << " and maximum iteration " << maxIter << endl;

           }

        }
        //Initial guess
        dofHand->nodeDOFs->setValue(dofHand->nodeDOFs0);
        hiperProbl->UpdateGhosts();

        //Newton-Raphson method
        //Newton-Raphson method
        bool converged = nonlinSolver->solve();

        nIter = userNum->iparam[1];
        //Check convergence
        if (converged) //converged
        {

            //PostProc(dofHand,cortexHand,deltat);
            //Save solution
            dofHand->nodeDOFs0->setValue(dofHand->nodeDOFs);


            //Modify time-step size
            if (nIter <= 8)
                 deltat /= adaptiveStepTime;
            else if (nIter > 8 )
                deltat *= adaptiveStepTime;
            if (deltat > maxDelt)
                deltat = maxDelt;

            //Update time variables
            timeStep++;
            time = time + deltat;





            //Write results and update solution
            if (timeStep % nPrint == 0)
            {

/*                 //PostProc                
                 for (int i=0;i<disMesh->loc_nPts();i++)
                 {                     
                    double x = dofHand->mesh->nodeCoord(i, 0, IndexType::Local); 
                    double y = dofHand->mesh->nodeCoord(i, 1, IndexType::Local);
  
                    if (y>r_out-1e-3)
                    {
                        dofHand->nodeAuxF->setValue(1,i,IndexType::Local,hiperProbl->globIntegral("Tx_top"));
                        dofHand->nodeAuxF->setValue(2,i,IndexType::Local,hiperProbl->globIntegral("Ty_top"));
                    }
                    if (y<1e-3)
                    {
                        dofHand->nodeAuxF->setValue(1,i,IndexType::Local,hiperProbl->globIntegral("Tx_bottom"));
                        dofHand->nodeAuxF->setValue(2,i,IndexType::Local,hiperProbl->globIntegral("Ty_bottom"));
                    }
                    if (x<1e-3)
                    {
                        dofHand->nodeAuxF->setValue(1,i,IndexType::Local,hiperProbl->globIntegral("Tx_left"));
                        dofHand->nodeAuxF->setValue(2,i,IndexType::Local,hiperProbl->globIntegral("Ty_left"));
                    }
                    if (x>r_out-1e-2)
                    {
                        dofHand->nodeAuxF->setValue(1,i,IndexType::Local,hiperProbl->globIntegral("Tx_right"));
                        dofHand->nodeAuxF->setValue(2,i,IndexType::Local,hiperProbl->globIntegral("Ty_right"));
                    } 
            
                 } 

*/
                //storing total time for each converged time step.    
                dofHand->nodeAuxF->setValue(0,time);                

                solname_v = oname + "." + to_string(timeStep);
                if (printVtk)
                    dofHand->printFileLegacyVtk(solname_v, true);
                if (timeStep % 800 == 0)
                    dofHand->printFile(solname_v, OutputMode::Text, true, time);




                //     userStr->dparam[9] = 0.
                //     userStr->dparam[26] = 0.
                //        std::cout<<time<<","<<hiperProbl->globIntegral("Border")<<"i+"<<std::endl;
                if (dofHand->myRank()==0)    myfile<<hiperProbl->globIntegral("biLinear")<<" "<<time<<" "<<hiperProbl->globIntegral("sigma_xx")/(r_out)<<" "<<hiperProbl->globIntegral("sigma_yy")/(r_out)<<" "<<hiperProbl->globIntegral("sigma_partial_xx")/(r_out)<<" "<<hiperProbl->globIntegral("sigma_partial_yy")/(r_out)<<" "<<hiperProbl->globIntegral("gradV_xx")/(r_out*r_out)<<" "<<hiperProbl->globIntegral("gradV_yy")/(r_out*r_out)<<" "<<hiperProbl->globIntegral("Tx")-hiperProbl->globIntegral("fric_xx")<<" "<<hiperProbl->globIntegral("Ty")-hiperProbl->globIntegral("fric_yy")<<endl;
                //   
//                }

            }        

        }
        else
        {
            //Reduce time-step size
            deltat *= adaptiveStepTime;
            //maxDelt = deltat;
        }

        //Update deltat
                
        if (timeStep==1000)
        {


    for (int i = 0; i < dofHand->mesh->loc_nPts(); i++)
    {
 //        double x = dofHand->mesh->nodeCoord(i, 0, IndexType::Local);
 //        double y = dofHand->mesh->nodeCoord(i, 1, IndexType::Local);       
 //        int crease = dofHand->mesh->nodeCrease(i, IndexType::Local); 
           dofHand->nodeDOFs->setValue(3, i, IndexType::Local,  noiseStrength*fRand(-1,1.0)) ;
           dofHand->nodeDOFs->setValue(4, i, IndexType::Local,  noiseStrength*fRand(-1,1.0)) ;
           //   dofHand->setConstraint(0, i, IndexType::Local, 0.0);         
     }



        }    
        userStr->dparam[0] = deltat;

    }

    //Print last time step
    solname_v = oname + "." + to_string(timeStep+1);
    if (printVtk)
        dofHand->printFileLegacyVtk(solname_v, true);
    if (printFile)
        dofHand->printFile(solname_v, OutputMode::Text, true, time+deltat);

    //Finalize
    MPI_Finalize();
    return 0;

}


