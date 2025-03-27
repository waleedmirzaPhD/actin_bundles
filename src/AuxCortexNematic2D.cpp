
/// C++ headers
#include <iostream>
#include <mpi.h>
#include <fstream>

/// Trilinos headers
#include <Teuchos_RCP.hpp>

/// hiperlife headers
#include "hl_TypeDefs.h"
#include "hl_GlobalBasisFunctions.h"
#include "hl_StructMeshGenerator.h"
#include "hl_DistributedMesh.h"
#include "hl_FillStructure.h"
#include "hl_DOFsHandler.h"
#include "hl_HiPerProblem.h"
#include "hl_LinearSolver_Direct_Amesos.h"
#include "hl_Tensor.h"
#include "hl_LinearSolver_Iterative_AztecOO.h"
#include "AuxCortexNematic2D.h"
  


                                                                                                                                           
void setCircleVelBC_uniform(Teuchos::RCP<hiperlife::DOFsHandler> dofHand, Teuchos::RCP<hiperlife::UserStructure> userStr)
{

    using namespace hiperlife;
    double phi = userStr->dparam[15];
    double r_in = userStr->dparam[16];
    double r_out = userStr->dparam[17];
    double noise_strength = userStr->dparam[24];    
    double noise_wave_length = 20.0;
    double thickness = 0.125;
    double endAng   =  0.5*phi*180.0/M_PI +  thickness;
    double startAng =  0.5*phi*180.0/M_PI -  thickness;
    
    for (int i = 0; i < dofHand->mesh->loc_nPts(); i++) {
        double x = dofHand->mesh->nodeCoord(i, 0, IndexType::Local);
        double y = dofHand->mesh->nodeCoord(i, 1, IndexType::Local);
        // SETTING AUX-DOFs
        dofHand->nodeAuxF->setValue(0, i, IndexType::Local, x);
        dofHand->nodeAuxF->setValue(1, i, IndexType::Local, y);
        int crease = dofHand->mesh->nodeCrease(i, IndexType::Local);
        double rad = sqrt(x * x + y * y);

        if (crease > 0)
        {
            double theeta = atan2(y, x);
            // theeta_ is the angle between 0 and 2pi
            double theeta_;
            //Since tan2 returns angles in interval [-pi,pi]
            if (theeta<0)
                theeta_ = theeta +2*M_PI + noise_strength*fRand(-1.0,1.0);    
            else
                theeta_ = theeta + noise_strength*fRand(-1.0,1.0);                        
            
            double n1 = cos(theeta_);
            double n2 = sin(theeta_);
            double Vx, Vy, Vr;
            Vr = 0./60.0; // + noise_strength*fRand(-1.,1.);
            double slope = 30.0;
            //Components of velocity in cartesian coordiantes
            Vx = -Vr*n1;//*(0.5*(1+tanh(slope*(theeta_*180.0/M_PI-endAng))));
            Vy = -Vr*n2;//*(0.5*(1+tanh(slope*(theeta_*180.0/M_PI-endAng))));
            
        }
    }
}                                                                                                                                                                                                                                                                          






void setCircleThiBC_uniform(Teuchos::RCP<hiperlife::DOFsHandler> dofHand, Teuchos::RCP<hiperlife::UserStructure> userStr)
{
    using namespace hiperlife;
    double r_in = userStr->dparam[16];
    double r_out = userStr->dparam[17];
    double noise_strength  = userStr->dparam[24]; 
    double h_boundary=0.2;
    double noise_wave_length = 20.0;

    for (int i = 0; i < dofHand->mesh->loc_nPts(); i++)
    {
        double x = dofHand->mesh->nodeCoord(i, 0, IndexType::Local);
        double y = dofHand->mesh->nodeCoord(i, 1, IndexType::Local);       
        int crease = dofHand->mesh->nodeCrease(i, IndexType::Local);
        double rad = sqrt(x * x + y * y);
        double theeta = atan2(y, x);


       double theeta_{};
       //Since tan2 returns angles in interval [-pi,p
       if (theeta<0)
           theeta_ = theeta +2*M_PI;    
       else
           theeta_ = theeta;                        

       if (crease > 0)
       {
    
            // In case of annulus making sure to apply BC only on the outer boundary for the outer boundary
            if (rad > r_out - 0.001)
            { 
                 dofHand->nodeDOFs->setValue(0, i, IndexType::Local, h_boundary +  noise_strength*fRand(0,1.0));               
                 dofHand->setConstraint(0, i, IndexType::Local, 0.0);
            }
            else
            {
                 dofHand->nodeDOFs->setValue(0, i, IndexType::Local, h_boundary +  noise_strength*fRand(0,1.0));
         //        dofHand->setConstraint(0, i, IndexType::Local, 0.0);  
            } 
            
       }
       else 
       {
             dofHand->nodeDOFs->setValue(0, i, IndexType::Local, h_boundary +   noise_strength*fRand(0,1.0));
        //     dofHand->setConstraint(0, i, IndexType::Local, 0.0);  
       }
    }
}                                                                                                                                                                                       

               


void setCircleNemBC_square(Teuchos::RCP<hiperlife::DOFsHandler> dofHand, Teuchos::RCP<hiperlife::UserStructure> userStr)
{
    using namespace hiperlife;
    double phi = userStr->dparam[15];
    double r_in = userStr->dparam[16];
    double r_out = userStr->dparam[17];
    double noise_strength =   userStr->dparam[24];
    double sus20  = userStr->dparam[1];
    double sus40  = userStr->dparam[2];
    double heqb    = userStr->dparam[22];  
    double lambda_rot    = userStr->dparam[9];


    
    for (int i = 0; i < dofHand->mesh->loc_nPts(); i++)
    {
    
       double x = dofHand->mesh->nodeCoord(i, 0, IndexType::Local);
       double y = dofHand->mesh->nodeCoord(i, 1, IndexType::Local);
       double rad = sqrt(x*x + y*y);

       double theeta = atan2(y, x);
       double Q1{}, Q2{};
       double n1 = 1.0;
       double n2 = 0.0;
       
       double S_boundary = 0.3 + noise_strength*fRand(0,1.0);
       Q1 = S_boundary*((n1 * n1 - n2 * n2)) / (n1 * n1 + n2 * n2);
       Q2 = S_boundary*(2 * n1 * n2) / (n1 * n1 + n2 * n2);       

       int crease = dofHand->mesh->nodeCrease(i, IndexType::Local);
       double  q_int = -sqrt(-(4*sus20+lambda_rot*heqb)/(32*sus40) );
     
       dofHand->nodeDOFs->setValue(1, i, IndexType::Local, noise_strength*fRand(-1,1.0)); //-0.1768
       dofHand->nodeDOFs->setValue(2, i, IndexType::Local, noise_strength*fRand(-1,1.0));   
       //dofHand->setConstraint(2, i, IndexType::Local, 0.0);  
       //dofHand->setConstraint(4, i, IndexType::Local, 0.0);  

       cout<<q_int<<endl; //4*sus20+lambda_rot*heqb<<endl;

    
    }
}                                                                                                                                     


void setCircleThiBC_square(Teuchos::RCP<hiperlife::DOFsHandler> dofHand, Teuchos::RCP<hiperlife::UserStructure> userStr)
{
    using namespace hiperlife;
//    double r_in = userStr->dparam[16];
    double r_out = userStr->dparam[17];
    double noise_strength  = userStr->dparam[24]; 
    double heqb = userStr->dparam[22];
    double testCaseStrainRate = userStr->dparam[32];
    double strainRate = userStr->dparam[31];
    double visc    = userStr->dparam[10];
                                                                 
    for (int i = 0; i < dofHand->mesh->loc_nPts(); i++)
    {
//        double x = dofHand->mesh->nodeCoord(i, 0, IndexType::Local);
//        double y = dofHand->mesh->nodeCoord(i, 1, IndexType::Local);       
//        int crease = dofHand->mesh->nodeCrease(i, IndexType::Local);
        double h_init{};
        if (testCaseStrainRate>0.999)
        {
            h_init=visc*heqb/(visc + strainRate*1./r_out);          

        }
        else
        {
            h_init = heqb;
        } 

                    
        dofHand->nodeDOFs->setValue(0, i, IndexType::Local, h_init + noise_strength*fRand(0,1.0)) ;
          //   dofHand->setConstraint(0, i, IndexType::Local, 0.0);         
    }
} 
                                                                                                               

void setCircleVelBC_square(Teuchos::RCP<hiperlife::DOFsHandler> dofHand, Teuchos::RCP<hiperlife::UserStructure> userStr)
{

    using namespace hiperlife;
//    double noise_strength = userStr->dparam[24];
    double r_out = userStr->dparam[17];
    double strainRate = userStr->dparam[31];
    double testCaseStrainRate = userStr->dparam[32];

    if (testCaseStrainRate<0.9999)
        strainRate=0;

    for (int i = 0; i < dofHand->mesh->loc_nPts(); i++) 
    {
//        double x = dofHand->mesh->nodeCoord(i, 0, IndexType::Local);
        double y = dofHand->mesh->nodeCoord(i, 1, IndexType::Local);
        // SETTING AUX-DOFs
//        int crease = dofHand->mesh->nodeCrease(i, IndexType::Local);
//        double rad = sqrt(x * x + y * y);
        
        dofHand->nodeDOFs->setValue(3, i, IndexType::Local, 0.0) ;          
        dofHand->nodeDOFs->setValue(4, i, IndexType::Local, strainRate*(y-0.5*r_out)/r_out) ;
        if((testCaseStrainRate>0.9999) && (y<1e-3 or y>r_out-1e-3))
        {
            dofHand->setConstraint(4, i, IndexType::Local, 0.0);     
        }   
    }

} 



                                                                                                                               
void LS_CortexNematic_2D(Teuchos::RCP<hiperlife::FillStructure> fillStr)
{
    using namespace hiperlife;
    using namespace hiperlife::Tensor;
    //-----------------------------------------------------------
    //[1] INPUT DATa

    double deltat = fillStr->userStr->dparam[0];
    double hcrit  = fillStr->userStr->dparam[3];
    double lambda_iso    = fillStr->userStr->dparam[7];
    double lambda_aniso    = fillStr->userStr->dparam[8];
    double visc    = fillStr->userStr->dparam[10];
    double fric    = fillStr->userStr->dparam[13];
    double stab    = fillStr->userStr->dparam[14];
    double r_in    = fillStr->userStr->dparam[16];
    double r_out    = fillStr->userStr->dparam[17];
    double heqb    = fillStr->userStr->dparam[22];
    double sus20  = fillStr->userStr->dparam[1];
    double sus40  = fillStr->userStr->dparam[2];
    double rvisc   = fillStr->userStr->dparam[11];
    double cvisc   = fillStr->userStr->dparam[12];
    double frank  = fillStr->userStr->dparam[6];
    double lambda_rot    = fillStr->userStr->dparam[9];
    double crossLinker = fillStr->userStr->dparam[26];                     
    double c = 2.;
    double strainRate = fillStr->userStr->dparam[31];
    double testCaseStrainRate = fillStr->userStr->dparam[32];
    double testCase = fillStr->userStr->dparam[30];
    // Scale the free energy parameters with time
    frank = frank/deltat; 
    sus20 = sus20/deltat;
    sus40 = sus40/deltat;                                       

    auto& subFill = (*fillStr)["cortexHand"];
    int nDim = subFill.nDim;
    int pDim = subFill.pDim;
    int eNN  = subFill.eNN;
    int numDOFs = subFill.numDOFs;
    int N_e  = numDOFs*eNN;
    double *bf  = subFill.getDer(0);
    double *Dbf_l = subFill.getDer(1);
    vector<double>& nborCoords  = subFill.nborCoords;
    vector<double>& nborDOFs  = subFill.nborDOFs;
    vector<double>& nborDOFs0 = subFill.nborDOFs0;               

    //[2] GEOMETRY
    double TMap[4], iTMap[4];
    Array::Fill(TMap, 4, 0.0);
    for (int i = 0; i < eNN; i++)
    {
        TMap[0] += nborCoords[nDim * i + 0] * Dbf_l[pDim * i + 0];
        TMap[1] += nborCoords[nDim * i + 0] * Dbf_l[pDim * i + 1];
        TMap[2] += nborCoords[nDim * i + 1] * Dbf_l[pDim * i + 0];
        TMap[3] += nborCoords[nDim * i + 1] * Dbf_l[pDim * i + 1];
    }

    Math::Invert2x2(iTMap, TMap);
    double jac= Math::DetMat2x2(TMap);


    //[3] VARIABLES
    double xCoords[3], fricTensor_along[3]{}, fricTensor_trans[3]{};
    double v[2],v0[2],  Dv[4],Dv0[4],  DDv[8], rodt[4], rodt0[4], omega[4], divv, divv0;
    double h{}, h0{}, Dh[2], Dh0[2];
    double q[2], q0[2], Dq[4], Dq0[4];


    Array::Fill(xCoords, 3, 0.0);
    Array::Fill(v, 2, 0.0);
    Array::Fill(v0, 2, 0.0);
    Array::Fill(Dv, 4, 0.0);
    Array::Fill(Dv0, 4, 0.0);
    Array::Fill(q, 2, 0.0);
    Array::Fill(q0, 2, 0.0);
    Array::Fill(Dq, 4, 0.0);
    Array::Fill(Dq0, 4, 0.0);
    Array::Fill(Dh, 2, 0.0);
    Array::Fill(Dh0, 2, 0.0);                                                                                          
                                           
    for (int i = 0; i < eNN; i++)
    {
        Math::AXPY(xCoords, nDim, bf[i], &nborCoords[nDim*i]);

        //Height
        h  += bf[i] * nborDOFs[numDOFs*i];
        h0 += bf[i] * nborDOFs0[numDOFs*i];

        //Velocities
        Math::AXPY(v, 2, bf[i], &nborDOFs[numDOFs * i+3]);
        Math::AXPY(v0, 2, bf[i], &nborDOFs0[numDOFs * i+3]);

        //Gradients
        double *Dbf_lI = &Dbf_l[2*i];
        double  Dbf_gI[2];
        Math::MatProduct(Dbf_gI, 1, 2, 2, Dbf_lI, iTMap);

        //Height gradients
        Dh[0]  += Dbf_gI[0]  * nborDOFs[numDOFs * i ];
        Dh[1]  += Dbf_gI[1]  * nborDOFs[numDOFs * i ];

        Dh0[0]  += Dbf_gI[0]  * nborDOFs0[numDOFs * i ];
        Dh0[1]  += Dbf_gI[1]  * nborDOFs0[numDOFs * i ];



        //Velocity gradients
        Dv[0]  += Dbf_gI[0] * nborDOFs[numDOFs * i + 3];
        Dv[1]  += Dbf_gI[1] * nborDOFs[numDOFs * i + 3];
        Dv[2]  += Dbf_gI[0] * nborDOFs[numDOFs * i + 4];
        Dv[3]  += Dbf_gI[1] * nborDOFs[numDOFs * i + 4];


        Dv0[0]  += Dbf_gI[0] * nborDOFs0[numDOFs * i + 3];
        Dv0[1]  += Dbf_gI[1] * nborDOFs0[numDOFs * i + 3];
        Dv0[2]  += Dbf_gI[0] * nborDOFs0[numDOFs * i + 4];
        Dv0[3]  += Dbf_gI[1] * nborDOFs0[numDOFs * i + 4];



        //Nematic
        Math::AXPY(q, 2, bf[i], &nborDOFs[numDOFs * i +1]);
        Math::AXPY(q0, 2, bf[i], &nborDOFs0[numDOFs * i +1 ]);


        //Nematic gradients
        Dq[0]  += Dbf_gI[0] * nborDOFs[numDOFs * i + 1];
        Dq[1]  += Dbf_gI[1] * nborDOFs[numDOFs * i + 1];
        Dq[2]  += Dbf_gI[0] * nborDOFs[numDOFs * i + 2];
        Dq[3]  += Dbf_gI[1] * nborDOFs[numDOFs * i + 2];

        Dq0[0] += Dbf_gI[0] * nborDOFs0[numDOFs * i + 1];
        Dq0[1] += Dbf_gI[1] * nborDOFs0[numDOFs * i + 1];
        Dq0[2] += Dbf_gI[0] * nborDOFs0[numDOFs * i + 2];
        Dq0[3] += Dbf_gI[1] * nborDOFs0[numDOFs * i + 2];

    }

    //Divergence of velocity
    divv = Dv[0] + Dv[3];
    divv0 = Dv0[0] + Dv0[3];
    // Rate of deformation tensor
    rodt[0] = Dv[0];
    rodt[1] = 0.5*(Dv[1]  + Dv[2]) ; 
    rodt[2] = 0.5*(Dv[2]  + Dv[1]) ;
    rodt[3] = Dv[3] ;                                

    rodt0[0] = Dv0[0];
    rodt0[1] = 0.5*(Dv0[1]  + Dv0[2]); 
    rodt0[2] = 0.5*(Dv0[2]  + Dv0[1]);
    rodt0[3] = Dv0[3];                              

    // Vorticity tensor
    omega[0] = 0.0;
    omega[1] = 0.5*(Dv[1]  - Dv[2]);
    omega[2] = 0.5*(Dv[2]  - Dv[1]);
    omega[3] = 0.0;

    // Jaumann derivative
    double Jq[2]{};
    Jq[0] = (q[0]-q0[0])/deltat + v[0]*Dq[0] + v[1]*Dq[1] - 2*omega[1]*q[1];
    Jq[1] = (q[1]-q0[1])/deltat + v[0]*Dq[2] + v[1]*Dq[3] + 2*omega[1]*q[0];

    double S = sqrt(4*(q0[0]*q0[0] + q0[1]*q0[1]));
    double S2  = 4*(q[0]*q[0] + q[1]*q[1] );
    double S20 = 4*(q0[0]*q0[0] + q0[1]*q0[1]);               

    // Check to make sure system stays contractile

    double Snn0[3], Snn0_[3];
    Snn0[0]  =  q0[0] + S/2.;
    Snn0[1]  =  q0[1];
    Snn0[2]  =  -q0[0] + S/2  ;
    // n_ is transverse to the filament orientations n
    Snn0_[0] =  Snn0[2];
    Snn0_[1] =  -Snn0[1];    
    Snn0_[2] = Snn0[0];                  

//    double h_ =  h0 + deltat*( -v[0]*Dh[0] -v[1]*Dh[1] - h*divv+ visc*(heqb-h) );
    // Terms used in frank constant
    double DqDq = 0.5*(Dq[0]*Dq[0] + Dq[1]*Dq[1] + Dq[2]*Dq[2] + Dq[3]*Dq[3]);


//    double rad = sqrt(xCoords[0]*xCoords[0] + xCoords[1]*xCoords[1]); 
//    double n1 = cos(theeta_);
//    double n2 = sin(theeta_);
    double Vx{}, Vy{};
    if (testCaseStrainRate>0.999)
    {       

        double h_init=visc*heqb/(visc + strainRate*1./r_out);       
        lambda_rot = lambda_rot*heqb/h_init;                       
        Vx = 0.;
        Vy =  strainRate*(xCoords[1]-0.5*r_out)/r_out;
    }                                                                                 


                                                                                                                                                                                                                                              
   //[3] OUTPUT PARAMETERS
    double rayleighian{};                                                                                                                                                                                                                                                                                     
    for (int i = 0; i < eNN; i++)
    {

       //Gradients
        double *Dbf_lI = &Dbf_l[2*i];
        double  Dbf_gI[2];
        Math::MatProduct(Dbf_gI, 1, 2, 2, Dbf_lI, iTMap);

        // SHEAR DISSIPATION
        double  dv1rodtI[4], dv2rodtI[4];

        dv1rodtI[0] =  Dbf_gI[0];
        dv1rodtI[1] =  0.5*Dbf_gI[1];
        dv1rodtI[2] =  0.5*Dbf_gI[1];
        dv1rodtI[3] =  0.0;

        dv2rodtI[0] = 0.0;
        dv2rodtI[1] = 0.5*Dbf_gI[0];
        dv2rodtI[2] = 0.5*Dbf_gI[0];
        dv2rodtI[3] = Dbf_gI[1];

        // CONTINUITY EQUATION
        double dhh_I       =  bf[i];
        double dv1h_I      =  0.0;
        double dv2h_I      =  0.0;


       // ROTATIONAL VISCOSITY
        double dq1dJq_I[2], dq2dJq_I[2];
        dq1dJq_I[0] = bf[i]/deltat;
        dq1dJq_I[1] = 0.0;
        dq2dJq_I[0] = 0.0;
        dq2dJq_I[1] = bf[i]/deltat;

        double dv1dJq_I[2], dv2dJq_I[2];

        dv1dJq_I[0] = bf[i]*Dq[0] - Dbf_gI[1]*q[1];
        dv1dJq_I[1] = bf[i]*Dq[2] + Dbf_gI[1]*q[0];

        dv2dJq_I[0] = bf[i]*Dq[1] + Dbf_gI[0]*q[1];
        dv2dJq_I[1] = bf[i]*Dq[3] - Dbf_gI[0]*q[0];

        //FREE ENERGY (SUSCEPTBILITY)
        double dq1dsus2_I, dq2dsus2_I;
        dq1dsus2_I =   4*bf[i]*q[0] ;
        dq2dsus2_I =   4*bf[i]*q[1] ;

        double dq1dsus4_I, dq2dsus4_I;
        dq1dsus4_I =  8*(q[0]*q[0] + q[1]*q[1] )*2*(bf[i]*q[0]);
        dq2dsus4_I =  8*(q[0]*q[0] + q[1]*q[1] )*2*(bf[i]*q[1]);

        //FREE ENERGY (FRANK)
        double dq1DqDq_I, dq2DqDq_I;

        dq1DqDq_I = Dbf_gI[0]*Dq[0] +  Dbf_gI[1]*Dq[1]; 
        dq2DqDq_I = Dbf_gI[0]*Dq[2] +  Dbf_gI[1]*Dq[3];

    
        double dv1dh_I,  dv2dh_I, dhdh_I;

        dv1dh_I = deltat*( -bf[i]*Dh[0]-h*Dbf_gI[0] ); 
        dv2dh_I = deltat*( -bf[i]*Dh[1]-h*Dbf_gI[1] );
        dhdh_I  = deltat*( -v[0]*Dbf_gI[0] -v[1]*Dbf_gI[1] -bf[i]*divv - visc*bf[i]  );
       

        double dv1divvI = Dbf_gI[0] ;
        double dv2divvI = Dbf_gI[1] ;





        // [9] SHEAR DISSIPATION
        // [9.2] RESIDUAL
        fillStr->Bk(0)[i*numDOFs + 3] += c*(jac *h) * (rodt[0] * dv1rodtI[0] +  rodt[1] * dv1rodtI[1] + rodt[2] * dv1rodtI[2] + rodt[3] * dv1rodtI[3]);
        fillStr->Bk(0)[i*numDOFs + 4] += c*(jac *h) * (rodt[0] * dv2rodtI[0] +  rodt[1] * dv2rodtI[1] + rodt[2] * dv2rodtI[2]+  rodt[3] * dv2rodtI[3]);

        fillStr->Bk(0)[i*numDOFs + 3] +=c*(jac *h) * (divv*Dbf_gI[0]) ;
        fillStr->Bk(0)[i*numDOFs + 4] +=c*(jac *h) * (divv*Dbf_gI[1]) ;




        //[10-11-12] CONTINUITY EQUATIONS 
        fillStr->Bk(0)[i*numDOFs + 0] += jac*( bf[i]*((h-h0)/deltat + Dh[0]*v[0] + Dh[1]*v[1] + h*divv -visc*(heqb-h)) +  (Dbf_gI[0]*Dh[0]  + Dbf_gI[1]*Dh[1])); 

        //[17]  ISOTROPIC MYOSIN TERM
        fillStr->Bk(0)[i*numDOFs + 3] += (jac * h * lambda_iso)*(1.0 - S) *Dbf_gI[0];
        fillStr->Bk(0)[i*numDOFs + 4] += (jac * h * lambda_iso)*(1.0 - S) *Dbf_gI[1];


        // FRICTION TERM
        fillStr->Bk(0)[i*numDOFs + 3] += (jac *h) *(v[0]-Vx)*bf[i];
        fillStr->Bk(0)[i*numDOFs + 4] += (jac *h) *(v[1]-Vy)*bf[i];

///      fillStr->Bk(0)[i*numDOFs + 3] += 0.5*(jac *h0) *(v[0]*bf[i]);
//        fillStr->Bk(0)[i*numDOFs + 4] += 0.5*(jac *h0) *((v[1]-0.0)*bf[i]);


       //[13] ROTATIONAL DISSIPATION
        fillStr->Bk(0)[i*numDOFs + 1] += (0.5 * jac * h * rvisc) *4*(dq1dJq_I[0]*Jq[0] + dq1dJq_I[1]*Jq[1]);
        fillStr->Bk(0)[i*numDOFs + 2] += (0.5 * jac * h * rvisc) *4*(dq2dJq_I[0]*Jq[0] + dq2dJq_I[1]*Jq[1]);
        fillStr->Bk(0)[i*numDOFs + 3] += (0.5 * jac * h * rvisc) *4*(dv1dJq_I[0]*Jq[0] + dv1dJq_I[1]*Jq[1]);
        fillStr->Bk(0)[i*numDOFs + 4] += (0.5 * jac * h * rvisc) *4*(dv2dJq_I[0]*Jq[0] + dv2dJq_I[1]*Jq[1]);


        //[14] COUPLED DISSIPATION
        fillStr->Bk(0)[i*numDOFs + 1] += (jac * h * cvisc) *(dq1dJq_I[0]*rodt[0] + 2*dq1dJq_I[1]*rodt[1]  -  dq1dJq_I[0]*rodt[3] );
        fillStr->Bk(0)[i*numDOFs + 2] += (jac * h * cvisc) *(dq2dJq_I[0]*rodt[0] + 2*dq2dJq_I[1]*rodt[1]  -  dq2dJq_I[0]*rodt[3] );
        fillStr->Bk(0)[i*numDOFs + 3] += (jac * h * cvisc) *((dv1dJq_I[0]*rodt[0] + 2*dv1dJq_I[1]*rodt[1] -  dv1dJq_I[0]*rodt[3] ) + (Jq[0]*dv1rodtI[0] + 2*Jq[1]*dv1rodtI[1] -  Jq[0]*dv1rodtI[3] ));
        fillStr->Bk(0)[i*numDOFs + 4] += (jac * h * cvisc) *((dv2dJq_I[0]*rodt[0] + 2*dv2dJq_I[1]*rodt[1] -  dv2dJq_I[0]*rodt[3] ) + (Jq[0]*dv2rodtI[0] + 2*Jq[1]*dv2rodtI[1] -  Jq[0]*dv2rodtI[3] ));


        //[15] FRANK CONSTANT
        fillStr->Bk(0)[i*numDOFs + 1] += (jac * h*frank) *(dq1DqDq_I);
        fillStr->Bk(0)[i*numDOFs + 2] += (jac * h*frank) *(dq2DqDq_I);
        fillStr->Bk(0)[i*numDOFs + 3] += (jac * dv1dh_I*frank) *DqDq;
        fillStr->Bk(0)[i*numDOFs + 4] += (jac * dv2dh_I*frank) *DqDq;         


        // POWER TERMS
        fillStr->Bk(0)[i*numDOFs + 3] += (jac * h*crossLinker) *(Dbf_gI[0]*Snn0_[0]  + Dbf_gI[1]*Snn0_[1]) ;
        fillStr->Bk(0)[i*numDOFs + 4] += (jac * h*crossLinker) *(Dbf_gI[1]*Snn0_[2]  + Snn0_[1]*Dbf_gI[0]);

        fillStr->Bk(0)[i*numDOFs + 3] += (jac * h*lambda_aniso) *(Dbf_gI[0]*Snn0[0]  + Dbf_gI[1]*Snn0[1]); 
        fillStr->Bk(0)[i*numDOFs + 4] += (jac * h*lambda_aniso) *(Dbf_gI[1]*Snn0[2]  + Snn0[1]*Dbf_gI[0]);                                                                                                                    

        for (int j = 0; j < eNN; j++)
        {
            // Indices of the Hessian tensor
             int row = i*numDOFs;
             int col = j*numDOFs;
            // Local and global gradients in terms of indices
            //Gradients
             double *Dbf_lJ = &Dbf_l[2*j];
             double  Dbf_gJ[2];
             Math::MatProduct(Dbf_gJ, 1, 2, 2, Dbf_lJ, iTMap);
    
        
            // SHEAR DISSIPATION
            double  dv1rodt_J[4], dv2rodt_J[4];
            dv1rodt_J[0] =  Dbf_gJ[0];
            dv1rodt_J[1] =  0.5*Dbf_gJ[1];
            dv1rodt_J[2] =  0.5*Dbf_gJ[1];
            dv1rodt_J[3] =  0.0;

            dv2rodt_J[0] = 0.0;
            dv2rodt_J[1] = 0.5*Dbf_gJ[0];
            dv2rodt_J[2] = 0.5*Dbf_gJ[0];
            dv2rodt_J[3] = Dbf_gJ[1];


            // CONSERVATION OF MASS

            double dv1h_J      =  0.0;
            double dv2h_J      =  0.0;
            double dhh_J       =  bf[j]; // + L*(DDbf_gJ[0] + DDbf_gJ[3])

//h0 + deltat*( -v[0]*Dh[0] -v[1]*Dh[1] + visc*(heqb-h) );

            double dv1dh_J = deltat*( -bf[j]*Dh[0]-h*Dbf_gJ[0] ); 
            double dv2dh_J = deltat*( -bf[j]*Dh[1]-h*Dbf_gJ[1] );
            double dhdh_J  = deltat*( -v[0]*Dbf_gJ[0] -v[1]*Dbf_gJ[1] -bf[j]*divv - visc*bf[j]  );
            

            double dhdv1dh_IJ  = deltat*( -bf[i]*Dbf_gJ[0] - bf[j]*Dbf_gI[0]); 
            double dhdv2dh_IJ  = deltat*( -bf[i]*Dbf_gJ[1] - bf[j]*Dbf_gI[1]);
            double dv1dhdh_IJ  = deltat*( -bf[j]*Dbf_gI[0] - bf[i]*Dbf_gJ[0]);
            double dv2dhdh_IJ  = deltat*( -bf[j]*Dbf_gI[1] - bf[i]*Dbf_gJ[1]);

            double dq1dsus4_J, dq2dsus4_J;
            dq1dsus4_J =  8*(q[0]*q[0] + q[1]*q[1] )*2*(bf[j]*q[0]);
            dq2dsus4_J =  8*(q[0]*q[0] + q[1]*q[1] )*2*(bf[j]*q[1]);



            double dq1dsus2_J, dq2dsus2_J, dq1dq1dsus2_IJ, dq2dq1dsus2_IJ, dq1dq2dsus2_IJ, dq2dq2dsus2_IJ ;
            dq1dsus2_J =   4*bf[j]*q[0] ;
            dq2dsus2_J =   4*bf[j]*q[1] ;
            dq1dq1dsus2_IJ =   4.0*bf[j]*bf[i];
            dq1dq2dsus2_IJ =   0.0;
            dq2dq1dsus2_IJ =   0.0;
            dq2dq2dsus2_IJ =   4.0*bf[j]*bf[i];                                                                    

            double dq1dq1dsus4_IJ, dq1dq2dsus4_IJ,dq2dq1dsus4_IJ, dq2dq2dsus4_IJ;;
            dq1dq1dsus4_IJ =  8*(q[0]*q[0] + q[1]*q[1] )*2*(bf[i]*bf[j]) +   8*(2*bf[j]*q[0])*2*(bf[i]*q[0]);
            dq1dq2dsus4_IJ =  8*(2.*bf[j]*q[1])*2*(bf[i]*q[0]);    
            dq2dq1dsus4_IJ =  8*(2.*bf[j]*q[0])*2*(bf[i]*q[1]);
            dq2dq2dsus4_IJ =  8*(q[0]*q[0] + q[1]*q[1] )*2*(bf[i]*bf[j]) +   8*(2*bf[j]*q[1])*2*(bf[i]*q[1]);


            // FREE ENERGY (FRANK) 
            double dq1DqDq_J, dq2DqDq_J;
            double dq1dq1DqDq_IJ, dq1dq2DqDq_IJ, dq2dq1DqDq_IJ, dq2dq2DqDq_IJ;
    
            dq1DqDq_J = (Dbf_gJ[0]*Dq[0] +  Dbf_gJ[1]*Dq[1]); ;
            dq2DqDq_J = (Dbf_gJ[0]*Dq[2] +  Dbf_gJ[1]*Dq[3]);

            dq1dq1DqDq_IJ = (Dbf_gI[0]*Dbf_gJ[0] +  Dbf_gI[1]*Dbf_gJ[1]); 
            dq2dq1DqDq_IJ = 0.0;
            dq1dq2DqDq_IJ = 0.0;
            dq2dq2DqDq_IJ = (Dbf_gI[0]*Dbf_gJ[0] +  Dbf_gI[1]*Dbf_gJ[1]);

                                                         

            double dq1dJq_J[2], dq2dJq_J[2];
            dq1dJq_J[0] = bf[j]/deltat + v[0]*Dbf_gJ[0] + v[1]*Dbf_gJ[1];
            dq1dJq_J[1] =  2*omega[1]*bf[j];
            dq2dJq_J[0] =  -2*omega[1]*bf[j];
            dq2dJq_J[1] = bf[j]/deltat+v[0]*Dbf_gJ[0] + v[1]*Dbf_gJ[1];

            double dv1dJq_J[2], dv2dJq_J[2];

            dv1dJq_J[0] = bf[j]*Dq[0] - Dbf_gJ[1]*q[1];
            dv1dJq_J[1] = bf[j]*Dq[2] + Dbf_gJ[1]*q[0];

            dv2dJq_J[0] = bf[j]*Dq[1] + Dbf_gJ[0]*q[1];
            dv2dJq_J[1] = bf[j]*Dq[3] - Dbf_gJ[0]*q[0];


            double dq1dv1dJq_IJ[2], dq1dv2dJq_IJ[2], dq2dv1dJq_IJ[2], dq2dv2dJq_IJ[2];

            dq1dv1dJq_IJ[0] = bf[i]*Dbf_gJ[0];
            dq1dv1dJq_IJ[1] = Dbf_gI[1]*bf[j];

            dq1dv2dJq_IJ[0] = bf[i]*Dbf_gJ[1];
            dq1dv2dJq_IJ[1] = - Dbf_gI[0]*bf[j];


            dq2dv1dJq_IJ[0] = - Dbf_gI[1]*bf[j];
            dq2dv1dJq_IJ[1] = bf[i]*Dbf_gJ[0];

            dq2dv2dJq_IJ[0] = Dbf_gI[0]*bf[j];
            dq2dv2dJq_IJ[1] = bf[i]*Dbf_gJ[1]; 


            double dv1divvJ = Dbf_gJ[0] ;
            double dv2divvJ = Dbf_gJ[1] ;


  //          double dbfmod_dv1_IJ =        stab*(Dbf_gI[0]*bf[j] +  Dbf_gJ[0]*bf[i]  );
  //          double dbfmod_dv2_IJ =        stab*(Dbf_gI[1]*bf[j] +  Dbf_gJ[1]*bf[i]  );

    


            // [9] SHEAR DISSIPATION
            // [9.2] HESSIAN
            fillStr->Ak(0, 0)[(row + 3) * N_e + col + 3] += c*(jac   *(h)) * (dv1rodt_J[0] * dv1rodtI[0] +  dv1rodt_J[1] * dv1rodtI[1] + dv1rodt_J[2] * dv1rodtI[2] + dv1rodt_J[3] * dv1rodtI[3]);
            fillStr->Ak(0, 0)[(row + 3) * N_e + col + 4] += c*(jac   *(h)) * (dv2rodt_J[0] * dv1rodtI[0] +  dv2rodt_J[1] * dv1rodtI[1] + dv2rodt_J[2] * dv1rodtI[2] + dv2rodt_J[3] * dv1rodtI[3]);
            fillStr->Ak(0, 0)[(row + 3) * N_e + col + 0] += c*(jac *dhh_J) * (rodt[0] * dv1rodtI[0] +  rodt[1] * dv1rodtI[1] + rodt[2] * dv1rodtI[2] + rodt[3] * dv1rodtI[3]);
            fillStr->Ak(0, 0)[(row + 3) * N_e + col + 3] += c*(jac *dv1h_J) *(rodt[0] * dv1rodtI[0] +  rodt[1] * dv1rodtI[1] + rodt[2] * dv1rodtI[2] + rodt[3] * dv1rodtI[3]); ;
            fillStr->Ak(0, 0)[(row + 3) * N_e + col + 4] += c*(jac *dv2h_J) * (rodt[0] * dv1rodtI[0] +  rodt[1] * dv1rodtI[1] + rodt[2] * dv1rodtI[2] + rodt[3] * dv1rodtI[3]);


            fillStr->Ak(0, 0)[(row + 4) * N_e + col + 3] += c*(jac   *(h)) * (dv1rodt_J[0] * dv2rodtI[0] +  dv1rodt_J[1] * dv2rodtI[1]+  dv1rodt_J[2] * dv2rodtI[2] + dv2rodt_J[3] * dv1rodtI[3]);
            fillStr->Ak(0, 0)[(row + 4) * N_e + col + 4] +=c* (jac   *(h)) * (dv2rodt_J[0] * dv2rodtI[0] +  dv2rodt_J[1] * dv2rodtI[1] + dv2rodt_J[2] * dv2rodtI[2] + dv2rodt_J[3] * dv2rodtI[3]);
            fillStr->Ak(0, 0)[(row + 4) * N_e + col + 0] += c*(jac *dhh_J) *(rodt[0] * dv2rodtI[0] +  rodt[1] * dv2rodtI[1] + rodt[2] * dv2rodtI[2]+  rodt[3] * dv2rodtI[3]);
            fillStr->Ak(0, 0)[(row + 4) * N_e + col + 3] += c*(jac *dv1h_J)*(rodt[0] * dv2rodtI[0] +  rodt[1] * dv2rodtI[1] + rodt[2] * dv2rodtI[2]+  rodt[3] * dv2rodtI[3]);
            fillStr->Ak(0, 0)[(row + 4) * N_e + col + 4] += c*(jac *dv2h_J)*(rodt[0] * dv2rodtI[0] +  rodt[1] * dv2rodtI[1] + rodt[2] * dv2rodtI[2]+  rodt[3] * dv2rodtI[3]);
                                                                                                                                                                                                                                                      

            fillStr->Ak(0, 0)[(row + 3) * N_e + col + 0] += c*(jac *bf[j]) * (divv*dv1divvI) ;
            fillStr->Ak(0, 0)[(row + 3) * N_e + col + 3] += c*(jac *h) * (dv1divvJ*dv1divvI) ;
            fillStr->Ak(0, 0)[(row + 3) * N_e + col + 4] += c*(jac *h) * (dv2divvJ*dv1divvI) ;

            fillStr->Ak(0, 0)[(row + 4) * N_e + col + 0] += c*(jac *bf[j]) * (divv*dv2divvI); 
            fillStr->Ak(0, 0)[(row + 4) * N_e + col + 3] += c*(jac *h)* (dv2divvI*dv1divvJ) ;
            fillStr->Ak(0, 0)[(row + 4) * N_e + col + 4] += c*(jac*h) * (dv2divvJ*dv2divvI)  ;



           // bf_mod_v_I       = bf[i] + stab*ledge*(Dbf_gI[0]*v[0]+ Dbf_gI[1]*v[1]);     
           
            fillStr->Ak(0, 0)[(row + 0) * N_e + col + 0] +=  jac*(bf[i]*  ( bf[j]/deltat + Dbf_gJ[0]*v[0] + Dbf_gJ[1]*v[1] + bf[j]*divv + visc*bf[j]) +  (Dbf_gI[0]*Dbf_gJ[0]  + Dbf_gI[1]*Dbf_gJ[1]))  ;//jac *( bf_mod_v_I*(bf[j] - dhh_J)  +  tol*deltat*(Dbf_gI[0]*Dbf_gJ[0]+Dbf_gI[1]*Dbf_gJ[1])); 
                                                            
            fillStr->Ak(0, 0)[(row + 0) * N_e + col + 3] +=  jac*bf[i]*( Dh[0]*bf[j]  + h*Dbf_gJ[0]) ; // + stab*ledge*( 2*v[0]*bf[j] )*(Dbf_gI[0]*Dh[0]  + Dbf_gI[1]*Dh[1] ) +   stab*ledge*divv*h*bf[j]*Dbf_gI[0] +  
                                                                    
            fillStr->Ak(0, 0)[(row + 0) * N_e + col + 4] +=  jac*bf[i]*( Dh[1]*bf[j]  + h*Dbf_gJ[1] );// +º  stab*ledge*( 2*v[1]*bf[j] )*(Dbf_gI[0]*Dh[0]  + Dbf_gI[1]*Dh[1] ) +   stab*ledge*divv*h*bf[j]*Dbf_gI[1] + stab*ledge*Dbf_gJ[1]*h*( v[0]*Dbf_gI[0] + v[1]*Dbf_gI[1])) ;
//            fillStr->Ak(0, 0)[(row + 0) * N_e + col + 0] +=  jac *( bf[i]*(bf[j] - dhh_J)  + tol*deltat*(Dbf_gI[0]*Dbf_gJ[0]+Dbf_gI[1]*Dbf_gJ[1]));
 //           fillStr->Ak(0, 0)[(row + 0) * N_e + col + 0] +=  0.5*jac*(bf[i]*  ( bf[j]/deltat));                                                                                                                               



            //FRICTION

            fillStr->Ak(0, 0)[(row + 3) * N_e + col + 3] +=  (jac *h) *(bf[j]*bf[i]) ;
            fillStr->Ak(0, 0)[(row + 3) * N_e + col + 0] +=  (jac *dhh_J) *(bf[i]*(v[0]-Vx)) ;

            fillStr->Ak(0, 0)[(row + 4) * N_e + col + 4] +=  (jac *h) *(bf[j]*bf[i]) ;
            fillStr->Ak(0, 0)[(row + 4) * N_e + col + 0] +=  (jac *dhh_J) *(bf[i]*(v[1]-Vy)) ;

            // SUSCEPTIBILITY 
//        fillStr->Bk(0)[i*numDOFs + 3] += (jac*dv1dh_I*sus20)*(2*(q[0]*q[0] + q[1]*q[1] ));   
/*
            fillStr->Ak(0, 0)[(row + 1) * N_e + col + 0] +=  (4.*jac* dhdh_J*sus20)*bf[i]*q[0];  
            fillStr->Ak(0, 0)[(row + 1) * N_e + col + 1] +=  (4.*jac*h_*sus20)*bf[i]*bf[j];    
            fillStr->Ak(0, 0)[(row + 1) * N_e + col + 2] +=  0.;  
            fillStr->Ak(0, 0)[(row + 1) * N_e + col + 3] +=  (4.*jac*dv1dh_J*sus20)*bf[i]*q[0];  
            fillStr->Ak(0, 0)[(row + 1) * N_e + col + 4] +=  (4.*jac*dv2dh_J*sus20)*bf[i]*q[0];  


            fillStr->Ak(0, 0)[(row + 2) * N_e + col + 0] +=  (4.*jac* dhdh_J*sus20)*bf[i]*q[1];  
            fillStr->Ak(0, 0)[(row + 2) * N_e + col + 1] +=  0.0;    
            fillStr->Ak(0, 0)[(row + 2) * N_e + col + 2] +=  (4.*jac*h_*sus20)*bf[i]*bf[j];       
            fillStr->Ak(0, 0)[(row + 2) * N_e + col + 3] +=  (4.*jac*dv1dh_J*sus20)*bf[i]*q[1];  
            fillStr->Ak(0, 0)[(row + 2) * N_e + col + 4] +=  (4.*jac*dv2dh_J*sus20)*bf[i]*q[1];  


            fillStr->Ak(0, 0)[(row + 3) * N_e + col + 0] +=  (jac*dhdv1dh_IJ*sus20)*(2*(q[0]*q[0] + q[1]*q[1] )); 
            fillStr->Ak(0, 0)[(row + 3) * N_e + col + 3] +=  0.0; 
            fillStr->Ak(0, 0)[(row + 3) * N_e + col + 4] +=  0.0;  
            fillStr->Ak(0, 0)[(row + 3) * N_e + col + 1] +=  (jac*dv1dh_I*sus20)*(4*(bf[j]*q[0]));  
            fillStr->Ak(0, 0)[(row + 3) * N_e + col + 2] +=  (jac*dv1dh_I*sus20)*(4*(bf[j]*q[1]));  


            fillStr->Ak(0, 0)[(row + 4) * N_e + col + 0] +=   (jac*dhdv2dh_IJ*sus20)*(2*(q[0]*q[0] + q[1]*q[1] ));   
            fillStr->Ak(0, 0)[(row + 4) * N_e + col + 3] +=   0.0;      
            fillStr->Ak(0, 0)[(row + 4) * N_e + col + 4] +=   0.0;   
            fillStr->Ak(0, 0)[(row + 4) * N_e + col + 1] +=   (jac*dv2dh_I*sus20)*(4*(bf[j]*q[0])); 
            fillStr->Ak(0, 0)[(row + 4) * N_e + col + 2] +=   (jac*dv2dh_I*sus20)*(4*(bf[j]*q[1]));                                           

*/


/*            fillStr->Ak(0, 0)[(row + 1) * N_e + col + 0] +=  (jac*dhdh_J*sus40)*dq1dsus4_I ; ;  
            fillStr->Ak(0, 0)[(row + 1) * N_e + col + 1] +=  (jac*h_*sus40)*dq1dq1dsus4_IJ ;  
            fillStr->Ak(0, 0)[(row + 1) * N_e + col + 2] +=  (jac*h_*sus40)*dq2dq1dsus4_IJ ;   
            fillStr->Ak(0, 0)[(row + 1) * N_e + col + 3] +=  (jac*dv1dh_J*sus40)*dq1dsus4_I ;  
            fillStr->Ak(0, 0)[(row + 1) * N_e + col + 4] +=  (jac*dv2dh_J*sus40)*dq1dsus4_I ;   


            fillStr->Ak(0, 0)[(row + 2) * N_e + col + 0] +=  (jac*dhdh_J*sus40)*dq2dsus4_I ;   
            fillStr->Ak(0, 0)[(row + 2) * N_e + col + 1] +=  (jac*h_*sus40)*dq2dq1dsus4_IJ ;     
            fillStr->Ak(0, 0)[(row + 2) * N_e + col + 2] +=  (jac*h_*sus40)*dq2dq2dsus4_IJ ;       
            fillStr->Ak(0, 0)[(row + 2) * N_e + col + 3] +=  (jac*dv1dh_J*sus40)*dq2dsus4_I ;  
            fillStr->Ak(0, 0)[(row + 2) * N_e + col + 4] +=  (jac*dv2dh_J*sus40)*dq2dsus4_I ;   


            fillStr->Ak(0, 0)[(row + 3) * N_e + col + 0] +=  (jac*dhdv1dh_IJ*sus40)*S4 ; 
            fillStr->Ak(0, 0)[(row + 3) * N_e + col + 1] +=  (jac*dv1dh_I*sus40)* dq1dsus4_J ; 
            fillStr->Ak(0, 0)[(row + 3) * N_e + col + 2] +=  (jac*dv1dh_I*sus40)* dq2dsus4_J ;  
            fillStr->Ak(0, 0)[(row + 3) * N_e + col + 3] +=  0.0; 
            fillStr->Ak(0, 0)[(row + 3) * N_e + col + 4] +=  0.0;


            fillStr->Ak(0, 0)[(row + 4) * N_e + col + 0] +=  (jac*dhdv2dh_IJ*sus40)*(2*(q[0]*q[0] + q[1]*q[1] ))*(2*(q[0]*q[0] + q[1]*q[1] )); ;   
            fillStr->Ak(0, 0)[(row + 4) * N_e + col + 1] +=  (jac*dv2dh_I*sus40)*2*(2*(q[0]*q[0] + q[1]*q[1] ))*(4*(q[0]*bf[j] ));       
            fillStr->Ak(0, 0)[(row + 4) * N_e + col + 2] +=  (jac*dv2dh_I*sus40)*2*(2*(q[0]*q[0] + q[1]*q[1] ))*(4*(q[1]*bf[j] ));    
            fillStr->Ak(0, 0)[(row + 4) * N_e + col + 3] +=  0.0 ; 
            fillStr->Ak(0, 0)[(row + 4) * N_e + col + 4] +=  0.0 ;                                           

*/


//            fillStr->Ak(0, 0)[(row + 1) * N_e + col + 0] +=  (jac*bf[j]*sus40)*dq1dsus4_I ;  
//            fillStr->Ak(0, 0)[(row + 1) * N_e + col + 1] +=  (jac*h*sus40)*dq1dq1dsus4_IJ ;  
//            fillStr->Ak(0, 0)[(row + 1) * N_e + col + 2] +=  (jac*h*sus40)*dq1dq2dsus4_IJ ;  


//            fillStr->Ak(0, 0)[(row + 2) * N_e + col + 0] +=  (jac*bf[j]*sus40)*dq2dsus4_I ;
//            fillStr->Ak(0, 0)[(row + 2) * N_e + col + 1] +=  (jac*h*sus40)*dq2dq1dsus4_IJ ; 
//            fillStr->Ak(0, 0)[(row + 2) * N_e + col + 2] +=  (jac*h*sus40)*dq2dq2dsus4_IJ ;          




           // ROTATIONAL VISCOSITY  

            fillStr->Ak(0, 0)[(row + 1) * N_e + col + 0] += (0.5 * jac * dhh_J * rvisc)*4*(dq1dJq_I[0]*Jq[0] + dq1dJq_I[1]*Jq[1]);
            fillStr->Ak(0, 0)[(row + 1) * N_e + col + 1] += (0.5 * jac * h * rvisc) *4*(dq1dJq_I[0]*dq1dJq_J[0] + dq1dJq_I[1]*dq1dJq_J[1]);
            fillStr->Ak(0, 0)[(row + 1) * N_e + col + 2] += (0.5 * jac * h * rvisc) *4*(dq1dJq_I[0]*dq2dJq_J[0] + dq1dJq_I[1]*dq2dJq_J[1]);
            fillStr->Ak(0, 0)[(row + 1) * N_e + col + 3] += (0.5 * jac * h * rvisc) *4*(dq1dJq_I[0]*dv1dJq_J[0] + dq1dJq_I[1]*dv1dJq_J[1]) + (0.5 * jac * dv1h_J * rvisc) *4*(dq1dJq_I[0]*Jq[0] + dq1dJq_I[1]*Jq[1]);;
            fillStr->Ak(0, 0)[(row + 1) * N_e + col + 4] += (0.5 * jac * h * rvisc) *4*(dq1dJq_I[0]*dv2dJq_J[0] + dq1dJq_I[1]*dv2dJq_J[1]) + (0.5 * jac * dv2h_J * rvisc) *4*(dq1dJq_I[0]*Jq[0] + dq1dJq_I[1]*Jq[1]);;


            fillStr->Ak(0, 0)[(row + 2) * N_e + col + 0] +=  (0.5 * jac * dhh_J * rvisc) *4*(dq2dJq_I[0]*Jq[0] + dq2dJq_I[1]*Jq[1]);  
            fillStr->Ak(0, 0)[(row + 2) * N_e + col + 1] +=  (0.5 * jac * h * rvisc) *4*(dq2dJq_I[0]*dq1dJq_J[0] + dq2dJq_I[1]*dq1dJq_J[1]);
            fillStr->Ak(0, 0)[(row + 2) * N_e + col + 2] +=  (0.5 * jac * h * rvisc) *4*(dq2dJq_I[0]*dq2dJq_J[0] + dq2dJq_I[1]*dq2dJq_J[1]);;
            fillStr->Ak(0, 0)[(row + 2) * N_e + col + 3] +=  (0.5 * jac * h * rvisc) *4*(dq2dJq_I[0]*dv1dJq_J[0] + dq2dJq_I[1]*dv1dJq_J[1]) + (0.5 * jac * dv1h_J * rvisc) *4*(dq2dJq_I[0]*Jq[0] + dq2dJq_I[1]*Jq[1]); ; 
            fillStr->Ak(0, 0)[(row + 2) * N_e + col + 4] +=  (0.5 * jac * h * rvisc) *4*(dq2dJq_I[0]*dv2dJq_J[0] + dq2dJq_I[1]*dv2dJq_J[1]) + (0.5 * jac * dv2h_J * rvisc) *4*(dq2dJq_I[0]*Jq[0] + dq2dJq_I[1]*Jq[1]); 


            fillStr->Ak(0, 0)[(row + 3) * N_e + col + 0] +=  (0.5 * jac * dhh_J * rvisc) *4*(dv1dJq_I[0]*Jq[0] + dv1dJq_I[1]*Jq[1]);
            fillStr->Ak(0, 0)[(row + 3) * N_e + col + 1] +=  (0.5 * jac * h * rvisc) *4*(dq1dv1dJq_IJ[0]*Jq[0] + dv1dJq_I[0]*dq1dJq_J[0]  + dq1dv1dJq_IJ[1]*Jq[1] + dv1dJq_I[1]*dq1dJq_J[1]);
            fillStr->Ak(0, 0)[(row + 3) * N_e + col + 2] +=  (0.5 * jac * h * rvisc) *4*(dq2dv1dJq_IJ[0]*Jq[0] + dv1dJq_I[0]*dq2dJq_J[0]  + dq2dv1dJq_IJ[1]*Jq[1] + dv1dJq_I[1]*dq2dJq_J[1]) ;  
            fillStr->Ak(0, 0)[(row + 3) * N_e + col + 3] +=  (0.5 * jac * h * rvisc) *4*(dv1dJq_I[0]*dv1dJq_J[0] + dv1dJq_I[1]*dv1dJq_J[1]) + (0.5 * jac * dv1h_J * rvisc) *4*(dv1dJq_I[0]*Jq[0] + dv1dJq_I[1]*Jq[1]);;  
            fillStr->Ak(0, 0)[(row + 3) * N_e + col + 4] +=  (0.5 * jac * h * rvisc) *4*(dv1dJq_I[0]*dv2dJq_J[0] + dv1dJq_I[1]*dv2dJq_J[1]) + (0.5 * jac * dv2h_J * rvisc) *4*(dv1dJq_I[0]*Jq[0] + dv1dJq_I[1]*Jq[1]);;              


            fillStr->Ak(0, 0)[(row + 4) * N_e + col + 0] +=  (0.5 * jac * dhh_J * rvisc) *4*(dv2dJq_I[0]*Jq[0] + dv2dJq_I[1]*Jq[1]);
            fillStr->Ak(0, 0)[(row + 4) * N_e + col + 1] +=  (0.5 * jac * h * rvisc) *4*(dq1dv2dJq_IJ[0]*Jq[0] + dv2dJq_I[0]*dq1dJq_J[0]  + dq1dv2dJq_IJ[1]*Jq[1] +  dv2dJq_I[1]*dq1dJq_J[1]);
            fillStr->Ak(0, 0)[(row + 4) * N_e + col + 2] +=  (0.5 * jac * h * rvisc) *4*(dq2dv2dJq_IJ[0]*Jq[0] + dv2dJq_I[0]*dq2dJq_J[0]  + dq2dv2dJq_IJ[1]*Jq[1] +  dv2dJq_I[1]*dq2dJq_J[1]);   
            fillStr->Ak(0, 0)[(row + 4) * N_e + col + 3] +=  (0.5 * jac * h * rvisc) *4*(dv2dJq_I[0]*dv1dJq_J[0] + dv2dJq_I[1]*dv1dJq_J[1]) + (0.5 * jac * dv1h_J * rvisc) *4*(dv2dJq_I[0]*Jq[0] + dv2dJq_I[1]*Jq[1]);  
            fillStr->Ak(0, 0)[(row + 4) * N_e + col + 4] +=  (0.5 * jac * h * rvisc) *4*(dv2dJq_I[0]*dv2dJq_J[0] + dv2dJq_I[1]*dv2dJq_J[1]) + (0.5 * jac * dv2h_J * rvisc) *4*(dv2dJq_I[0]*Jq[0] + dv2dJq_I[1]*Jq[1]);               




           // COUPLED TERM 
            fillStr->Ak(0, 0)[(row + 1) * N_e + col + 0] +=  (jac * dhh_J * cvisc) *(dq1dJq_I[0]*rodt[0] + 2*dq1dJq_I[1]*rodt[1]  -  dq1dJq_I[0]*rodt[3] );
            fillStr->Ak(0, 0)[(row + 1) * N_e + col + 3] +=  (jac * h * cvisc) *(dq1dJq_I[0]*dv1rodt_J[0] + 2*dq1dJq_I[1]*dv1rodt_J[1]  -  dq1dJq_I[0]*dv1rodt_J[3] )+(jac * dv1h_J * cvisc) *(dq1dJq_I[0]*rodt[0] + 2*dq1dJq_I[1]*rodt[1]  -  dq1dJq_I[0]*rodt[3] ); 
            fillStr->Ak(0, 0)[(row + 1) * N_e + col + 4] +=  (jac * h * cvisc) *(dq1dJq_I[0]*dv2rodt_J[0] + 2*dq1dJq_I[1]*dv2rodt_J[1]  -  dq1dJq_I[0]*dv2rodt_J[3] )+(jac * dv2h_J * cvisc) *(dq1dJq_I[0]*rodt[0] + 2*dq1dJq_I[1]*rodt[1]-dq1dJq_I[0]*rodt[3] ); ;  


            fillStr->Ak(0, 0)[(row + 2) * N_e + col + 0] +=  (jac * dhh_J * cvisc) *(dq2dJq_I[0]*rodt[0] + 2*dq2dJq_I[1]*rodt[1]  -  dq2dJq_I[0]*rodt[3] );
            fillStr->Ak(0, 0)[(row + 2) * N_e + col + 3] +=  (jac * h * cvisc) *(dq2dJq_I[0]*dv1rodt_J[0] + 2*dq2dJq_I[1]*dv1rodt_J[1]  -  dq2dJq_I[0]*dv1rodt_J[3] ) + (jac * dv1h_J * cvisc) *(dq2dJq_I[0]*rodt[0] + 2*dq2dJq_I[1]*rodt[1]  -  dq2dJq_I[0]*rodt[3] ); ;  
            fillStr->Ak(0, 0)[(row + 2) * N_e + col + 4] +=  (jac * h * cvisc) *(dq2dJq_I[0]*dv2rodt_J[0] + 2*dq2dJq_I[1]*dv2rodt_J[1]  -  dq2dJq_I[0]*dv2rodt_J[3] ) + (jac * dv2h_J * cvisc) *(dq2dJq_I[0]*rodt[0] + 2*dq2dJq_I[1]*rodt[1]  -  dq2dJq_I[0]*rodt[3] ); ;  


            fillStr->Ak(0, 0)[(row + 3) * N_e + col + 0] +=  (jac * dhh_J * cvisc) *((dv1dJq_I[0]*rodt[0] + 2*dv1dJq_I[1]*rodt[1] -  dv1dJq_I[0]*rodt[3] ) + (Jq[0]*dv1rodtI[0] + 2*Jq[1]*dv1rodtI[1] -  Jq[0]*dv1rodtI[3] ));
            fillStr->Ak(0, 0)[(row + 3) * N_e + col + 1] +=  (jac * h * cvisc) *((dq1dv1dJq_IJ[0]*rodt[0] + 2*dq1dv1dJq_IJ[1]*rodt[1] -  dq1dv1dJq_IJ[0]*rodt[3] )  + (dq1dJq_J[0]*dv1rodtI[0] + 2*dq1dJq_J[1]*dv1rodtI[1] -  dq1dJq_J[0]*dv1rodtI[3] ));
            fillStr->Ak(0, 0)[(row + 3) * N_e + col + 2] +=  (jac * h * cvisc) *((dq2dv1dJq_IJ[0]*rodt[0] + 2*dq2dv1dJq_IJ[1]*rodt[1] -  dq2dv1dJq_IJ[0]*rodt[3] )+  (dq2dJq_J[0]*dv1rodtI[0] + 2*dq2dJq_J[1]*dv1rodtI[1]  -  dq2dJq_J[0]*dv1rodtI[3] ));
            fillStr->Ak(0, 0)[(row + 3) * N_e + col + 3] +=  (jac * h * cvisc) *((dv1dJq_I[0]*dv1rodt_J[0] + 2*dv1dJq_I[1]*dv1rodt_J[1] -  dv1dJq_I[0]*dv1rodt_J[3] ) + (dv1dJq_J[0]*dv1rodtI[0] + 2*dv1dJq_J[1]*dv1rodtI[1] - dv1dJq_J[0]*dv1rodtI[3] )) + 
                                                             (jac * dv1h_J * cvisc) *((dv1dJq_I[0]*rodt[0] + 2*dv1dJq_I[1]*rodt[1] -  dv1dJq_I[0]*rodt[3] ) + (Jq[0]*dv1rodtI[0] + 2*Jq[1]*dv1rodtI[1] -  Jq[0]*dv1rodtI[3] ));  
            fillStr->Ak(0, 0)[(row + 3) * N_e + col + 4] +=  (jac * h * cvisc) *((dv1dJq_I[0]*dv2rodt_J[0] + 2*dv1dJq_I[1]*dv2rodt_J[1] -  dv1dJq_I[0]*dv2rodt_J[3] ) + (dv2dJq_J[0]*dv1rodtI[0] + 2*dv2dJq_J[1]*dv1rodtI[1] -  dv2dJq_J[0]*dv1rodtI[3] ))  + 
                                                             (jac * dv2h_J * cvisc) *((dv1dJq_I[0]*rodt[0] + 2*dv1dJq_I[1]*rodt[1] -  dv1dJq_I[0]*rodt[3] ) + (Jq[0]*dv1rodtI[0] + 2*Jq[1]*dv1rodtI[1] -  Jq[0]*dv1rodtI[3] ));    
                                                                                               

            fillStr->Ak(0, 0)[(row + 4) * N_e + col + 0] +=  (jac * dhh_J * cvisc) *((dv2dJq_I[0]*rodt[0] + 2*dv2dJq_I[1]*rodt[1] -  dv2dJq_I[0]*rodt[3] ) + (Jq[0]*dv2rodtI[0] + 2*Jq[1]*dv2rodtI[1] -  Jq[0]*dv2rodtI[3] ));
            fillStr->Ak(0, 0)[(row + 4) * N_e + col + 1] +=  (jac * h * cvisc) *((dq1dv2dJq_IJ[0]*rodt[0] + 2*dq1dv2dJq_IJ[1]*rodt[1] -  dq1dv2dJq_IJ[0]*rodt[3] ) + (dq1dJq_J[0]*dv2rodtI[0] + 2*dq1dJq_J[1]*dv2rodtI[1] -  dq1dJq_J[0]*dv2rodtI[3] ));
            fillStr->Ak(0, 0)[(row + 4) * N_e + col + 2] +=  (jac * h * cvisc) *((dq2dv2dJq_IJ[0]*rodt[0] + 2*dq2dv2dJq_IJ[1]*rodt[1] - dq2dv2dJq_IJ[0]*rodt[3] ) + (dq2dJq_J[0]*dv2rodtI[0] + 2*dq2dJq_J[1]*dv2rodtI[1] -  dq2dJq_J[0]*dv2rodtI[3] ));
            fillStr->Ak(0, 0)[(row + 4) * N_e + col + 3] +=  (jac * h * cvisc) *((dv2dJq_I[0]*dv1rodt_J[0] + 2*dv2dJq_I[1]*dv1rodt_J[1] -  dv2dJq_I[0]*dv1rodt_J[3] ) + (dv1dJq_J[0]*dv2rodtI[0] + 2*dv1dJq_J[1]*dv2rodtI[1] -  dv1dJq_J[0]*dv2rodtI[3] )) + 
                                                             (jac * dv1h_J * cvisc) *((dv2dJq_I[0]*rodt[0] + 2*dv2dJq_I[1]*rodt[1] -  dv2dJq_I[0]*rodt[3] ) + (Jq[0]*dv2rodtI[0] + 2*Jq[1]*dv2rodtI[1] -  Jq[0]*dv2rodtI[3] ));
 
            fillStr->Ak(0, 0)[(row + 4) * N_e + col + 4] +=  (jac * h * cvisc) *((dv2dJq_I[0]*dv2rodt_J[0] + 2*dv2dJq_I[1]*dv2rodt_J[1] -  dv2dJq_I[0]*dv2rodt_J[3] ) + (dv2dJq_J[0]*dv2rodtI[0] + 2*dv2dJq_J[1]*dv2rodtI[1] -  dv2dJq_J[0]*dv2rodtI[3] )) + 
                                                             (jac * dv2h_J * cvisc) *((dv2dJq_I[0]*rodt[0] + 2*dv2dJq_I[1]*rodt[1] -  dv2dJq_I[0]*rodt[3] ) + (Jq[0]*dv2rodtI[0] + 2*Jq[1]*dv2rodtI[1] -  Jq[0]*dv2rodtI[3] ));;     
                                                                                                                                                                                                                                                                           
           



           // FRANK CONSTANT



           fillStr->Ak(0, 0)[(row + 1) * N_e + col + 0] += (jac * bf[j]*frank) *(dq1DqDq_I);;
           fillStr->Ak(0, 0)[(row + 1) * N_e + col + 1] += (jac * h*frank) *(dq1dq1DqDq_IJ); 
           fillStr->Ak(0, 0)[(row + 1) * N_e + col + 2] += 0.0;
//           fillStr->Ak(0, 0)[(row + 1) * N_e + col + 3] += (jac * dv1dh_J*frank) *(dq1DqDq_I);
//           fillStr->Ak(0, 0)[(row + 1) * N_e + col + 4] += (jac * dv2dh_J*frank) *(dq1DqDq_I);

           fillStr->Ak(0, 0)[(row + 2) * N_e + col + 0] +=  (jac * bf[j]*frank) *dq2DqDq_I;
           fillStr->Ak(0, 0)[(row + 2) * N_e + col + 1] +=  (jac * h*frank) *dq1dq2DqDq_IJ;     
           fillStr->Ak(0, 0)[(row + 2) * N_e + col + 2] +=  (jac * h*frank) *dq2dq2DqDq_IJ;     
//           fillStr->Ak(0, 0)[(row + 2) * N_e + col + 3] +=  (jac * dv1dh_J*frank) *dq2DqDq_I;
//           fillStr->Ak(0, 0)[(row + 2) * N_e + col + 4] +=  (jac * dv2dh_J*frank) *dq2DqDq_I;


           fillStr->Ak(0, 0)[(row + 3) * N_e + col + 0] +=  (jac * dhdv1dh_IJ*frank) *DqDq; ;
           fillStr->Ak(0, 0)[(row + 3) * N_e + col + 1] +=  (jac * dv1dh_I*frank) *dq1DqDq_J; ;  
           fillStr->Ak(0, 0)[(row + 3) * N_e + col + 2] +=  (jac * dv1dh_I*frank) *dq2DqDq_J ;  



           fillStr->Ak(0, 0)[(row + 4) * N_e + col + 0] +=  (jac * dhdv2dh_IJ*frank) *DqDq; ;
           fillStr->Ak(0, 0)[(row + 4) * N_e + col + 1] +=  (jac * dv2dh_I*frank) *dq1DqDq_J; ;  
           fillStr->Ak(0, 0)[(row + 4) * N_e + col + 2] +=  (jac * dv2dh_I*frank) *dq2DqDq_J ;  




/*
           fillStr->Ak(0, 0)[(row + 1) * N_e + col + 0] += h*(jac * dhdh_J*lambda_rot) *2*(q[0]*dq1dJq_I[0]+ q[1]*dq1dJq_I[1])  + bf[j]*(jac * h_*lambda_rot) *2*(q[0]*dq1dJq_I[0]+ q[1]*dq1dJq_I[1]);;
           fillStr->Ak(0, 0)[(row + 1) * N_e + col + 1] += h*(jac *h_*lambda_rot) *2*(bf[j]*dq1dJq_I[0]); 
           fillStr->Ak(0, 0)[(row + 1) * N_e + col + 2] += h*(jac *h_*lambda_rot) *2*(bf[j]*dq1dJq_I[1]);                                                       
           fillStr->Ak(0, 0)[(row + 1) * N_e + col + 3] += h*(jac * dv1dh_J*lambda_rot) *2*(q[0]*dq1dJq_I[0]+ q[1]*dq1dJq_I[1]);
           fillStr->Ak(0, 0)[(row + 1) * N_e + col + 4] += h*(jac * dv2dh_J*lambda_rot) *2*(q[0]*dq1dJq_I[0]+ q[1]*dq1dJq_I[1]);


           fillStr->Ak(0, 0)[(row + 2) * N_e + col + 0] += h*(jac * dhdh_J*lambda_rot) *2*(q[0]*dq2dJq_I[0]+ q[1]*dq2dJq_I[1]) + bf[j]*(jac * h_*lambda_rot) *2*(q[0]*dq2dJq_I[0]+ q[1]*dq2dJq_I[1]);;
           fillStr->Ak(0, 0)[(row + 2) * N_e + col + 1] += h*(jac *h_*lambda_rot) *2*(bf[j]*dq2dJq_I[0]); 
           fillStr->Ak(0, 0)[(row + 2) * N_e + col + 2] += h*(jac *h_*lambda_rot) *2*(bf[j]*dq2dJq_I[1]);                                
           fillStr->Ak(0, 0)[(row + 2) * N_e + col + 3] += h*(jac * dv1dh_J*lambda_rot) *2*(q[0]*dq2dJq_I[0]+ q[1]*dq2dJq_I[1]);
           fillStr->Ak(0, 0)[(row + 2) * N_e + col + 4] += h*(jac * dv2dh_J*lambda_rot) *2*(q[0]*dq2dJq_I[0]+ q[1]*dq2dJq_I[1]);



           fillStr->Ak(0, 0)[(row + 3) * N_e + col + 0] +=h*(jac *dhdh_J*lambda_rot) *2*(q[0]*dv1dJq_I[0]+ q[1]*dv1dJq_I[1]) + bf[j]*(jac *h_*lambda_rot) *2*(q[0]*dv1dJq_I[0]+ q[1]*dv1dJq_I[1]);;  
           fillStr->Ak(0, 0)[(row + 3) * N_e + col + 1] +=h*(jac * h_*lambda_rot) *2*(q[0]*dq1dv1dJq_IJ[0]+ q[1]*dq1dv1dJq_IJ[1]  + bf[j]*dv1dJq_I[0]); 
           fillStr->Ak(0, 0)[(row + 3) * N_e + col + 2] +=h*(jac * h_*lambda_rot) *2*(q[0]*dq2dv1dJq_IJ[0]+ q[1]*dq2dv1dJq_IJ[1]  + bf[j]*dv1dJq_I[1]);
           fillStr->Ak(0, 0)[(row + 3) * N_e + col + 3] +=h*(jac * dv1dh_J*lambda_rot) *2*(q[0]*dv1dJq_I[0]+ q[1]*dv1dJq_I[1]);  
           fillStr->Ak(0, 0)[(row + 3) * N_e + col + 4] +=h*(jac * dv2dh_J*lambda_rot) *2*(q[0]*dv1dJq_I[0]+ q[1]*dv1dJq_I[1]);  




           fillStr->Ak(0, 0)[(row + 4) * N_e + col + 0] +=h*(jac * dhdh_J*lambda_rot) *2*(q[0]*dv2dJq_I[0]+ q[1]*dv2dJq_I[1]) + bf[j]*(jac * h_*lambda_rot) *2*(q[0]*dv2dJq_I[0]+ q[1]*dv2dJq_I[1]);;  
           fillStr->Ak(0, 0)[(row + 4) * N_e + col + 1] +=h*(jac * h_*lambda_rot) *2*(q[0]*dq1dv2dJq_IJ[0]+ q[1]*dq1dv2dJq_IJ[1] + bf[j]*dv2dJq_I[0]);
           fillStr->Ak(0, 0)[(row + 4) * N_e + col + 2] +=h*(jac * h_*lambda_rot) *2*(q[0]*dq2dv2dJq_IJ[0]+ q[1]*dq2dv2dJq_IJ[1] + bf[j]*dv2dJq_I[1]);
           fillStr->Ak(0, 0)[(row + 4) * N_e + col + 3] +=h*(jac * dv1dh_J*lambda_rot) *2*(q[0]*dv2dJq_I[0]+ q[1]*dv2dJq_I[1]);  
           fillStr->Ak(0, 0)[(row + 4) * N_e + col + 4] +=h*(jac * dv2dh_J*lambda_rot) *2*(q[0]*dv2dJq_I[0]+ q[1]*dv2dJq_I[1]);                              
                                                                                                                                                                                                                  

*/

           //Power Terms
           fillStr->Ak(0, 0)[(row + 3) * N_e + col + 0]+= (jac * bf[j] * lambda_iso)*(1.0 - S) *Dbf_gI[0];
           fillStr->Ak(0, 0)[(row + 4) * N_e + col + 0]+= (jac * bf[j] * lambda_iso)*(1.0 - S) *Dbf_gI[1];
           fillStr->Ak(0, 0)[(row + 3) * N_e + col + 0]+= (jac * bf[j]*crossLinker) *(Dbf_gI[0]*Snn0_[0]  + Dbf_gI[1]*Snn0_[1]) ;
           fillStr->Ak(0, 0)[(row + 4) * N_e + col + 0]+= (jac * bf[j]*crossLinker) *(Dbf_gI[1]*Snn0_[2]  + Snn0_[1]*Dbf_gI[0]);
           fillStr->Ak(0, 0)[(row + 3) * N_e + col + 0]+= (jac * bf[j]*lambda_aniso) *(Dbf_gI[0]*Snn0[0]  + Dbf_gI[1]*Snn0[1]); 
           fillStr->Ak(0, 0)[(row + 4) * N_e + col + 0]+= (jac * bf[j]*lambda_aniso) *(Dbf_gI[1]*Snn0[2]  + Snn0[1]*Dbf_gI[0]);                








         }



      }  








    tensor<double,2> nborCoords_(subFill.nborCoords.data(),eNN,nDim);
    tensor<double,2> nborDOFs_(subFill.nborDOFs.data(),eNN,numDOFs);
    tensor<double,2> nborDOFs0_(subFill.nborDOFs0.data(),eNN,numDOFs);

    tensor<double,1> bf_(subFill.getDer(0),eNN);
    tensor<double,2> Dbf_l_(subFill.getDer(1),eNN,pDim);
    //tensor<double,2> Dbf_g(eNN,pDim);



    tensor<double,2> Bk(fillStr->Bk(0).data(),eNN,numDOFs);
    tensor<double,4> Ak(fillStr->Ak(0,0).data(),eNN,numDOFs,eNN,numDOFs);

    //[3] GEOMETRY
    tensor<double,1> x_ = bf_ * nborCoords_;
    //Transformation matrix (from reference to spatial)
    tensor<double,2> T = nborCoords_(all,range(0,1)).T() * Dbf_l_;
    //Jacobian of the transformation (for integration)

//    double jac = T.det();

    //Global derivatives (wrt x and y ) of the basis functions
//    GlobalBasisFunctions::hessians(Dbf_g, jac, DDbf_g, Lapbf_g, eNN, pDim, nDim, nborCoords, Dbf_l, DDbf_l);
    tensor<double,2>  Dbf_g_ = Dbf_l_*T.inv();
    //-----------------------------------------------------------
    //[4] VARIABLES

    //[4.1] AUXILIARY
    tensor<double,1> nbor_h   = nborDOFs_(all,0);
    tensor<double,1> nbor_h0  = nborDOFs0_(all,0);
    tensor<double,2> nbor_q   = nborDOFs_(all,range(1,2));
    tensor<double,2> nbor_q0  = nborDOFs0_(all,range(1,2));
    tensor<double,2> nbor_v   = nborDOFs_(all,range(3,4));

    tensor<double,2> id2d = {{1.0,0.0},{0.0,1.0}};
    tensor<double,3> voigt = {{{1,0},{0,1}},
                              {{0,1},{-1,0}}};
    tensor<double,4> bf_voigt = outer(bf_,voigt.transpose({2,0,1}));
    tensor<double,5> Dbf_voigt = outer(Dbf_g_,voigt.transpose({2,0,1})).transpose({0,2,3,4,1});

    //[4.2] STATE VARIABLES

    tensor<double,1> Dh_   = nbor_h * Dbf_g_;
    tensor<double,1> Dh0_  = nbor_h0 * Dbf_g_;


    tensor<double,2> q_  = product(bf_voigt,nbor_q, {{0,0},{1,1}});
    tensor<double,2> q0_ = product(bf_voigt,nbor_q0,{{0,0},{1,1}});
    tensor<double,3> Dq_ = product(nbor_q,Dbf_voigt, {{0,0},{1,1}});
    tensor<double,3> Dq_0_ = product(nbor_q0,Dbf_voigt, {{0,0},{1,1}});



    //[4.3] VELOCITY
    tensor<double,1> v_  = bf_ * nbor_v;
    //Velocity gradient
    tensor<double,2> Dv_ = product(nbor_v,Dbf_g_,{{0,0}});
    tensor<double,4> dvDv_ =outer(Dbf_g_,id2d).transpose({0,3,1,2});
    double divv_ = trace(Dv_);



    tensor<double,2> omega_ = 0.5*(Dv_ - Dv_.transpose({1,0}));
    tensor<double,4> dvomega_ = 0.5 * ( outer(Dbf_g_,id2d).transpose({0,3,2,1}) - outer(Dbf_g_,id2d).transpose({0,3,1,2}) );



    tensor<double,2> Jq_ = (q_-q0_)/deltat + Dq_0_ * v_ - 2.0 * product(omega_,q0_,{{1,0}});
    tensor<double,4> dqJq = bf_voigt/deltat ;
    tensor<double,4> dvJq = outer(bf_,Dq_.transpose({2,0,1})) - 2.0 * product(dvomega_,q_,{{3,0}}) ;

    tensor<double,4> dvJq_ = outer(bf_,Dq_0_.transpose({2,0,1})) - 2.0 * product(dvomega_,q0_,{{3,0}}) ;


    double S2_ = 2.0*product(q_,q_,{{0,0},{1,1}});
    double S2_0_ = 2.0*product(q0_,q0_,{{0,0},{1,1}});
    double S_ = sqrt(S2_);
    // Defining a spatial field for the polymerization and de-polymerization rate




    //[5] FREE ENERGY (SUSCEPTIBILITY)
    double sus2{}, dsus2{}, ddsus2{};
    //computeSus2(sus2, dsus2, ddsus2, h_, sus20, hcrit);
    // double ccrit = hcrit ;
//    computeSus2(sus2, dsus2, ddsus2, 1., sus20, 1.0);

    double sus4{}, dsus4{}, ddsus4{};
//    computeSus4(sus4 , dsus4, ddsus4, 1., sus40, 1.);
    sus2=sus20;
    sus4=sus40;
//    cout<<sqrt(- (4*sus2*deltat+2*lambda_rot*h_)/ 32*sus4*deltat)<<endl;



    //computeSus4(sus4, dsus4, ddsus4, h_, sus40, hcrit);
    tensor<double,2> dvh_ = -deltat * ( outer(bf_,Dh_) + h * Dbf_g_ );


    tensor<double,1> dhh_ = -deltat * ( Dbf_g_ * v_ + bf_ * (divv_ + visc));
    tensor<double,3> dvdhh_ = -deltat * ( outer(bf_,Dbf_g_.T()) + outer(Dbf_g_,bf_) );
    tensor<double,3> dhdvh_ = dvdhh_.transpose({2,0,1});

    tensor<double,2> qbfvoigt = product(q_,bf_voigt,{{0,2},{1,3}});
    //[5.2] GRADIENT
    Bk(all,range(1,2)) += 4.0 * (jac * h * (sus2 + 2.0 * sus4 * S2_))  * qbfvoigt;
    Bk(all,range(3,4)) += 2.0*(0.5 * jac * S2_ * ( ( sus2 + sus4 * S2_) + h * (dsus2 + dsus4 * S2_))) * dvh_;


    //[5.3] HESSIAN
    Ak(all,range(1,2),all,range(1,2)) += 4.0 * jac * h * ((sus2 + 2.0 * sus4 * S2_) * product(bf_voigt,bf_voigt,{{2,2},{3,3}}) + 8.0 * sus4 * outer(qbfvoigt,qbfvoigt));
//    Ak(all,range(1,2),all,range(3,4)) += 4.0 * jac * ((sus2 + 2.0 * sus4 * S2_) + h_ * (dsus2 + 2.0 * dsus4 * S2_)) * outer(qbfvoigt, dvh_);
    Ak(all,range(3,4),all,range(1,2)) += 4.0 * jac * ((sus2 + 2.0 * sus4 * S2_) + h * (dsus2 + 2.0 * dsus4 * S2_)) * outer(dvh_,qbfvoigt);
//    Ak(all,range(3,4),all,range(3,4)) += 2.*(0.5 * jac * S2_ * ( 2.0 * ( dsus2 + dsus4 * S2_) + h_ * (ddsus2 + ddsus4 * S2_))) * outer(dvh_,dvh_);
    Ak(all,range(1,2),all,0)          += 4.0 * jac * ((sus2 + 2.0 * sus4 * S2_) + h * (dsus2 + 2.0 * dsus4 * S2_)) * outer(qbfvoigt, bf_);
    Ak(all,range(3,4),all,0)          += 2.0*(0.5 * jac * S2_ * ( 2.0 * ( dsus2 + dsus4 * S2_) + h * (ddsus2 + ddsus4 * S2_))) * outer(dvh_,bf_) + 2.0*(0.5 * jac * S2_ * ( ( sus2 + sus4 * S2_) + h * (dsus2 + dsus4 * S2_))) *  dvdhh_;


//    Bk(all,range(1,2)) += (jac *h* h_ * lambda_rot)* product(q_,dqJq,{{0,2},{1,3}});
//    Bk(all,range(3,4)) += (jac *h* h_ * lambda_rot)* product(q_,dvJq,{{0,2},{1,3}});                                                                                                                                                                                

    Bk(all,range(1,2)) += (jac *h* h * lambda_rot)* product(q0_,dqJq,{{0,2},{1,3}});
    Bk(all,range(3,4)) += (jac *h* h * lambda_rot)* product(q0_,dvJq_,{{0,2},{1,3}});         


    Ak(all,range(1,2),all,0)    += (jac *2* h * lambda_rot)* outer(product(q0_,dqJq,{{0,2},{1,3}}),bf_);
    Ak(all,range(3,4),all,0)    += (jac *2* h * lambda_rot)* outer(product(q0_,dvJq_,{{0,2},{1,3}}),bf_);        



//   Bk(all,range(3,4)) += c*((jac * h0 ) * divv_) * Dbf_g_;


//   Ak(all,range(3,4),all,range(3,4)) +=  c*(jac * h0 ) * outer(Dbf_g_,Dbf_g_);

    //postProcess
    fillStr->addContribGlobInteg("biLinear", product(outer(Dh_,Dh_),q_,{{0,0},{1,1}})*jac);   

    tensor<double,2> id2d_ = {{1.0,0.0},{0.0,1.0}};
    tensor<double,2> rodt_ = 0.5*(Dv_ + Dv_.transpose({1,0}));
    tensor<double,2> _Snn = q_ + 0.5 * sqrt(S2)*id2d;
    tensor<double,2> _Snn_ = -q_ + 0.5 * sqrt(S2)*id2d;
    //Components of antisymmetric stress
    tensor<double,2> nablaq = (rvisc*Jq_ + cvisc*rodt_ + lambda_rot*h*q_ + 4.0 * (sus20 + 2.*sus40 * S2)*q_ )/frank;
    tensor<double,2> Ten_sym =  jac*h*( 4*rodt_  + lambda_iso*(1-S)*id2d_ + lambda_aniso*_Snn  +  crossLinker*_Snn_  + cvisc*Jq_ - frank*product(Dq_,Dq_,{{0,0},{1,1}} ));                                                                                         
    tensor<double,2> Ten_antisym = 2*jac*h*frank*( product(q_,nablaq,{{1,1}}) -   product(nablaq,q_,{{1,1}})); 
    tensor<double,2> Ten = Ten_sym + Ten_antisym;
    tensor<double,2> Ten_partial = jac*h*( 4*rodt_  + lambda_iso*(1-S)*id2d_ +  crossLinker*_Snn_); 

    fillStr->addContribGlobInteg("sigma_xx", Ten_sym(0,0));   
    fillStr->addContribGlobInteg("sigma_yy", Ten_sym(1,1));   
    fillStr->addContribGlobInteg("sigma_partial_xx", Ten_partial(0,0)); 
    fillStr->addContribGlobInteg("sigma_partial_yy", Ten_partial(1,1)); 
    fillStr->addContribGlobInteg("gradV_xx", 4*h*jac*rodt_(0,0)); 
    fillStr->addContribGlobInteg("gradV_yy", 4*h*jac*rodt_(1,1)); 
    fillStr->addContribGlobInteg("fric_x",jac*v[0]*xCoords[0]/r_out); 
    fillStr->addContribGlobInteg("fric_y",jac*v[1]*xCoords[1]/r_out); 
    
    double d = (0.999-S2);

    tensor<double,2>  dqd   = -4*qbfvoigt; // (1./(S+1E-10))*qbfvoigt;
    tensor<double,4>  dqdqd = -4*product(bf_voigt,bf_voigt,{{2,2},{3,3}}) ; //(1./(S+1E-10))*product(bf_voigt,bf_voigt,{{2,2},{3,3}}); 
    //dqdqd +=   -(4./(S+1E-10)*(S+1E-10))*outer(qbfvoigt,qbfvoigt);
//4*qbfvoigt;
//    tensor<double,4> dqdqd  =4*product(bf_voigt,bf_voigt,{{2,2},{3,3}}) ;
double k_lim =1E-4/deltat;
//    double f_lim = k_lim/(d*d);
    double f_lim = -k_lim*log(d);

    tensor<double,2> dqf_lim     = -(k_lim/d)*dqd;
    tensor<double,4> dqdqf_lim   =  -(k_lim/d)*dqdqd + k_lim/(d*d)*outer(dqd,dqd);
//    tensor<double,2> dqf_lim     = -2*k_lim/(d*d*d)*dqd;
//    tensor<double,4> dqdqf_lim   =  6*(k_lim/(d*d*d*d))*outer(dqd,dqd)-2*(k_lim/(d*d*d))*dqdqd;
   
    Bk(all,range(1,2))          +=  jac*h*dqf_lim  ;
    Bk(all,range(3,4))          +=  jac*h*f_lim*dvh_;

    Ak(all,range(1,2),all,0)    +=  jac*outer(dqf_lim,bf_);  
//    Ak(all,range(1,2),all,range(1,2)) += (0.05*0.001)*(-1 *8* ( jac * h *(1/(d*d*d)) ) * product(bf_voigt,bf_voigt,{{2,2},{3,3}}) + 3 *4* 8*( jac * h 
    Ak(all,range(1,2),all,range(1,2)) += jac*h*dqdqf_lim ;                                                     
 //jac*h*dqdqf_lim;

//0.001*(-1 *4* ( jac * h *(1/(d*d*d)) ) * product(bf_voigt,bf_voigt,{{2,2},{3,3}}) + 3 *4* 4*( jac * h *(1/(d*d*d*d)) ) * outer(qbfvoigt,qbfvoigt))/del

    Ak(all,range(3,4),all,0)          +=  jac*f_lim*(outer(dvh_,bf_)+h*dvdhh_);
    Ak(all,range(3,4),all,range(1,2)) +=  jac*h*outer(dvh_,dqf_lim);                                                                                   
    


    /*
    double d = (S2_-1.);
   
    Bk(all,range(1,2))          += -0.001*1 *4* ( jac * h *(1/(d*d*d)) ) * qbfvoigt/deltat  ;

    Ak(all,range(1,2),all,0)    += -0.001*1 *4* ( jac *  (1/(d*d*d)) ) * outer(qbfvoigt/deltat,bf_)  ;  
    Ak(all,range(1,2),all,range(1,2)) += 0.001*(-1 *4* ( jac * h *(1/(d*d*d)) ) * product(bf_voigt,bf_voigt,{{2,2},{3,3}}) + 3 *4* 4*( jac * h *(1/(d*d*d*d)) ) * outer(qbfvoigt,qbfvoigt))/deltat ;      */

}

void LS_CortexNematic_Border(Teuchos::RCP<hiperlife::FillStructure> fillStr)
{
    using namespace hiperlife;
    using namespace hiperlife::Tensor;

    auto& subFill = (*fillStr)["cortexHand"];
    int nDim = subFill.nDim;
    int pDim = subFill.pDim;
    int eNN  = subFill.eNN;
    int numDOFs = subFill.numDOFs;

    //Material parameters
    double deltat = fillStr->userStr->dparam[0];
    double lambda_iso    = fillStr->userStr->dparam[7];
    double lambda_aniso    = fillStr->userStr->dparam[8];
    double visc    = fillStr->userStr->dparam[10];
    double fric    = fillStr->userStr->dparam[13];
    double r_in    = fillStr->userStr->dparam[16];
    double heqb    = fillStr->userStr->dparam[22];
    double sus20  = fillStr->userStr->dparam[1];
    double sus40  = fillStr->userStr->dparam[2];
    double rvisc   = fillStr->userStr->dparam[11];
    double cvisc   = fillStr->userStr->dparam[12];
    double frank  = fillStr->userStr->dparam[6];
    double lambda_rot    = fillStr->userStr->dparam[9];
    double lambda_trans = fillStr->userStr->dparam[26];                                  
    double r_out = fillStr->userStr->dparam[17];
    vector<double> tr = subFill.tangentsBoundaryRef();
    tensor<double,1> tangentRef(tr.data(),2);
    tensor<double,2> nborCoords(subFill.nborCoords.data(),eNN,nDim);
    tensor<double,2> nborDOFs0(subFill.nborDOFs0.data(),eNN,numDOFs);
    tensor<double,2> nborDOFs(subFill.nborDOFs.data(),eNN,numDOFs);

    tensor<double,1>  bf(subFill.getDer(0),eNN);
    tensor<double,2> Dbf_l(subFill.getDer(1),eNN,pDim);


    //[2] GEOMETRY
    tensor<double,1> x = bf * nborCoords;
    tensor<double,2> T = nborCoords(all,range(0,1)).T() * Dbf_l;
//    double jac = T.det();
    tensor<double,1> tangent = T * tangentRef;
    double normt = sqrt(tangent(0)*tangent(0)+tangent(1)*tangent(1));
    tensor<double,1> bnormal = {tangent(1),-tangent(0)};
    bnormal /= normt;
    tensor<double,2>  Dbf_g = Dbf_l*T.inv();
//[3] VARIABLES
    tensor<double,1> nbor_h  = nborDOFs(all,0);
    tensor<double,1> nbor_h0 = nborDOFs0(all,0);
    tensor<double,2> nbor_q  = nborDOFs(all,range(1,2));
    tensor<double,2> nbor_q0 = nborDOFs0(all,range(1,2));
    tensor<double,2> nbor_v = nborDOFs(all,range(3,4));
    tensor<double,2> id2d = {{1.0,0.0},{0.0,1.0}};
    tensor<double,3> voigt = {{{1,0},{0,1}},
                              {{0,1},{-1,0}}};
    tensor<double,4> bf_voigt = outer(bf,voigt.transpose({2,0,1}));
    tensor<double,5> Dbf_voigt = outer(Dbf_g,voigt.transpose({2,0,1})).transpose({0,2,3,4,1});


//[3.2] STATE/ PROCESS VARIABLES
    double h            = bf * nbor_h;
    tensor<double,1> v  = bf * nbor_v;
    tensor<double,2> Dv = product(nbor_v,Dbf_g,{{0,0}});
    tensor<double,2> rodt = 0.5*(Dv + Dv.transpose({1,0}));
    tensor<double,2> omega = 0.5*(Dv - Dv.transpose({1,0}));

    tensor<double,2> q  = product(bf_voigt,nbor_q, {{0,0},{1,1}});
    tensor<double,3> Dq = product(nbor_q,Dbf_voigt, {{0,0},{1,1}});
    tensor<double,2> q0 = product(bf_voigt,nbor_q0,{{0,0},{1,1}});
    double S2 = 2.0*product(q,q,{{0,0},{1,1}});
    double S  = sqrt(S2); 
    tensor<double,2> Snn = q + 0.5 * sqrt(S2)*id2d;
    tensor<double,2> Snn_ = -q + 0.5 * sqrt(S2)*id2d;
    tensor<double,2> Jq = (q-q0)/deltat + Dq * v - 2.0 * product(omega,q,{{1,0}}); 
    //Components of antisymmetric stress
    tensor<double,2> nablaq = (rvisc*Jq + cvisc*rodt + lambda_rot*h*q + 4.0 * (sus20 + 2.*sus40 * S2)*q )/frank;
  
    tensor<double,1> Ten_sym = normt*h*( product(rodt,bnormal,{{1,0}})  + lambda_iso*(1-S)*product(id2d,bnormal,{{1,0}}) + lambda_aniso*product(Snn,bnormal,{{1,0}})  +  lambda_trans*product(Snn_,bnormal,{{1,0}})   + cvisc*product(Jq,bnormal,{{1,0}} ) - frank*product(product(Dq,Dq,{{0,0},{1,1}} ),bnormal,{{1,0}} ));                                                                                           
    //cout<<"Hello"<<endl;
    //tensor<double,1> Ten_antisym = 2*normt*h*frank*( product(product(q,nablaq,{{1,1}}) -   product(nablaq,q,{{1,1}}),bnormal,{{1,0}})); 
    //tensor<double,1> Ten = Ten_sym + Ten_antisym;
    //if (x(0)>r_out-1e-3 and x(0)<1E-3)
        fillStr->addContribGlobInteg("Tx", 0.5*Ten_sym(0));
    //if (x(1)>r_out-1e-3 and x(1)<1E-3)
        fillStr->addContribGlobInteg("Ty", 0.5*Ten_sym(1) );

    //fillStr->addContribGlobInteg("Border", normt );   
}                                                                                                                                                                        

