#ifndef AuxCortexNematic2D
#define AuxCortexNematic2D

#include <Teuchos_RCP.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include "hl_FillStructure.h"
#include "hl_HiPerProblem.h"
#include "hl_Tensor.h"
#include "hl_DOFsHandler.h"
#include "hl_UserStructure.h"
#include <time.h>       /* time */

using namespace hiperlife;
using namespace hiperlife::Tensor;


void LS_CortexNematic_2D(Teuchos::RCP<hiperlife::FillStructure> fillStr);
void setCircleVelBC(Teuchos::RCP<hiperlife::DOFsHandler> dofHand, Teuchos::RCP<hiperlife::UserStructure> userStr);
void setCircleNemBC(Teuchos::RCP<hiperlife::DOFsHandler> dofHand, Teuchos::RCP<hiperlife::UserStructure> userStr);
void setCircleThiBC(Teuchos::RCP<hiperlife::DOFsHandler> dofHand, Teuchos::RCP<hiperlife::UserStructure> userStr);
void setLinearConstraints(std::vector<double>& dparam, Teuchos::RCP<hiperlife::HiPerProblem> hiperProbl);
void setCircleVelBC_square(Teuchos::RCP<hiperlife::DOFsHandler> dofHand, Teuchos::RCP<hiperlife::UserStructure> userStr);
void setCircleThiBC_square(Teuchos::RCP<hiperlife::DOFsHandler> dofHand, Teuchos::RCP<hiperlife::UserStructure> userStr);
void setCircleNemBC_square(Teuchos::RCP<hiperlife::DOFsHandler> dofHand, Teuchos::RCP<hiperlife::UserStructure> userStr);
void LS_CortexNematic_Border(Teuchos::RCP<hiperlife::FillStructure> fillStr);


inline double fRand(double fMin, double fMax)
{
    //srand (time(NULL));
    double f = (double)std::rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}                                                    

inline void computeSus2(double &sus2, double &dsus2, double &ddsus2, double sus20)
{
    sus2  = sus20;// * (hcrit - h);
    dsus2 = 0.0;//;
    ddsus2 = 0.0;
}

inline void computeSus4(double &sus4, double &dsus4, double &ddsus4, double sus40)
{

    sus4  = sus40;
    dsus4 = 0.0;
    ddsus4 = 0.0;
}                                                                                                                    





#endif 

