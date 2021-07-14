#include "LambdaBstructures.h"

ClassImp(RLambdaB)
ClassImp(RLambdaB3O2)
ClassImp(RLambdaB3KF)
ClassImp(SLambdaB<RLambdaB3KF>)
ClassImp(SLambdaB<RLambdaB3O2>)

SLambdaB3O2 __dummy_instanceO2__()
{
  SLambdaB3O2 a;
  a.gPt = -12;
  return a;
}

SLambdaB3KF __dummy_instanceKF__()
{
  SLambdaB3KF a;
  a.gPt = -12;
  return a;
}
