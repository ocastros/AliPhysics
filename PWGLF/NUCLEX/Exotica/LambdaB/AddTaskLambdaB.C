#include "AliAnalysisTaskLambdaB.h"
#include "TString.h"

AliAnalysisTaskLambdaB *AddTaskLambdaB(bool isMC = false, TString suffix="") {
  return AliAnalysisTaskLambdaB::AddTask(isMC, suffix);
}
