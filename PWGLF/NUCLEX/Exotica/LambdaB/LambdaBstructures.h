#ifndef LambdaB3STRUCTURES_H
#define LambdaB3STRUCTURES_H

#include <Rtypes.h>

struct RLambdaB {
  virtual ~RLambdaB() = default;
  float centrality = -1.;
  float pt = -999.f;
  float phi = -999.f;
  float pz = -999.f;
  float ct = -1.f;
  float r = -1.f;
  float cosPA = -2.f;
  float m = -1;
  float cosPA_Lambda = -2.; 
  Double32_t momDstar = 0;   //[0,10.22,8]
  Double32_t cosTheta_AProton1Aproton2 = 1.; //[-1,1,8]
  Double32_t cosThetaStar = 1.; //[-1,1,8]
  Double32_t mAp1Ap2_vert = -1.; //[1.059,1.1352,8]
  Double32_t mAp1Ap2 = -1.; //[1.150,1.531,8]
  Double32_t mHeAp2 = -1.; //[4.0,4.381,8]
  Double32_t dca_lambda_hyper = -1.0; //[0.0,8.0,8]
  Double32_t dca_He = -1.0; //[0.0,8.0,8]
  Double32_t dca_Apr1 = -1.0; //[0.0,8.0,8]
  Double32_t dca_Apr2 = -1.0; //[0.0,40.96,12]
  Double32_t tpcNsig_He = -4.0; //[-4.0,4.0,8]
  Double32_t tpcNsig_Apr1 = -4.0; //[-4.0,4.0,8]
  Double32_t tpcNsig_Apr2 = -4.0; //[-4.0,4.0,8]
  Double32_t tofNsig_He = -4.0; //[-4.0,4.0,8]
  Double32_t tofNsig_Apr1 = -4.0; //[-4.0,4.0,8]
  Double32_t tofNsig_Apr2 = -4.0; //[-4.0,4.0,8]
  Double32_t dca_He_Apr1 = -4.0; //[0.0,8.0,8]
  Double32_t dca_He_Apr2 = -4.0; //[0.0,8.0,8]
  Double32_t dca_Apr1_Apr2 = -4.0; //[0.0,8.0,8]
  UChar_t tpcClus_He = 0u;
  UChar_t tpcClus_Apr1 = 0u;
  UChar_t tpcClus_Apr2 = 0u;
  UChar_t its_clusmap_He = 0u;
  UChar_t its_clusmap_Apr1 = 0u;
  UChar_t its_clusmap_Apr2 = 0u;
  UChar_t candidates = 0u;
  UChar_t trigger = 0u;
  bool hasTOF_He = false;
  bool hasTOF_Apr1 = false;
  bool hasTOF_Apr2 = false;
  bool is_ITSrefit_He = false;
  bool is_ITSrefit_Apr1 = false;
  bool is_ITSrefit_Apr2 = false;
  bool positive = false;
  ClassDef(RLambdaB, 5)
};

struct RLambdaB3O2 : public RLambdaB {
  RLambdaB3O2() : RLambdaB{} {}
  virtual ~RLambdaB3O2() = default;
  Double32_t dca_He_sv = -4.0; //[0.0,8.0,8]
  Double32_t dca_Apr1_sv = -4.0; //[0.0,8.0,8]
  Double32_t dca_Apr2_sv = -4.0; //[0.0,8.0,8]
  Double32_t chi2 = -1.f;      //[0.0,32.,16]
  char rotation = 0;
  ClassDef(RLambdaB3O2,3)
};

struct RLambdaB3KF : public RLambdaB {
  RLambdaB3KF() : RLambdaB{} {}
  virtual ~RLambdaB3KF() = default;
  float chi2_HeApr1 = -1.f;
  float chi2_3prongs = -1.f;
  float chi2_topology = -1.f;
  ClassDef(RLambdaB3KF,1)
};

template<class Hyper>
struct SLambdaB : public Hyper {
  SLambdaB() : Hyper{} {}
  SLambdaB(const Hyper& other) : Hyper{other} {}
  virtual ~SLambdaB() = default;
  float gPt = -999.f;
  float gPhi = -999.f;
  float gPz = -999.f;
  float gCt = -1.f;
  float gT = -1.f;
  bool  gPositive = false;
  bool  gReconstructed = false;
  ClassDef(SLambdaB,1)
};

using SLambdaB3KF = SLambdaB<RLambdaB3KF>;
using SLambdaB3O2 = SLambdaB<RLambdaB3O2>;

#endif
