#include "AliAnalysisTaskLambdaB.h"

#include "AliAnalysisDataContainer.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliInputEventHandler.h"
#include "AliMCEvent.h"
#include "AliPDG.h"
#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliVVertex.h"

#include <Riostream.h>
#include <TChain.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TList.h>
#include <TRandom3.h>

#include <array>
#include <cmath>
#include <vector>
#include <unordered_map>

#include "Math/GenVector/Boost.h"
#include "Math/Vector3Dfwd.h"
#include "Math/Vector3D.h"
#include "Math/LorentzVector.h"

#include "AliDataFile.h"
#include <TFile.h>
#include <TSpline.h>

#include "Track.h"
#include <memory>

#define HomogeneousField
#include "KFParticle.h"
#include "KFVertex.h"
#include "KFPVertex.h"
#include "KFPTrack.h"

ClassImp(AliAnalysisTaskLambdaB);

namespace
{

  using lVector = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>>;

  struct HelperParticle
  {
    o2::track::TrackParCov *track = nullptr;
    int index = -1;
    float nSigmaTPC = -1.f;
    float nSigmaTOF = -1.f;
  };

  constexpr float kHeMass{1.87561};
  constexpr float kAP2Mass{0.938272};
  constexpr float kAP1Mass{0.938272};
  constexpr float kMasses[3]{kHeMass, kAP1Mass, kAP2Mass};
  constexpr AliPID::EParticleType kAliPID[3]{AliPID::kHe3, AliPID::kAProton, AliPID::kAProton};
  const int kPDGs[3]{AliPID::ParticleCode(kAliPID[0]), AliPID::ParticleCode(kAliPID[1]), AliPID::ParticleCode(kAliPID[2])};

  bool IsLambdaB(const AliVParticle *vPart, AliMCEvent *mcEvent)
  {
    int nDaughters = 0;

    int vPartPDG = vPart->PdgCode();
    int vPartLabel = vPart->GetLabel();

    if (!mcEvent->IsPhysicalPrimary(vPartLabel) || (std::abs(vPartPDG) != 5122))
      return false;

    for (int iD = vPart->GetDaughterFirst(); iD <= vPart->GetDaughterLast(); iD++)
    {
      AliVParticle *dPart = mcEvent->GetTrack(iD);

      int dPartPDG = dPart->PdgCode();
      if (std::abs(dPartPDG) != 11)
        nDaughters++;
    }
    if (nDaughters == 3)
      return true;
    return false;
  }

  int IsTrueLambdaBCandidate(AliESDtrack *t1, AliESDtrack *t2, AliESDtrack *t3, AliMCEvent *mcEvent)
  {
    if (!mcEvent)
      return 0;

    int lab1 = std::abs(t1->GetLabel());
    int lab2 = std::abs(t2->GetLabel());
    int lab3 = std::abs(t3->GetLabel());

    if (mcEvent->IsPhysicalPrimary(lab1))
      return -1;
    if (mcEvent->IsPhysicalPrimary(lab2))
      return -1;
    if (mcEvent->IsPhysicalPrimary(lab3))
      return -1;

    AliVParticle *part1 = mcEvent->GetTrack(lab1);
    AliVParticle *part2 = mcEvent->GetTrack(lab2);
    AliVParticle *part3 = mcEvent->GetTrack(lab3);

    if (!part1 || !part2 || !part3)
      return -1;

    int mom1 = part1->GetMother();
    int mom2 = part2->GetMother();
    int mom3 = part3->GetMother();

    if (mom1 != mom2 || mom1 != mom3 || mom2 != mom3)
      return -1;

    AliVParticle *mom = mcEvent->GetTrack(mom1);
    if (!mom)
      return -1;

    return (IsLambdaB(mom, mcEvent)) ? mom1 : -1;
  }

  bool HasTOF(AliVTrack *track)
  {
    const bool hasTOFout = track->GetStatus() & AliVTrack::kTOFout;
    const bool hasTOFtime = track->GetStatus() & AliVTrack::kTIME;
    return hasTOFout && hasTOFtime;
  }

} // namespace

AliAnalysisTaskLambdaB::~AliAnalysisTaskLambdaB(bool mc, std::string name)
    : AliAnalysisTaskSE(name.data()), fEventCuts{}, fVertexer{}, fVertexerLambda{}, fMC{mc}
{
  fTrackCuts.SetMinNClustersTPC(0);
  fTrackCuts.SetEtaRange(-0.9, 0.9);
  /// Settings for the custom vertexer

  /// Standard output
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class()); // Basic Histograms
  DefineOutput(2, TTree::Class()); // Hypertriton Candidates Tree output
}

AliAnalysisTaskLambdaB::~AliAnalysisTaskLambdaB()
{
  if (fListHist)
  {
    delete fListHist;
    fListHist = nullptr;
  }

  if (fTreeHyp3)
  {
    delete fTreeHyp3;
    fTreeHyp3 = nullptr;
  }

  if (fCosPAsplineFile)
  {
    delete fCosPAsplineFile;
  }

  if (fGenHypO2)
  {
    delete fGenHypO2;
  }
  else if (fRecHyp)
  {
    delete fRecHyp;
  }
}

void AliAnalysisTaskLambdaB::UserCreateOutputObjects()
{
  fCounter = 0;

  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  fInputHandler = (AliInputEventHandler *)(man->GetInputEventHandler());
  fPIDResponse = fInputHandler->GetPIDResponse();

  fInputHandler->SetNeedField();

  fListHist = new TList();
  fListHist->SetOwner(true);
  fEventCuts.AddQAplotsToList(fListHist);

  fHistNSigmaHelium = new TH2D("fHistNSigmaHelium", ";#it{p}_{T} (GeV/#it{c});n_{#sigma} TPC Helium; Counts", 100, 0., 10.,
                            80, -5.0, 5.0);
  fHistNSigmaAP1 =
      new TH2D("fHistNSigmaAP1", ";#it{p}_{T} (GeV/#it{c});n_{#sigma} TPC AntiProton 1; Counts", 100, 0., 10., 80, -5.0, 5.0);
  fHistNSigmaAP2 =
      new TH2D("fHistNSigmaAP2", ";#it{p}_{T} (GeV/#it{c});n_{#sigma} TPC AntiProton 2; Counts", 100, 0., 10., 80, -5.0, 5.0);

  fHistInvMass =
      new TH2D("fHistInvMass", ";m_{dp#pi}(GeV/#it{c^2}); #it{p}_{T} (GeV/#it{c}); Counts", 30, 2.96, 3.05, 100, 0, 10);

  fHistDecVertexRes =
      new TH1D("fHistDecVertexRes", "; Resoultion(cm); Counts", 40, -1, 1);
  fListHist->Add(fHistNSigmaHelium);
  fListHist->Add(fHistNSigmaAP1);
  fListHist->Add(fHistNSigmaAP2);
  fListHist->Add(fHistInvMass);
  fListHist->Add(fHistDecVertexRes);

  OpenFile(2);
  fTreeHyp3 = new TTree("LambdaB", "LambdaB 3 Body with the O2 Vertexer");

//!!
  if (fMC && man->GetMCtruthEventHandler())
  {
    fGenHypO2 = new SHyperTriton3O2;
    fRecHyp = (RHyperTriton3O2 *)fGenHypO2;
    fTreeHyp3->Branch("SHyperTriton", fGenHypO2);
  }
  else
  {
    fRecHyp = new RHyperTriton3O2;
    fTreeHyp3->Branch("RHyperTriton", static_cast<RHyperTriton3O2 *>(fRecHyp));
  }
  fCosPAsplineFile = TFile::Open(AliDataFile::GetFileName(fCosPAsplineName).data());
  if (fCosPAsplineFile)
  {
    fCosPAspline = (TSpline3 *)fCosPAsplineFile->Get("cutSpline");
  }

  PostData(1, fListHist);
  PostData(2, fTreeHyp3);

  AliPDG::AddParticlesToPdgDataBase();

} /// end UserCreateOutputObjects

void AliAnalysisTaskLambdaB::UserExec(Option_t *)
{

  AliESDEvent *esdEvent = dynamic_cast<AliESDEvent *>(InputEvent());
  if (!esdEvent)
  {
    ::Fatal("AliAnalysisTaskLambdaB::UserExec", "AliESDEvent not found.");
    return;
  }

  AliMCEvent *mcEvent = MCEvent();
  if (!mcEvent && fMC)
  {
    ::Fatal("AliAnalysisTaskLambdaB::UserExec", "Could not retrieve MC event");
    return;
  }

  if (!fEventCuts.AcceptEvent(esdEvent))
  {
    PostData(1, fListHist);
    PostData(2, fTreeHyp3);
    return;
  }

  if (!fMC && fDownscaling)
  {
    if (gRandom->Rndm() > fDownscalingFactorByEvent)
      return;
  }

  double pvPos[3], pvCov[6];
  fEventCuts.GetPrimaryVertex()->GetXYZ(pvPos);
  fEventCuts.GetPrimaryVertex()->GetCovarianceMatrix(pvCov);
  fRecHyp->centrality = fEventCuts.GetCentrality();

  fRecHyp->trigger = 0u;
  if (fInputHandler->IsEventSelected() & AliVEvent::kINT7)
    fRecHyp->trigger |= kINT7;
  if (fInputHandler->IsEventSelected() & AliVEvent::kCentral)
    fRecHyp->trigger |= kCentral;
  if (fInputHandler->IsEventSelected() & AliVEvent::kSemiCentral)
    fRecHyp->trigger |= kSemiCentral;
  if (fInputHandler->IsEventSelected() & AliVEvent::kHighMultV0)
    fRecHyp->trigger |= kHighMultV0;
  fRecHyp->trigger |= esdEvent->GetMagneticField() > 0 ? kPositiveB : 0;

  std::vector<HelperParticle> helpers[3][2];
  std::vector<AliESDtrack *> He3AP2Tracks[2][2];
  std::vector<AliESDtrack *> AP1AP2Tracks[2][2];
  std::vector<EventMixingTrack> heliumsForMixing;

  for (int iTrack = 0; iTrack < esdEvent->GetNumberOfTracks(); iTrack++)
  {
    AliESDtrack *track = esdEvent->GetTrack(iTrack);

    if (!track)
      continue;

    if (!fTrackCuts.AcceptTrack(track))
      continue;

    if (fMC && fOnlyTrueCandidates)
    {
      int lab = std::abs(track->GetLabel());
      if (!mcEvent->IsSecondaryFromWeakDecay(lab))
        continue;
      AliVParticle *part = mcEvent->GetTrack(lab);
      AliVParticle *moth = mcEvent->GetTrack(part->GetMother());
      if (std::abs(moth->PdgCode()) != 5122)
        continue;
    }

    bool candidate[3]{false, false, false};
    float nSigmasTPC[3]{-1., -1., -1.}, nSigmasTOF[3]{-1., -1., -1.};
    bool hasTOF{HasTOF(track)};
    float dca[2];
    track->GetImpactParameters(dca[0], dca[1]);
    double dcaNorm = std::hypot(dca[0], dca[1]);

    if (fUseCovarianceCut)
    {
      float cyy = track->GetSigmaY2(), czz = track->GetSigmaZ2(), cyz = track->GetSigmaZY();
      float detYZ = cyy * czz - cyz * cyz;
      if (detYZ < 0.)
        continue;
    }

    for (int iT{0}; iT < 3; ++iT)
    {
      nSigmasTPC[iT] = fPIDResponse->NumberOfSigmasTPC(track, kAliPID[iT]);
      nSigmasTOF[iT] = fPIDResponse->NumberOfSigmasTOF(track, kAliPID[iT]);
      bool requireTOFpid = track->P() > fRequireTOFpid[iT];
      if (std::abs(nSigmasTPC[iT]) < fTPCsigmas[iT] && dcaNorm > fMinTrackDCA[iT] && track->Pt() < fTrackPtRange[iT][1] &&
          track->Pt() > fTrackPtRange[iT][0] && track->GetTPCsignalN() >= fMinTPCpidClusters[iT])
        candidate[iT] = (std::abs(nSigmasTOF[iT]) < fTOFsigmas[iT]) || (!hasTOF && !requireTOFpid);
    }

    if (candidate[0] || candidate[1] || candidate[2])
    {
      HelperParticle helper;
      helper.track = static_cast<o2::track::TrackParCov *>((AliExternalTrackParam *)track);
      for (int iT{0}; iT < 3; ++iT)
      {
        if (candidate[iT])
        {
          int chargeIndex = (fSwapSign && iT == fMixingTrack) ? track->GetSigned1Pt() < 0 : track->GetSigned1Pt() > 0;
          helper.nSigmaTPC = nSigmasTPC[iT];
          helper.nSigmaTOF = nSigmasTOF[iT];
          if (iT == fMixingTrack && fEnableEventMixing)
            heliumsForMixing.emplace_back(track, nSigmasTPC[iT], nSigmasTOF[iT], 0);
          else
          {
            helpers[iT][chargeIndex].push_back(helper);
          }
        }
      }
    }
  }

  if (fEnableEventMixing)
  {
    auto mixingHeliums = GetEventMixingTracks(fEventCuts.GetCentrality(), pvPos[2]);
    for (auto mixTrack : mixingHeliums)
    {
      HelperParticle helper;
      AliESDtrack *track = &(mixTrack->track);
      helper.track = static_cast<o2::track::TrackParCov *>((AliExternalTrackParam *)track);
      int chargeIndex = track->GetSigned1Pt() > 0;
      helper.nSigmaTPC = mixTrack->nSigmaTPC;
      helper.nSigmaTOF = mixTrack->nSigmaTOF;
      helpers[fMixingTrack][chargeIndex].push_back(helper);
    }
  }

  lVector lambdab;
  ROOT::Math::XYZVectorF decayVtx, decayVtxLambda;
  lVector lproL, lpiL;
  std::unordered_map<int, int> mcMap;
  auto fillTreeInfo = [&](std::array<AliESDtrack *, 3> tracks, std::array<float, 3> nSigmaTPC, std::array<float, 3> nSigmaTOF) {
    const float mass = lambdab.mass();
    if (mass < fMassWindow[0] || mass > fMassWindow[1])
      return false;

    const float totalMom = lambdab.P();
    const float len = std::sqrt(decayVtx.Mag2());
    fRecHyp->cosPA = lambdab.Vect().Dot(decayVtx) / (totalMom * len);
    const float cosPA = fUseAbsCosPAcut ? std::abs(fRecHyp->cosPA) : fRecHyp->cosPA;
    fRecHyp->ct = len * kLambdaBMass / totalMom;
    if (fRecHyp->ct < fCandidateCtRange[0] || fRecHyp->ct > fCandidateCtRange[1])
      return false;
    if (fCosPAspline)
    {
      if (cosPA < fCosPAspline->Eval(fRecHyp->ct))
        return false;
    }
    else if (cosPA < fMinCosPA)
    {
      return false;
    }
    fRecHyp->r = decayVtx.Rho();
    fRecHyp->positive = tracks[0]->Charge() > 0;
    fRecHyp->pt = lambdab.pt();
    fRecHyp->phi = lambdab.phi();
    fRecHyp->pz = lambdab.pz();
    fRecHyp->m = mass;

    float dca[2], bCov[3];
    tracks[0]->GetImpactParameters(dca, bCov);
    fRecHyp->dca_de = std::hypot(dca[0], dca[1]);
    tracks[1]->GetImpactParameters(dca, bCov);
    fRecHyp->dca_pr = std::hypot(dca[0], dca[1]);
    tracks[2]->GetImpactParameters(dca, bCov);
    fRecHyp->dca_pi = std::hypot(dca[0], dca[1]);

    fRecHyp->hasTOF_de = HasTOF(tracks[0]);
    fRecHyp->hasTOF_pr = HasTOF(tracks[1]);
    fRecHyp->hasTOF_pi = HasTOF(tracks[2]);

    fRecHyp->tofNsig_de = nSigmaTOF[0];
    fRecHyp->tofNsig_pr = nSigmaTOF[1];
    fRecHyp->tofNsig_pi = nSigmaTOF[2];

    fRecHyp->tpcNsig_de = nSigmaTPC[0];
    fRecHyp->tpcNsig_pr = nSigmaTPC[1];
    fRecHyp->tpcNsig_pi = nSigmaTPC[2];

    fRecHyp->tpcClus_de = tracks[0]->GetTPCsignalN();
    fRecHyp->tpcClus_pr = tracks[1]->GetTPCsignalN();
    fRecHyp->tpcClus_pi = tracks[2]->GetTPCsignalN();

    fRecHyp->its_clusmap_de = tracks[0]->GetITSClusterMap();
    fRecHyp->its_clusmap_pr = tracks[1]->GetITSClusterMap();
    fRecHyp->its_clusmap_pi = tracks[2]->GetITSClusterMap();

    fRecHyp->is_ITSrefit_de = tracks[0]->GetStatus() & AliVTrack::kITSrefit;
    fRecHyp->is_ITSrefit_pr = tracks[1]->GetStatus() & AliVTrack::kITSrefit;
    fRecHyp->is_ITSrefit_pi = tracks[2]->GetStatus() & AliVTrack::kITSrefit;

    return true;
  };

  fVertexer.setBz(esdEvent->GetMagneticField());
  fVertexerLambda.setBz(esdEvent->GetMagneticField());
  int indices[2][3]{{1, 1, 0}, {0, 0, 1}};

  //!!
  RHyperTriton3O2 &o2RecHyp = *(RHyperTriton3O2 *)fRecHyp;

  for (int idx{0}; idx < 2; ++idx)
  {
    for (const auto &he3 : helpers[kDeuteron][indices[idx][0]])
    {
      int rotations{0};
      auto He3TrackSnapshot = *he3.track;
      double alpha = he3.track->GetAlpha();
      do
      {
        o2RecHyp.rotation = rotations;
        if (rotations)
        {
          double deltaAngle{rotations * TMath::TwoPi() / (fTrackRotations + 1)};
          He3TrackSnapshot.SetParamOnly(he3.track->GetX(), alpha + deltaAngle, he3.track->getParams());
        }
        for (const auto &Ap1 : helpers[kAProton][indices[idx][1]])
        {
          if (he3.track == Ap1.track)
            continue;

          for (const auto &Ap2 : helpers[kAProton][indices[idx][2]])
          {
            if (Ap1.track == Ap2.track || he3.track == Ap2.track || he3.track == Ap1.track)
              continue;

            ROOT::Math::SVector<double, 3U> vert;
            lVector lhe3, lApro1, lApro2;
            int nVert{0};

            try
            {
              nVert = fVertexer.process(He3TrackSnapshot, *Ap1.track, *Ap2.track);
            }
            catch (std::runtime_error &e)
            {
            }
            if (!nVert)
              continue;

            fVertexer.propagateTracksToVertex();
            auto &he3Track = fVertexer.getTrack(0);
            auto &Apr1Track = fVertexer.getTrack(1);
            auto &Apr2Track = fVertexer.getTrack(2);
            lhe3.SetCoordinates((float)he3Track.Pt(), (float)he3Track.Eta(), (float)he3Track.Phi(), kDeuMass);
            lApro1.SetCoordinates((float)Apr1Track.Pt(), (float)Apr1Track.Eta(), (float)Apr1Track.Phi(), kPMass);
            lApro2.SetCoordinates((float)Apr2Track.Pt(), (float)Apr2Track.Eta(), (float)Apr2Track.Phi(), kPiMass);

            lambdab = lhe3 + lApro1 + lApro2;

            o2RecHyp.mppi = (lApro1 + lApro2).mass2();
            o2RecHyp.mdpi = (lhe3 + lApro2).mass2();
            ROOT::Math::Boost boostHyper{lambdab.BoostToCM()};
            auto d{boostHyper(lhe3).Vect()};
            // auto lambda{boostHyper(lApro1 + lApro2).Vect()};
            auto pV{boostHyper(lApro1).Vect()};
            auto piV{boostHyper(lApro2).Vect()};
            o2RecHyp.momDstar = std::sqrt(d.Mag2());
            o2RecHyp.cosThetaStar = d.Dot(lambdab.Vect()) / (o2RecHyp.momDstar * lambdab.P());
            o2RecHyp.cosTheta_ProtonPiH = pV.Dot(piV) / std::sqrt(pV.Mag2() * piV.Mag2());
            vert = fVertexer.getPCACandidate();
            decayVtx.SetCoordinates((float)(vert[0] - pvPos[0]), (float)(vert[1] - pvPos[1]), (float)(vert[2] - pvPos[2]));
            o2RecHyp.candidates = nVert;

            double He3Pos[3], Apro1Pos[3], Apro2Pos[3];
            he3Track.GetXYZ(He3Pos);
            Apr1Track.GetXYZ(Apro1Pos);
            Apr2Track.GetXYZ(Apro2Pos);

            o2RecHyp.dca_He3_Apr1 = Hypot(He3Pos[0] - Apro1Pos[0], He3Pos[1] - Apro1Pos[1], He3Pos[2] - Apro1Pos[2]);
            if (o2RecHyp.dca_He3_Apr1 > fMaxTrack2TrackDCA[0])
              continue;
            o2RecHyp.dca_He3_Apr2 = Hypot(He3Pos[0] - Apro2Pos[0], He3Pos[1] - Apro2Pos[1], He3Pos[2] - Apro2Pos[2]);
            if (o2RecHyp.dca_He3_Apr2 > fMaxTrack2TrackDCA[1])
              continue;

            o2RecHyp.dca_Apr1_Apr2 = Hypot(Apro1Pos[0] - Apro2Pos[0], Apro1Pos[1] - Apro2Pos[1], Apro1Pos[2] - Apro2Pos[2]);
            if (o2RecHyp.dca_Apr1_Apr2 > fMaxTrack2TrackDCA[2])
              continue;

            o2RecHyp.dca_He3_sv = Hypot(He3Pos[0] - vert[0], He3Pos[1] - vert[1], He3Pos[2] - vert[2]);
            if (o2RecHyp.dca_He3_sv > fMaxTrack2SVDCA[0])
              continue;
            o2RecHyp.dca_Apr1_sv = Hypot(Apro1Pos[0] - vert[0], Apro1Pos[1] - vert[1], Apro1Pos[2] - vert[2]);
            if (o2RecHyp.dca_Apr1_sv > fMaxTrack2SVDCA[1])
              continue;
            o2RecHyp.dca_Apr2_sv = Hypot(Apro2Pos[0] - vert[0], Apro2Pos[1] - vert[1], Apro2Pos[2] - vert[2]);
            if (o2RecHyp.dca_Apr2_sv > fMaxTrack2SVDCA[2])
              continue;

            o2RecHyp.chi2 = fVertexer.getChi2AtPCACandidate();

            std::array<AliESDtrack *, 3> tracks{(AliESDtrack *)he3.track, (AliESDtrack *)Apr1.track, (AliESDtrack *)Apr2.track};
            std::array<float, 3> nSigmaTPC{he3.nSigmaTPC, Apr1.nSigmaTPC, Apr2.nSigmaTPC};
            std::array<float, 3> nSigmaTOF{he3.nSigmaTOF, Apr1.nSigmaTOF, Apr2.nSigmaTOF};
            if (!fillTreeInfo(tracks, nSigmaTPC, nSigmaTOF))
              continue;

            bool record{!fMC || !fOnlyTrueCandidates};
            if (fMC)
            {
              int momId = IsTrueLambdaBCandidate((AliESDtrack *)he3.track, (AliESDtrack *)p.track, (AliESDtrack *)pi.track, mcEvent);
              record = record || momId >= 0;
              if (record)
              {
                FillGenHypertriton(fGenHypO2, momId, true, mcEvent);
                mcMap[momId] = 1;
              }
            }
            if (record)
            {
              fTreeHyp3->Fill();
            }
          }
        }
      } while (rotations++ < fTrackRotations);
    }
  }
  if (fEnableEventMixing)
    FillEventMixingPool(fEventCuts.GetCentrality(), pvPos[2], heliumsForMixing);

  if (fMC)
  {
    RHyperTriton3O2 rec;
    rec.centrality = fRecHyp->centrality;
    rec.trigger = fRecHyp->trigger;
    *fRecHyp = rec;
    for (int iTrack = 0; iTrack < mcEvent->GetNumberOfTracks(); iTrack++)
    {
      AliVParticle *part = mcEvent->GetTrack(iTrack);
      if (!part)
      {
        ::Warning("AliAnalysisTaskLambdaB::UserExec",
                  "Generated loop %d - MC TParticle pointer to current stack particle = 0x0 ! Skipping.", iTrack);
        continue;
      }
      if (std::abs(part->Y()) > 1.)
        continue;
      if (!IsLambdaB(part, mcEvent))
        continue;
      if (mcMap.find(iTrack) != mcMap.end())
        continue;
      FillGenHypertriton(fGenHypO2, iTrack, false, mcEvent);
      fTreeHyp3->Fill();
    }
  }

  PostData(1, fListHist);
  PostData(2, fTreeHyp3);
}

int AliAnalysisTaskLambdaB::FindEventMixingCentBin(const float centrality)
{
  if (centrality >= 100.)
    return -999;
  return static_cast<int>(centrality / 10);
}

int AliAnalysisTaskLambdaB::FindEventMixingZBin(const float zvtx)
{
  if (zvtx > 10. || zvtx < -10.)
    return -999.;
  return static_cast<int>((zvtx + 10.) / 2);
}

int AliAnalysisTaskLambdaB::CheckPionCharge(std::vector<AliESDtrack *> tracks[2][2], AliESDv0 v0)
{

  double pP[3], nP[3];
  v0.GetPPxPyPz(pP[0], pP[1], pP[2]);
  v0.GetNPxPyPz(nP[0], nP[1], nP[2]);
  int isPiPositive = -1;
  if (tracks[1][1].size() <= int(v0.GetPindex()) && tracks[0][1].size() <= int(v0.GetNindex()))
    return isPiPositive;
  isPiPositive = (tracks[1][1].size() <= int(v0.GetPindex())) ? 0 : -1;
  isPiPositive = (tracks[0][1].size() <= int(v0.GetNindex())) ? 1 : -1;

  if (isPiPositive == -1)
  {
    double posDiff = std::abs(tracks[1][1][v0.GetPindex()]->Px() - pP[0]);
    double negDiff = std::abs(tracks[0][1][v0.GetNindex()]->Px() - nP[0]);
    isPiPositive = posDiff < negDiff ? 1 : 0;
  }
  return isPiPositive;
}

void AliAnalysisTaskLambdaB::FillEventMixingPool(const float centrality, const float zvtx,
                                                      const std::vector<EventMixingTrack> &tracks)
{
  int centBin = FindEventMixingCentBin(centrality);
  int zBin = FindEventMixingZBin(zvtx);

  auto &trackVector = fEventMixingPool[centBin][zBin];

  for (auto &t : tracks)
    trackVector.emplace_back(t);

  while (trackVector.size() > fEventMixingPoolDepth)
    trackVector.pop_front();

  return;
}

std::vector<EventMixingTrack *> AliAnalysisTaskLambdaB::GetEventMixingTracks(const float centrality,
                                                                                  const float zvtx)
{
  int centBin = FindEventMixingCentBin(centrality);
  int zBin = FindEventMixingZBin(zvtx);

  std::vector<EventMixingTrack *> tmpVector;

  for (auto &v : fEventMixingPool[centBin][zBin])
  {
    if (v.used >= fEventMixingPoolMaxReuse)
      continue;
    tmpVector.emplace_back(&(v));
    v.used++;
  }

  return tmpVector;
}

AliAnalysisTaskLambdaB *AliAnalysisTaskLambdaB::AddTask(bool isMC, TString suffix)
{
  // Get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskLambdaB2BodyML", "No analysis manager found.");
    return nullptr;
  }

  // Check the analysis type using the event handlers connected to the analysis
  // manager.
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskLambdaB", "This task requires an input event handler");
    return nullptr;
  }

  TString tskname = "AliAnalysisTaskLambdaB";
  tskname.Append(suffix.Data());
  AliAnalysisTaskLambdaB *task = new AliAnalysisTaskLambdaB(isMC, tskname.Data());

  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(
      Form("%s_summary", tskname.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, "AnalysisResults.root");

  AliAnalysisDataContainer *coutput2 =
      mgr->CreateContainer(Form("HyperTritonTree%s", suffix.Data()), TTree::Class(),
                           AliAnalysisManager::kOutputContainer, Form("HyperTritonTree3.root:%s", suffix.Data()));
  coutput2->SetSpecialOutput();

  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, coutput1);
  mgr->ConnectOutput(task, 2, coutput2);
  return task;
}
