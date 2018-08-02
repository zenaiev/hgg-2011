// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// >>>>>>>>>>>> Helper for ttbar event selection >>>>>>>>>>>>>>>>>>>>>>>
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// Consult analysis documentation (papers, description-ttbar.pdf) for 
// better description of applied cuts etc.

// additional files from this analysis 
#include "tree.h"
// C++ library or ROOT header files
#include <TLorentzVector.h>

bool gFlagDebug = 0;

// photon classes
int PhotonClass(const double eta, const double r9)
{
  if(TMath::Abs(eta) < 1.44 && r9 > 0.94)
    return 3;
  if(TMath::Abs(eta) < 1.44 && r9 < 0.94)
    return 4;
  if(TMath::Abs(eta) > 1.57 && TMath::Abs(eta) < 2.5 && r9 > 0.94)
    return 5;
  if(TMath::Abs(eta) > 1.57 && TMath::Abs(eta) < 2.5 && r9 < 0.94)
    return 6;

  /*if(TMath::Abs(eta) <= 1.44)
    return 1;
  if(TMath::Abs(eta) >= 1.57 && TMath::Abs(eta) <= 2.5)
    return 2;*/

  printf("Error: cannot determine photon in barrel or endcap, return 0\n");
  return 0;
}

// Routine for muon selection
// Arguments:
//   const ZTree* preselTree: input tree (see tree.h), GetEntry() should be done already
//   const int el: muon candidate
// Returns true for selected muon, false otherwise.
// See tree.h for ZTree variables description.
double SelectPh11(const int eventClass, const ZTree* preselTree, const int ph)
{
  
  // relative combined isolation using selected event vertex 5.4.1
  if(gFlagDebug) printf("5.4.1\n");
  const double aEff = 0.17;
  double relCombIso = preselTree->phTrkSumPtHollowConeDR03[ph] + preselTree->phEcalRecHitSumEtConeDR03[ph] + preselTree->phHcalTowerSumEtConeDR04[ph];
  //double relCombIso = preselTree->phTrkSumPtHollowConeDR03[ph] + preselTree->phEcalRecHitSumEtConeDR03[ph] + preselTree->phHcalTowerSumEtConeDR03[ph];
  relCombIso -= aEff * preselTree->rho;
  relCombIso /= (TMath::Abs(preselTree->phPt[ph]) / 50.0);
  if( (relCombIso > 3.8 && eventClass == 3) ||
      (relCombIso > 2.2 && eventClass == 4) ||
      (relCombIso > 1.77 && eventClass == 5) ||
      (relCombIso > 1.29 && eventClass == 6) )
    return 0;

  // relative track isolation using selected event vertex 5.4.3
  if(gFlagDebug) printf("5.4.3\n");
  double relTrackIso = preselTree->phTrkSumPtHollowConeDR03[ph];
  relTrackIso /= (TMath::Abs(preselTree->phPt[ph]) / 50.0);
  if( (relTrackIso > 3.5 && eventClass == 3) ||
      (relTrackIso > 2.2 && eventClass == 4) ||
      (relTrackIso > 2.3 && eventClass == 5) ||
      (relTrackIso > 1.45 && eventClass == 6) )
    return 0;

  // H/E 5.4.4
  if(gFlagDebug) printf("5.4.4\n");
  if( (preselTree->phHadronicOverEm[ph] > 0.082 && eventClass == 3) ||
      (preselTree->phHadronicOverEm[ph] > 0.062 && eventClass == 4) ||
      (preselTree->phHadronicOverEm[ph] > 0.065 && eventClass == 5) ||
      (preselTree->phHadronicOverEm[ph] > 0.048 && eventClass == 6) )
    return 0;

  // covietaieta 5.4.5
  if(gFlagDebug) printf("5.4.5\n");
  double sigmaIEtaIEta = preselTree->phSigmaIetaIeta[ph] * preselTree->phSigmaIetaIeta[ph];
  if( (sigmaIEtaIEta > 0.0106 && eventClass == 3) ||
      (sigmaIEtaIEta > 0.0097 && eventClass == 4) ||
      (sigmaIEtaIEta > 0.028 && eventClass == 5) ||
      (sigmaIEtaIEta > 0.027 && eventClass == 6) )
    return 0;

  // r9 5.4.6
  if(gFlagDebug) printf("5.4.6\n");
  if( (preselTree->phR9[ph] < 0.94 && eventClass == 3) ||
      (preselTree->phR9[ph] < 0.36 && eventClass == 4) ||
      (preselTree->phR9[ph] < 0.94 && eventClass == 5) ||
      (preselTree->phR9[ph] < 0.32 && eventClass == 6) )
    return 0;
  if(gFlagDebug) printf("PASSED\n");
  
  //electron veto (if eventclass 4 is not rejected, number increases because deltaR is too small)
  if( preselTree->phNumElectronsSuperCluster[ph] > 0)
    return 0;
    
  
  // matching
  if(preselTree->_flagMC && preselTree->phMatch[ph] > 0.1)
    return 0;
  
  // all cuts passed: return true
  return 1;
}

double SelectPh12(const int eventClass, const ZTree* preselTree, const int ph)
{
  //event preselection on the photon
  double EtCorrEcalIso = preselTree->phEcalRecHitSumEtConeDR03[ph] - 0.012 * preselTree->phPt[ph];
  double EtCorrHcalIso = preselTree->phHcalTowerSumEtConeDR03[ph] - 0.005 * preselTree->phPt[ph];
  double EtCorrTrkIso = preselTree->phTrkSumPtHollowConeDR03[ph] - 0.002 * preselTree->phPt[ph];
  double ChargedPFIso = preselTree->phChargedHadronIso[ph];
  //R9 <= 0.9
  if(preselTree->phR9[ph] <= 0.9)
  {
    //for both barrel and endcap
    if(EtCorrEcalIso > 4 || EtCorrHcalIso > 4 || EtCorrTrkIso > 4  || ChargedPFIso > 4)
      return 0;
    //for barrel
    if((eventClass == 3 || eventClass == 4)  && (preselTree->phHadronicOverEm[ph] > 0.075 || (preselTree->phSigmaIetaIeta[ph]*preselTree->phSigmaIetaIeta[ph])  > 0.014) )
      return 0;
    //for endcap
    if((eventClass == 5 || eventClass == 6) && (preselTree->phHadronicOverEm[ph] > 0.075 || (preselTree->phSigmaIetaIeta[ph]*preselTree->phSigmaIetaIeta[ph])  > 0.034))
      return 0;
  }
  //R9 > 0.9
  if(preselTree->phR9[ph] > 0.9)
  {
    //for both barrel and endcap
    if(EtCorrEcalIso > 50 ||  EtCorrHCalIso > 50 || EtCorrTrkIso > 50 || ChargedPFIso > 4)
      return 0;
    //for barrel
    if((eventClass == 3 || eventClass == 4)  && (preselTree->phHadronicOverEm[ph] > 0.082 || (preselTree->phSigmaIetaIeta[ph]*preselTree->phSigmaIetaIeta[ph]) > 0.014) )
      return 0;
    //for endcap
    if((eventClass == 5 || eventClass == 6) && (preselTree->phHadronicOverEm[ph] > 0.075 ||  (preselTree->phSigmaIetaIeta[ph]*preselTree->phSigmaIetaIeta[ph]) > 0.034))
      return 0;
  }
  //PFlow isolation 
  if(gFlagDebug) printf("PFlow isolation sum\n");
  double isoSum = preselTree->phIsolationSum[ph]; //add pile-up correction
  if( (isoSum > 6 && eventClass == 3) ||
      (isoSum > 4.7 && eventClass == 4) ||
      (isoSum > 5.6 && eventClass == 5) ||
      (isoSum > 3.6 && eventClass == 6) )
    return 0;
  //PFlow charged hadron isolation
  if(gFlagDebug) printf("PFlow charged hadron iso\n");
  double isoCharged = preselTree->phChargedHadronIso[ph]; //add pile-up correction
  if( (isoCharged > 3.8 && eventClass == 3) ||
      (isoCharged > 2.5 && eventClass == 4) ||
      (isoCharged > 3.1 && eventClass == 5) ||
      (isoCharged > 2.2 && eventClass == 6) )
    return 0;
  //PFlow isolation worst vertex
  double isoSumWrongVtx = preselTree->phIsolationSumWrongVtx[ph]; //add pile-up correction
  if(gFlagDebug) printf("Pflow photon iso wrong vertex\n");
  if( (isoSumWrongVtx > 10 && eventClass == 3) ||
      (isoSumWrongVtx > 6.5 && eventClass == 4) ||
      (isoSumWrongVtx > 5.6 && eventClass == 5) ||
      (isoSumWrongVtx > 4.4 && eventClass == 6) )
    return 0;
  
  // H/E 5.4.4
  if(gFlagDebug) printf("5.4.4\n");
  if( (preselTree->phHadronicOverEm[ph] > 0.124 && eventClass == 3) ||
      (preselTree->phHadronicOverEm[ph] > 0.092 && eventClass == 4) ||
      (preselTree->phHadronicOverEm[ph] > 0.142 && eventClass == 5) ||
      (preselTree->phHadronicOverEm[ph] > 0.063 && eventClass == 6) )
    return 0;

  // covietaieta 5.4.5
  if(gFlagDebug) printf("5.4.5\n");
  double sigmaIEtaIEta = preselTree->phSigmaIetaIeta[ph] * preselTree->phSigmaIetaIeta[ph];
  if( (sigmaIEtaIEta > 0.0108 && eventClass == 3) ||
      (sigmaIEtaIEta > 0.0102 && eventClass == 4) ||
      (sigmaIEtaIEta > 0.028 && eventClass == 5) ||
      (sigmaIEtaIEta > 0.028 && eventClass == 6) )
    return 0;

  // r9 5.4.6
  if(gFlagDebug) printf("5.4.6\n");
  if( (preselTree->phR9[ph] < 0.94 && eventClass == 3) ||
      (preselTree->phR9[ph] < 0.298 && eventClass == 4) ||
      (preselTree->phR9[ph] < 0.94 && eventClass == 5) ||
      (preselTree->phR9[ph] < 0.24 && eventClass == 6) )
    return 0;
  if(gFlagDebug) printf("PASSED\n");
  
  //electron veto
  if(gFlagDebug) printf("electron veto \n");
  if( preselTree->phNumElectronsSuperCluster[ph] > 0)
    return 0;
  //matching
  if(preselTree->_flagMC && preselTree->phMatch[ph] > 0.1)
    return 0;
  
  // all cuts passed: return true
  return 1;
}

// routine for electron-muon pair selection
// (select best e-mu pair in the event, with highest pT)
// Arguments:
//   const ZTree* preselTree: input tree (see tree.h), GetEntry() should be done already
//   TLorentzVector& vecLepM: selected lepton- (output)
//   TLorentzVector& vecLepP: selected lepton+ (output)
//   double& maxPtDiLep: transverse momentum of the selected dilepton pair (output)
// If no dilepton pair is selected, maxPtDiLep remains unchanged 
int SelectHgg(const ZTree* preselTree, int reqEventClass, TLorentzVector& momPh1, TLorentzVector& momPh2, TString name)
{
  int phClass[2];
  double phR9[2];
  TLorentzVector higgs, ph[2];
  double maxPtSum = -1.0;
  int eventClass = -1;

  // double loop over photons
  for(int ph1 = 0; ph1 < preselTree->Nph; ph1++)
  {
	//photon position 
    if(TMath::Abs(preselTree->phEta[ph1]) > 2.5)
      return 0;
    if(TMath::Abs(preselTree->phEta[ph1]) > 1.44 && TMath::Abs(preselTree->phEta[ph1]) < 1.57)
      return 0;
    //get Photon class
    phClass[0] = PhotonClass(preselTree->phEta[ph1], preselTree->phR9[ph1]);
    //PhotonID selection
    if(name == "data2011")
    {
      if(!SelectPh11(phClass[0], preselTree, ph1))
        continue;
    }
    if(name == "data2012")
    {
      if(!SelectPh12(phClass[0], preselTree, ph1))
        continue;
    }
    ph[0].SetPtEtaPhiM(preselTree->phPt[ph1], preselTree->phEta[ph1], preselTree->phPhi[ph1], 0.0);

    for(int ph2 = ph1 + 1; ph2 < preselTree->Nph; ph2++)
    {
	  //photon position
      if(TMath::Abs(preselTree->phEta[ph2]) > 2.5)
        return 0;
      if(TMath::Abs(preselTree->phEta[ph2]) > 1.44 && TMath::Abs(preselTree->phEta[ph2]) < 1.57)
        return 0;
      //get photon class
      phClass[1] = PhotonClass(preselTree->phEta[ph2], preselTree->phR9[ph2]);
      //PhotonID selection
      if(name == "data2011")
      {
        if(!SelectPh11(phClass[1], preselTree, ph2))
          continue;
      }
      if(name == "data2012")
      {
        if(!SelectPh12(phClass[1], preselTree, ph2))
          continue;
      }
      ph[1].SetPtEtaPhiM(preselTree->phPt[ph2], preselTree->phEta[ph2], preselTree->phPhi[ph2], 0.0);

      // determine event class
      double r9Min = std::min(preselTree->phR9[ph1], preselTree->phR9[ph2]);
      bool flagBothInBarrel = ((phClass[0] == 3 || phClass[0] == 4) && (phClass[1] == 3 || phClass[1] == 4));
      if(flagBothInBarrel)
      {
        if(r9Min > 0.94)
          eventClass = 3;
        else
          eventClass = 4;
      }
      else
      {
        if(r9Min > 0.94)
          eventClass = 5;
        else
          eventClass = 6;
      }


      if(reqEventClass != 1 && reqEventClass != eventClass)
        continue;

      // calculate higgs invariant mass and requirements on photon pT
      higgs = ph[0] + ph[1];
      double mgg = higgs.M();
      if(ph[0].Pt() > ph[1].Pt())
      {
        if(ph[0].Pt() < (mgg / 3.0) || ph[1].Pt() < (mgg / 4.0))
          continue;
      }
      else
      {
        if(ph[0].Pt() < (mgg / 4.0) || ph[1].Pt() < (mgg / 3.0))
          continue;
      }

      double sumPt = ph[0].Pt() + ph[1].Pt();
      if(sumPt < maxPtSum)
        continue;
      maxPtSum = sumPt;
      momPh1 = ph[0];
      momPh2 = ph[1];
    }
  }

  if(maxPtSum < 0.0)
    return 0;

  return eventClass;
}
