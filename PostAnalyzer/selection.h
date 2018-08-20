// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// >>>>>>>>>>>> Helper for gamma gamma event selection >>>>>>>>>>>>>>>>>>>>>>>
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// Consult analysis documentation (papers, description-gammagamma.pdf) for
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

  printf("Error: cannot determine photon in barrel or endcap, return 0\n");
  return 0;
}

// Routine for photon 2011 selection
// Arguments:
//	 const int eventClass	: EventClass of the event
//   const ZTree* preselTree: input tree (see tree.h), GetEntry() should be done already
//   const int ph			: photon candidate
// Returns true for selected photon, false otherwise.
// See tree.h for ZTree variables description.
double SelectPh11(const int phClass, const ZTree* preselTree, const int ph)
{
  // relative combined isolation using selected event vertex 5.4.1
  if(gFlagDebug) printf("5.4.1\n");
  const double aEff = 0.17;
  double relCombIso = preselTree->phTrkSumPtHollowConeDR03[ph] + preselTree->phEcalRecHitSumEtConeDR03[ph] + preselTree->phHcalTowerSumEtConeDR04[ph];
  relCombIso -= aEff * preselTree->rho;
  relCombIso /= (TMath::Abs(preselTree->phPt[ph]) / 50.0);
  if( (relCombIso > 3.8 && phClass == 3) ||
      (relCombIso > 2.2 && phClass == 4) ||
      (relCombIso > 1.77 && phClass == 5) ||
      (relCombIso > 1.29 && phClass == 6) )
    return 0;

  // relative track isolation using selected event vertex 5.4.3
  if(gFlagDebug) printf("5.4.3\n");
  double relTrackIso = preselTree->phTrkSumPtHollowConeDR03[ph];
  relTrackIso /= (TMath::Abs(preselTree->phPt[ph]) / 50.0);
  if( (relTrackIso > 3.5 && phClass == 3) ||
      (relTrackIso > 2.2 && phClass == 4) ||
      (relTrackIso > 2.3 && phClass == 5) ||
      (relTrackIso > 1.45 && phClass == 6) )
    return 0;

  // H/E 5.4.4
  if(gFlagDebug) printf("5.4.4\n");
  if( (preselTree->phHadronicOverEm[ph] > 0.082 && phClass == 3) ||
      (preselTree->phHadronicOverEm[ph] > 0.062 && phClass == 4) ||
      (preselTree->phHadronicOverEm[ph] > 0.065 && phClass == 5) ||
      (preselTree->phHadronicOverEm[ph] > 0.048 && phClass == 6) )
    return 0;

  // sigmaietaieta 5.4.5
  if(gFlagDebug) printf("5.4.5\n");
  double sigmaIEtaIEta = preselTree->phSigmaIetaIeta[ph] * preselTree->phSigmaIetaIeta[ph];
  if( (sigmaIEtaIEta > 0.0106 && phClass == 3) ||
      (sigmaIEtaIEta > 0.0097 && phClass == 4) ||
      (sigmaIEtaIEta > 0.028 && phClass == 5) ||
      (sigmaIEtaIEta > 0.027 && phClass == 6) )
    return 0;

  // r9 5.4.6
  if(gFlagDebug) printf("5.4.6\n");
  if( (preselTree->phR9[ph] < 0.94 && phClass == 3) ||
      (preselTree->phR9[ph] < 0.36 && phClass == 4) ||
      (preselTree->phR9[ph] < 0.94 && phClass == 5) ||
      (preselTree->phR9[ph] < 0.32 && phClass == 6) )
    return 0;
  if(gFlagDebug) printf("PASSED\n");
  
  //electron veto
  if( preselTree->phNumElectronsSuperCluster[ph] > 0)
    return 0;

  // matching
  if(preselTree->_flagMC && preselTree->phMatch[ph] > 0.1)
    return 0;
  
  // all cuts passed: return true
  return 1;
}


// Routine for photon 2012 selection
// Arguments:
//	 const int eventClass	: EventClass of the event
//   const ZTree* preselTree: input tree (see tree.h), GetEntry() should be done already
//   const int ph			: photon candidate
// Returns true for selected photon, false otherwise.
// See tree.h for ZTree variables description.
double SelectPh12(const int phClass, const ZTree* preselTree, const int ph)
{
  //trigger check
  if(preselTree->Triggers == 0)
  {
    printf("*******No Trigger fired!!!*******");
    return 0;
  }
  //effective area for pile up (could not find exact determined A_eff in analysis note/paper) 
  const double aEff = 0.17; //from 2011 -> 2012 should be higher
  //event preselection on the photon
  double EtCorrHcalIso = preselTree->phHcalTowerSumEtConeDR04[ph] - 0.005 * preselTree->phPt[ph];
  double EtCorrTrkIso = preselTree->phTrkSumPtHollowConeDR04[ph] - 0.002 * preselTree->phPt[ph];
  double ChargedPFIso = preselTree->phChargedHadronIsoDR02[ph];
  //R9 <= 0.9
  if(preselTree->phR9[ph] <= 0.9)
  {
    //for both barrel and endcap
    if(EtCorrHcalIso > 4 || EtCorrTrkIso > 4  || ChargedPFIso > 4)
      return 0;
    //for barrel
    if((phClass == 3 || phClass == 4)  && (preselTree->phHadronicOverEm[ph] > 0.075 || (preselTree->phSigmaIetaIeta[ph]*preselTree->phSigmaIetaIeta[ph])  > 0.014) )
      return 0;
    //for endcap
    if((phClass == 5 || phClass == 6) && (preselTree->phHadronicOverEm[ph] > 0.075 || (preselTree->phSigmaIetaIeta[ph]*preselTree->phSigmaIetaIeta[ph])  > 0.034))
      return 0;
  }
  //R9 > 0.9
  if(preselTree->phR9[ph] > 0.9)
  {
    //for both barrel and endcap
    if(EtCorrHcalIso > 50 || EtCorrTrkIso > 50 || ChargedPFIso > 4)
      return 0;
    //for barrel
    if((phClass == 3 || phClass == 4)  && (preselTree->phHadronicOverEm[ph] > 0.082 || (preselTree->phSigmaIetaIeta[ph]*preselTree->phSigmaIetaIeta[ph]) > 0.014) )
      return 0;
    //for endcap
    if((phClass == 5 || phClass == 6) && (preselTree->phHadronicOverEm[ph] > 0.075 ||  (preselTree->phSigmaIetaIeta[ph]*preselTree->phSigmaIetaIeta[ph]) > 0.034))
      return 0;
  }
  //PFlow isolation 
  if(gFlagDebug) printf("PFlow isolation sum\n");
  double isoSum = preselTree->phIsolationSumDR04[ph] - aEff * preselTree->rho; 
  if( (isoSum > 6 && phClass == 3) ||
      (isoSum > 4.7 && phClass == 4) ||
      (isoSum > 5.6 && phClass == 5) ||
      (isoSum > 3.6 && phClass == 6) )
    return 0;
  //PFlow charged hadron isolation
  if(gFlagDebug) printf("PFlow charged hadron iso\n");
  double isoCharged = preselTree->phChargedHadronIsoDR04[ph] - aEff * preselTree->rho;
  if( (isoCharged > 3.8 && phClass == 3) ||
      (isoCharged > 2.5 && phClass == 4) ||
      (isoCharged > 3.1 && phClass == 5) ||
      (isoCharged > 2.2 && phClass == 6) )
    return 0;
  //PFlow isolation worst vertex
  double isoSumWrongVtx = preselTree->phIsolationSumWrongVtxDR04[ph] - aEff * preselTree->rho;
  if(gFlagDebug) printf("Pflow photon iso wrong vertex\n");
  if( (isoSumWrongVtx > 10 && phClass == 3) ||
      (isoSumWrongVtx > 6.5 && phClass == 4) ||
      (isoSumWrongVtx > 5.6 && phClass == 5) ||
      (isoSumWrongVtx > 4.4 && phClass == 6) )
    return 0;
  
  // H/E 5.4.4
  if(gFlagDebug) printf("5.4.4\n");
  if( (preselTree->phHadronicOverEm[ph] > 0.124 && phClass == 3) ||
      (preselTree->phHadronicOverEm[ph] > 0.094 && phClass == 4) ||
      (preselTree->phHadronicOverEm[ph] > 0.142 && phClass == 5) ||
      (preselTree->phHadronicOverEm[ph] > 0.063 && phClass == 6) )
    return 0;

  // covietaieta 5.4.5
  if(gFlagDebug) printf("5.4.5\n");
  double sigmaIEtaIEta = preselTree->phSigmaIetaIeta[ph] * preselTree->phSigmaIetaIeta[ph];
  if( (sigmaIEtaIEta > 0.0108 && phClass == 3) ||
      (sigmaIEtaIEta > 0.0102 && phClass == 4) ||
      (sigmaIEtaIEta > 0.028 && phClass == 5) ||
      (sigmaIEtaIEta > 0.028 && phClass == 6) )
    return 0;

  // r9 5.4.6
  if(gFlagDebug) printf("5.4.6\n");
  if( (preselTree->phR9[ph] < 0.94 && phClass == 3) ||
      (preselTree->phR9[ph] < 0.298 && phClass == 4) ||
      (preselTree->phR9[ph] < 0.94 && phClass == 5) ||
      (preselTree->phR9[ph] < 0.24 && phClass == 6) )
    return 0;
  if(gFlagDebug) printf("PASSED\n");
  
  //electron veto
  if(gFlagDebug) printf("electron veto \n");
  if( preselTree->phNumElectronsSuperCluster[ph] > 0)
    return 0;

  // all cuts passed: return true
  return 1;
}

// routine for photon-photon pair selection
// (select best ph-ph pair in the event, with highest pT)
// Arguments:
//   const ZTree* preselTree: input tree (see tree.h), GetEntry() should be done already
//	 int reqEventClass 		: for debugging purposes to only get one class
//   TLorentzVector& momPh1	: selected photon (output)
//   TLorentzVector& momPh2	: selected photon (output)
//   TSring name			: Name of the ntuple file to be analyzed
// Returns event class of chosen gamma-gamma event
int SelectHgg(const ZTree* preselTree, int reqEventClass, TLorentzVector& momPh1, TLorentzVector& momPh2, TString name)
{
  int phClass[2];
  double phR9[2];
  TLorentzVector higgs, ph[2];
  double maxPtSum = -1.0;
  int eventClass = -1;

  //loop over photons
  for(int ph1 = 0; ph1 < preselTree->Nph; ph1++)
  {
	//photon position 
    if(TMath::Abs(preselTree->phEta[ph1]) > 2.5)
      return 0;
    if(TMath::Abs(preselTree->phEta[ph1]) > 1.44 && TMath::Abs(preselTree->phEta[ph1]) < 1.57)
      return 0;

    //get Photon class
    phClass[0] = PhotonClass(preselTree->phEta[ph1], preselTree->phR9[ph1]);

    //2011 selection
    if(name == "data2011_10GeV" || name == "data2011_15GeV")
    {
      if(!SelectPh11(phClass[0], preselTree, ph1))
        continue;
    }

    //2012 selection
    if(name == "data2012_10GeV" || name == "data2012_15GeV")
    {
      if(!SelectPh12(phClass[0], preselTree, ph1))
        continue;
    }
    //Set the TLorentzVector
    ph[0].SetPtEtaPhiM(preselTree->phPt[ph1], preselTree->phEta[ph1], preselTree->phPhi[ph1], 0.0);

    //loop over second photon
    for(int ph2 = ph1 + 1; ph2 < preselTree->Nph; ph2++)
    {
	  //photon position
      if(TMath::Abs(preselTree->phEta[ph2]) > 2.5)
        return 0;
      if(TMath::Abs(preselTree->phEta[ph2]) > 1.44 && TMath::Abs(preselTree->phEta[ph2]) < 1.57)
        return 0;
      //get photon class
      phClass[1] = PhotonClass(preselTree->phEta[ph2], preselTree->phR9[ph2]);
      //2011 selection
      if(name == "data2011_10GeV" || name == "data2011_15GeV")
      {
        if(!SelectPh11(phClass[1], preselTree, ph2))
          continue;
      }
      //2012 selection
      if(name == "data2012_10GeV" || name == "data2012_15GeV")
      {
        if(!SelectPh12(phClass[1], preselTree, ph2))
          continue;
      }
      //Set TLorentzVector
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
