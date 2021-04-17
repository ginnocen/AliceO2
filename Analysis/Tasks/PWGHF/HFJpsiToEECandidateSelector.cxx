// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file HFJpsiToEECandidateSelector.cxx
/// \brief Jpsi selection task.
/// \author Biao Zhang <biao.zhang@cern.ch>, CCNU
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>, CERN

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "ALICE3Analysis/RICH.h"
#include "AnalysisCore/HFSelectorCuts.h"
#include "AnalysisDataModel/HFSecondaryVertex.h"
#include "AnalysisDataModel/HFCandidateSelectionTables.h"
using namespace o2;
using namespace o2::framework;
using namespace o2::aod::hf_cand_prong2;
using namespace o2::aod::alice3rich;
using namespace o2::analysis;
using namespace o2::analysis::hf_cuts_jpsi_toee;

/// Struct for applying Jpsi selection cuts

struct HFJpsiToEECandidateSelector {

  Produces<aod::HFSelJpsiToEECandidate> hfSelJpsiToEECandidate;

  Configurable<double> d_pTCandMin{"d_pTCandMin", 0., "Lower bound of candidate pT"};
  Configurable<double> d_pTCandMax{"d_pTCandMax", 50., "Upper bound of candidate pT"};

  Configurable<double> d_pidTPCMinpT{"d_pidTPCMinpT", 0.15, "Lower bound of track pT for TPC PID"};
  Configurable<double> d_pidTPCMaxpT{"d_pidTPCMaxpT", 10., "Upper bound of track pT for TPC PID"};
  Configurable<double> d_pidTOFMinpT{"d_pidTOFMinpT", 0.15, "Lower bound of track pT for TOF PID"};
  Configurable<double> d_pidTOFMaxpT{"d_pidTOFMaxpT", 5., "Upper bound of track pT for TOF PID"};

  Configurable<double> d_pidRHICMinpT{"d_pidRHICMinpT", 0.15, "Lower bound of track pT for RHIC PID"};
  Configurable<double> d_pidRHICMaxpT{"d_pidRHICMaxpT", 10., "Upper bound of track pT for RHICPID"};

  Configurable<double> d_TPCNClsFindablePIDCut{"d_TPCNClsFindablePIDCut", 70., "Lower bound of TPC findable clusters for good PID"};
  Configurable<double> d_nSigmaTPC{"d_nSigmaTPC", 3., "Nsigma cut on TPC only"};
  Configurable<double> d_nSigmaTOF{"d_nSigmaTOF", 3., "Nsigma cut on TOF only"};
  Configurable<double> d_nSigmaRHIC{"d_nSigmaRHIC", 3., "Nsigma cut on RHIC only"};
  Configurable<std::vector<double>> pTBins{"pTBins", std::vector<double>{hf_cuts_jpsi_toee::pTBins_v}, "pT bin limits"};
  Configurable<LabeledArray<double>> cuts{"Jpsi_to_ee_cuts", {hf_cuts_jpsi_toee::cuts[0], npTBins, nCutVars, pTBinLabels, cutVarLabels}, "Jpsi candidate selection per pT bin"};

  /// Selection on goodness of daughter tracks
  /// \note should be applied at candidate selection
  /// \param track is daughter track
  /// \return true if track is good
  template <typename T>
  bool daughterSelection(const T& track)
  {
    /*if (track.tpcNClsFound() == 0) {
      return false; //is it clusters findable or found - need to check
      }*/
    return true;
  }

  /// Conjugate independent toplogical cuts
  /// \param hfCandProng2 is candidate
  /// \param trackPositron is the track with the positron hypothesis
  /// \param trackElectron is the track with the electron hypothesis
  /// \return true if candidate passes all cuts
  template <typename T1, typename T2>
  bool selectionTopol(const T1& hfCandProng2, const T2& trackPositron, const T2& trackElectron)
  {
    auto candpT = hfCandProng2.pt();
    auto pTBin = findBin(pTBins, candpT);
    if (pTBin == -1) {
      return false;
    }

    if (candpT < d_pTCandMin || candpT >= d_pTCandMax) {
      return false; //check that the candidate pT is within the analysis range
    }

    if (TMath::Abs(InvMassJpsiToEE(hfCandProng2) - RecoDecay::getMassPDG(pdg::code::kJpsi)) > cuts->get(pTBin, "m")) {
      return false;
    }

    if (trackElectron.pt() < cuts->get(pTBin, "pT El") || trackPositron.pt() < cuts->get(pTBin, "pT El")) {
      return false; //cut on daughter pT
    }
    if (TMath::Abs(trackElectron.dcaPrim0()) > cuts->get(pTBin, "DCA_xy") || TMath::Abs(trackPositron.dcaPrim0()) > cuts->get(pTBin, "DCA_xy")) {
      return false; //cut on daughter dca - need to add secondary vertex constraint here
    }
    if (TMath::Abs(trackElectron.dcaPrim1()) > cuts->get(pTBin, "DCA_z") || TMath::Abs(trackPositron.dcaPrim1()) > cuts->get(pTBin, "DCA_z")) {
      return false; //cut on daughter dca - need to add secondary vertex constraint here
    }

    return true;
  }

  /// Check if track is ok for TPC PID
  /// \param track is the track
  /// \note function to be expanded
  /// \return true if track is ok for TPC PID
  template <typename T>
  bool validTPCPID(const T& track)
  {
    if (TMath::Abs(track.pt()) < d_pidTPCMinpT || TMath::Abs(track.pt()) >= d_pidTPCMaxpT) {
      return false;
    }
    //if (track.TPCNClsFindable() < d_TPCNClsFindablePIDCut) return false;
    return true;
  }

  /// Check if track is ok for TOF PID
  /// \param track is the track
  /// \note function to be expanded
  /// \return true if track is ok for TOF PID
  template <typename T>
  bool validTOFPID(const T& track)
  {
    if (TMath::Abs(track.pt()) < d_pidTOFMinpT || TMath::Abs(track.pt()) >= d_pidTOFMaxpT) {
      return false;
    }
    return true;
  }

  /// Check if track is ok for TOF PID
  /// \param track is the track
  /// \note function to be expanded
  /// \return true if track is ok for TOF PID
  template <typename T>
  bool validRHICPID(const T& track)
  {
    if (TMath::Abs(track.pt()) < d_pidRHICMinpT || TMath::Abs(track.pt()) >= d_pidRHICMaxpT) {
      return false;
    }
    return true;
  }

  /// Check if track is compatible with given TPC Nsigma cut for the electron hypothesis
  /// \param track is the track
  /// \param nSigmaCut is the nsigma threshold to test against
  /// \return true if track satisfies TPC PID hypothesis for given Nsigma cut
  template <typename T>
  bool selectionPIDTPC(const T& track, int nSigmaCut)
  {
    if (nSigmaCut > 999.) {
      return true;
    }
    return track.tpcNSigmaEl() < nSigmaCut;
  }

  /// Check if track is compatible with given TOF NSigma cut for a given flavour hypothesis
  /// \param track is the track
  /// \param nPDG is the flavour hypothesis PDG number
  /// \param nSigmaCut is the nSigma threshold to test against
  /// \note nPDG=211 pion  nPDG=321 kaon
  /// \return true if track satisfies TOF PID hypothesis for given NSigma cut
  template <typename T>
  bool selectionPIDTOF(const T& track, double nSigmaCut)
  {
    if (nSigmaCut > 999.) {
      return true;
    }
    return track.tofNSigmaEl() < nSigmaCut;
  }

  /// Check if track is compatible with given TOF NSigma cut for a given flavour hypothesis
  /// \param track is the track
  /// \param nPDG is the flavour hypothesis PDG number
  /// \param nSigmaCut is the nSigma threshold to test against
  /// \note nPDG=211 pion  nPDG=321 kaon
  /// \return true if track satisfies TOF PID hypothesis for given NSigma cut
  template <typename T>
  bool selectionPIDRICH(const T& track, double nSigmaCut)
  {
    if (nSigmaCut > 999.) {
      return true;
    }
    Printf("track.richNsigmaEl(): %f track.globalIndex(): %lld",track.richNsigmaEl(),track.globalIndex());
    return track.richNsigmaEl() < nSigmaCut;
  }
  /// PID selection on daughter track
  /// \param track is the daughter track
  /// \return 1 if successful PID match, 0 if successful PID rejection, -1 if no PID info
  template <typename T>
  int selectionPID(const T& track)
  {
    int statusTPC = -1;
    int statusTOF = -1;

    /*   if (validTPCPID(track)) {
      if (!selectionPIDTPC(track, d_nSigmaTPC)) {

        statusTPC = 0; //rejected by PID
      } else {
        statusTPC = 1; //positive PID
      }
    } else {
      statusTPC = -1; //no PID info
    }
*/
    if (validTOFPID(track)) {
      if (!selectionPIDTOF(track, d_nSigmaTOF)) {

        statusTOF = 0; //rejected by PID
      } else {
        statusTOF = 1; //positive PID
      }
    } else {
      statusTOF = -1; //no PID info
    }

    /* if (statusTPC == 1 && statusTOF == 1) {
      return 1; //what if we have 2 && 0 ?
    } else if (statusTPC == 0 || statusTOF == 0) {
      return 0;
    } else {
      return -1;
    }
    */
    if (statusTOF == 1) {
      return 1; //what if we have 2 && 0 ?
    } else {
      return 0;
    }
  }

    using Trks = soa::Join<aod::Tracks,aod::BigTracksPID, aod::TracksExtra, aod::RICHs>;

  //void process(aod::HfCandProng2 const& hfCandProng2s, aod::BigTracksPID const& tracks, aod::RICHs const& tracks_rich)
  void process(aod::HfCandProng2 const& hfCandProng2s, Trks const& tracks)
    {

    for (auto& hfCandProng2 : hfCandProng2s) { //looping over 2 prong candidates

      auto trackPos = hfCandProng2.index0_as<Trks>(); //positive daughter
      auto trackNeg = hfCandProng2.index1_as<Trks>(); //negative daughter

      if (!(hfCandProng2.hfflag() & 1 << JpsiToEE)) {
        hfSelJpsiToEECandidate(0);
        continue;
      }

      // daughter track validity selection
      if (!daughterSelection(trackPos) || !daughterSelection(trackNeg)) {
        hfSelJpsiToEECandidate(0);
        continue;
      }

      //implement filter bit 4 cut - should be done before this task at the track selection level
      //need to add special cuts (additional cuts on decay length and d0 norm)

      if (!selectionTopol(hfCandProng2, trackPos, trackNeg)) {
        hfSelJpsiToEECandidate(0);
        continue;
      }
      
      selectionPIDRICH(trackPos, d_nSigmaRHIC);
      //selectionPIDRICH(trackNeg, d_nSigmaRHIC);
      
      bool pidrichPos = false;
      bool pidrichNeg = false;
      bool pidrich = false;
     /* for (const auto& rich : tracks_rich) {

        const auto track = rich.track();
        if (trackPos.globalIndex() == track.globalIndex() && validRHICPID(trackPos) && selectionPIDRICH(rich, d_nSigmaRHIC)) {
          pidrichPos = true;
        }
        if (trackNeg.globalIndex() == track.globalIndex() && validRHICPID(trackNeg) && selectionPIDRICH(rich, d_nSigmaRHIC)) {
          pidrichNeg = true;
	}
        if (pidrichNeg == false || pidrichPos == false) {
          pidrich = false;
        }
      }

      if ((selectionPID(trackPos) == 0 || selectionPID(trackNeg) == 0) || pidrich == false) {
        hfSelJpsiToEECandidate(0);
        continue;
      }
*/
 /* if (selectionPID(trackPos) == 0 || selectionPID(trackNeg) == 0)
	{
	hfSelJpsiToEECandidate(0);
	continue;
	}
*/	
      hfSelJpsiToEECandidate(1);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HFJpsiToEECandidateSelector>(cfgc, TaskName{"hf-jpsi-toee-candidate-selector"})};
}
