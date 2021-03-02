// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file taskBplus.cxx
/// \brief B+ analysis task
///
/// \author Gian Michele Innocenti <gian.michele.innocenti@cern.ch>, CERN
/// \author Vít Kučera <vit.kucera@cern.ch>, CERN
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>, CERN
/// \author Deepa Thomas <deepa.thomas@cern.ch>, UT Austin

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisDataModel.h"
#include "AnalysisDataModel/HFSecondaryVertex.h"
#include "AnalysisDataModel/HFCandidateSelectionTables.h"
#include "AnalysisCore/trackUtilities.h"
#include "DetectorsVertexing/DCAFitterN.h"
#include "ReconstructionDataFormats/DCA.h"
#include "AnalysisDataModel/TrackSelectionTables.h"
#include "ReconstructionDataFormats/V0.h"


using namespace o2;
using namespace o2::framework;
using namespace o2::aod::hf_cand_prong2;
using namespace o2::framework::expressions;

void customize(std::vector<o2::framework::ConfigParamSpec>& workflowOptions)
{
  ConfigParamSpec optionDoMC{"doMC", VariantType::Bool, false, {"Fill MC histograms."}};
  workflowOptions.push_back(optionDoMC);
}

#include "Framework/runDataProcessing.h"

namespace o2::aod
{
  namespace extra
  {
    DECLARE_SOA_INDEX_COLUMN(Collision, collision);
  }
  DECLARE_SOA_TABLE(Colls, "AOD", "COLLSID", o2::aod::extra::CollisionId);
} // namespace o2::aod
struct AddCollisionId {
  Produces<o2::aod::Colls> colls;
  void process(aod::HfCandProng2 const& candidates, aod::Tracks const&)
  {
    for (auto& candidate : candidates) {
      colls(candidate.index0_as<aod::Tracks>().collisionId());
    }
  }
};

/// B+ analysis task
struct TaskBplus {
  HistogramRegistry registry{
    "registry",
      {{"hmassD0", "2-prong candidates; inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{500, 0., 5.}}}},
        {"hptD0prong0", "2-prong candidates; prong0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 10.}}}},
        {"hptD0prong1", "2-prong candidates; prong1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 10.}}}},
        {"hEtaD0", "2-prong candidates; candidate #it{#eta};entries", {HistType::kTH1F, {{100, -2, 2.}}}},
        {"hptD0","D0 Candidates; #it{p}_{T} (GeV/#it{c});entries",{HistType::kTH1F, {{100, 0., 10.}}}},
        {"hptBachTrk","pi track; #it{p}_{T} (GeV/#it{c});entries",{HistType::kTH1F, {{100, 0., 10.}}}},
        {"hptD0inTrkLoop","D0 Cand(inside track loop); #it{p}_{T} (GeV/#it{c});entries",{HistType::kTH1F, {{100, 0., 10.}}}},
        {"hptcand", "D0+pi candidates; #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 10.}}}},
        {"hpTBCand","B Candidates; #it{p}_{T} (GeV/#it{c});entries",{HistType::kTH1F, {{100, 0., 10.}}}},
        {"hpTBprongD0","B prong D0; #it{p}_{T} (GeV/#it{c});entries",{HistType::kTH1F, {{100, 0., 10.}}}},
        {"hpTBprongpi","B prong pi; #it{p}_{T} (GeV/#it{c});entries",{HistType::kTH1F, {{100, 0., 10.}}}},
        {"hBCandInvmass","B candidates; Invmass (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{1000, 0., 10.}}}}
      }
  };

  Configurable<int> d_selectionFlagD0{"d_selectionFlagD0", 1, "Selection Flag for D0"};
  Configurable<int> d_selectionFlagD0bar{"d_selectionFlagD0bar", 1, "Selection Flag for D0bar"};
  Configurable<double> cutEtaCandMax{"cutEtaCandMax", 1., "max. cand. pseudorapidity"};

  // vertexing parameters
  Configurable<double> d_bz{"d_bz", 5., "magnetic field"};
  Configurable<bool> b_propdca{"b_propdca", true, "create tracks version propagated to PCA"};
  Configurable<double> d_maxr{"d_maxr", 200., "reject PCA's above this radius"};
  Configurable<double> d_maxdzini{"d_maxdzini", 4., "reject (if>0) PCA candidate if tracks DZ exceeds threshold"};
  Configurable<double> d_minparamchange{"d_minparamchange", 1.e-3, "stop iterations if largest change of any X is smaller than this"};
  Configurable<double> d_minrelchi2change{"d_minrelchi2change", 0.9, "stop iterations if chi2/chi2old > this"};
  Configurable<bool> d_UseAbsDCA{"d_UseAbsDCA", true, "Use Abs DCAs"};

  Filter filterSelectCandidates = (aod::hf_selcandidate_d0::isSelD0 >= d_selectionFlagD0 || aod::hf_selcandidate_d0::isSelD0bar >= d_selectionFlagD0bar);
  // Partition<aod::BigTracks> positiveTracks = aod::track::signed1Pt >= 0.f;
  // Partition<aod::BigTracks> negativeTracks = aod::track::signed1Pt < 0.f;

  double massPi = RecoDecay::getMassPDG(kPiPlus);
  double massD0 = RecoDecay::getMassPDG(421);

  void process(aod::Collision const& collisions, aod::BigTracks const& tracks, soa::Filtered<soa::Join<aod::HfCandProng2, aod::HFSelD0Candidate, aod::Colls>> const& candidates)
  {
    //Initialise fitter for B vertex
    o2::vertexing::DCAFitterN<2> bfitter;
    bfitter.setBz(d_bz);
    bfitter.setPropagateToPCA(b_propdca);
    bfitter.setMaxR(d_maxr);
    bfitter.setMinParamChange(d_minparamchange);
    bfitter.setMinRelChi2Change(d_minrelchi2change);
    bfitter.setUseAbsDCA(d_UseAbsDCA);

    //Initial fitter to redo D-vertex to get extrapolated daughter tracks
    o2::vertexing::DCAFitterN<2> df;
    df.setBz(d_bz);
    df.setPropagateToPCA(b_propdca);
    df.setMaxR(d_maxr);
    df.setMinParamChange(d_minparamchange);
    df.setMinRelChi2Change(d_minrelchi2change);
    df.setUseAbsDCA(d_UseAbsDCA);

    //Loop over D0 candi
    for (auto& candidate : candidates) {
      if (!(candidate.hfflag() & 1 << D0ToPiK)) {
        continue;
      }
      if (std::abs(candidate.eta()) > cutEtaCandMax) {
        continue;
      }
      if ((candidate.isSelD0bar() < d_selectionFlagD0bar) && (candidate.isSelD0() < d_selectionFlagD0))
        continue;

      if (candidate.isSelD0bar() >= d_selectionFlagD0bar)
        registry.fill(HIST("hmassD0"), InvMassD0bar(candidate));
      if(candidate.isSelD0() >= d_selectionFlagD0)
        registry.fill(HIST("hmassD0"), InvMassD0(candidate));
      registry.fill(HIST("hptD0prong0"), candidate.ptProng0());
      registry.fill(HIST("hptD0prong1"), candidate.ptProng1());
      registry.fill(HIST("hEtaD0"), candidate.eta());
      registry.fill(HIST("hptD0"),candidate.pt());

      const std::array<float, 3> vertexD0 = {candidate.xSecondaryVertex(), candidate.ySecondaryVertex(), candidate.zSecondaryVertex()};
      const std::array<float, 3> momentumD0 = {candidate.px(), candidate.py(), candidate.pz()};

      auto prong0 = candidate.index0_as<aod::BigTracks>();
      auto prong1 = candidate.index1_as<aod::BigTracks>();
      auto prong0TrackParCov = getTrackParCov(prong0);
      auto prong1TrackParCov = getTrackParCov(prong1);

      LOGF(INFO, "All track: %d (prong0); %d (prong1)", candidate.index0().globalIndex(), candidate.index1().globalIndex());
      LOGF(INFO, "All track pT: %f (prong0); %f (prong1)", prong0.pt(), prong1.pt());

      // reconstruct D0 secondary vertex
      if (df.process(prong0TrackParCov, prong1TrackParCov) == 0) {
        continue;
      }

      //Propogate prong tracks to extrapolated track position.
      prong0TrackParCov.propagateTo(df.getTrack(0).getX(), d_bz);
      prong1TrackParCov.propagateTo(df.getTrack(1).getX(), d_bz);

      // build a D0 neutral track
      auto trackD0 = o2::dataformats::V0(vertexD0, momentumD0, prong0TrackParCov, prong1TrackParCov, {0, 0}, {0, 0});

      //loop over tracks for pi selection
      auto count = 0;
      for (auto &track : tracks) {
        if (count % 100 == 0) {
          LOGF(INFO, "Col: %d (cand); %d (track)", candidate.collisionId(), track.collisionId());
          count++;
        }

        if(candidate.isSelD0() >= d_selectionFlagD0 && track.signed1Pt() > 0)
          continue; //to select D0pi- pair
        if(candidate.isSelD0bar() >= d_selectionFlagD0bar && track.signed1Pt() < 0)
          continue; //to select D0(bar)pi+ pair

        if(candidate.index0Id() == track.globalIndex() || candidate.index1Id() == track.globalIndex())
          continue; //daughter track id and bachelor track id not the same

        registry.fill(HIST("hptBachTrk"),track.pt());
        registry.fill(HIST("hptD0inTrkLoop"),candidate.pt());
        registry.fill(HIST("hptcand"), candidate.pt() + track.pt());

        auto bachTrack = getTrackParCov(track);

        std::array<float, 3> pvecD0 = {0., 0., 0.};
        std::array<float, 3> pvecbach = {0., 0., 0.};
        std::array<float, 3> pvecBCand = {0., 0., 0.};

        //find the DCA between the D0 and the bachelor track, for B+
        int nCand = bfitter.process(trackD0, bachTrack); //Plot nCand
        if(nCand == 0) continue;

        bfitter.propagateTracksToVertex();        // propagate the bachelor and D0 to the B+ vertex
        bfitter.getTrack(0).getPxPyPzGlo(pvecD0); //momentum of D0 at the B+ vertex
        bfitter.getTrack(1).getPxPyPzGlo(pvecbach); //momentum of pi+ at the B+ vertex

        pvecBCand = array{pvecbach[0] + pvecD0[0],
          pvecbach[1] + pvecD0[1],
          pvecbach[2] + pvecD0[2]};

        //B-candidate properties
        auto arrMom = array{pvecbach, pvecD0};
        double massBP = RecoDecay::M(arrMom, array{massPi, massD0}); //B+ mass
        double pTBP = RecoDecay::Pt(pvecBCand[0],pvecBCand[1]); //pT B+
        double pTBprongD0 = RecoDecay::Pt(pvecD0[0], pvecD0[1]); //pT of B decay prong D0
        double pTBprongpi = RecoDecay::Pt(pvecbach[0], pvecbach[1]); //pT of B decay prong pi

        registry.fill(HIST("hBCandInvmass"), massBP);
        registry.fill(HIST("hpTBCand"), pTBP);
        registry.fill(HIST("hpTBprongD0"), pTBprongD0);
        registry.fill(HIST("hpTBprongpi"), pTBprongpi);

      } //track loop
    }//D0 cand loop
  } //process
}; //struct

WorkflowSpec defineDataProcessing(ConfigContext const&)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<AddCollisionId>("hf-task-add-collisionId"),
      adaptAnalysisTask<TaskBplus>("hf-task-bplus"),
  };
  return workflow;
}
