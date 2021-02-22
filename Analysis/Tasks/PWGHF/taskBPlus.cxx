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
      {{"hmassD0", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{500, 0., 5.}}}},
        {"hptcand", "B+ candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 10.}}}},
        {"hptD0prong0", "2-prong candidates;prong0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 10.}}}},
        {"hptD0prong1", "2-prong candidates;prong1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 10.}}}},
        {"hEtaD0", "2-prong candidates;candidate #it{#eta};entries", {HistType::kTH1F, {{100, -2, 2.}}}}
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

    //Define o2 fitter, 2-prong
    o2::vertexing::DCAFitterN<2> fitter;
    fitter.setBz(d_bz);
    fitter.setPropagateToPCA(b_propdca);
    fitter.setMaxR(d_maxr);
    fitter.setMinParamChange(d_minparamchange);
    fitter.setMinRelChi2Change(d_minrelchi2change);
    fitter.setUseAbsDCA(d_UseAbsDCA);

    //Loop over D0 candi
    for (auto& candidate : candidates) {
      if (!(candidate.hfflag() & 1 << D0ToPiK)) {
        continue;
      }
      if (std::abs(candidate.eta()) > cutEtaCandMax) {
        continue;
      }

      if (candidate.isSelD0bar() >= d_selectionFlagD0bar || candidate.isSelD0() >= d_selectionFlagD0)
      {
        if (candidate.isSelD0bar() >= d_selectionFlagD0bar)
          registry.fill(HIST("hmassD0"), InvMassD0bar(candidate));
        if(candidate.isSelD0() >= d_selectionFlagD0)
          registry.fill(HIST("hmassD0"), InvMassD0(candidate));
        registry.fill(HIST("hptD0prong0"), candidate.ptProng0());
        registry.fill(HIST("hptD0prong1"), candidate.ptProng1());
        registry.fill(HIST("hEtaD0"), candidate.eta());

        const std::array<float, 3> vertexD0 = {candidate.xSecondaryVertex(), candidate.ySecondaryVertex(), candidate.zSecondaryVertex()};
        const std::array<float, 3> momentumD0 = {candidate.px(), candidate.py(), candidate.pz()};

        auto prong0 = candidate.index0_as<aod::BigTracks>();
        auto prong1 = candidate.index1_as<aod::BigTracks>();

       // LOGF(INFO, "All track: %d (prong0); %d (prong1)", candidate.index0().globalIndex(), candidate.index1().globalIndex());
       // LOGF(INFO, "All track pT: %f (prong0); %f (prong1)", prong0.pt(), prong1.pt());

        auto prong0TrackParCov = getTrackParCov(prong0);
        prong0TrackParCov.propagateTo(prong0.x(), d_bz); //Lc has fitter.getTrack(0).getX(). Not the same as prong.x(). Need to include info to D0 table?
        auto prong1TrackParCov = getTrackParCov(prong1);
        prong1TrackParCov.propagateTo(prong1.x(), d_bz); //Lc has fitter.getTrack(0).getX(). Not the same as prong.x(). Need to include info to D0 table?

        //loop over tracks for pi selection
        auto count = 0;
        for (auto &track : tracks) {
          if (count % 100 == 0) {
            LOGF(INFO, "Col: %d (cand); %d (track)", candidate.collisionId(), track.collisionId());
            count++;
          }
          if ((candidate.isSelD0() >= d_selectionFlagD0 && track.signed1Pt() < 0) || candidate.isSelD0bar() >= d_selectionFlagD0bar && track.signed1Pt() > 0){ //D0pi- or D0(bar)pi+
            if(candidate.index0Id() == track.globalIndex() || candidate.index1Id() == track.globalIndex())
              continue; //daughter track id and bachelor track id not the same

            registry.fill(HIST("hptcand"), candidate.pt() + track.pt());

            auto bachTrack = getTrackParCov(track);

            // build the neutral track to then build the B
            auto trackD0 = o2::dataformats::V0(vertexD0, momentumD0, prong0TrackParCov, prong1TrackParCov, prong0.globalIndex(), prong1.globalIndex());

            std::array<float, 3> pvecD0 = {0., 0., 0.};
            std::array<float, 3> pvecbach = {0., 0., 0.};
            std::array<float, 3> pvecBCand = {0., 0., 0.};

            //find the DCA between the D0 and the bachelor track, for B+
            int nCand = fitter.process(trackD0, bachTrack); //Plot nCand

            if(nCand == 0) continue;

            fitter.propagateTracksToVertex();        // propagate the bachelor and D0 to the B+ vertex
            fitter.getTrack(0).getPxPyPzGlo(pvecD0); //momentum of D0 at the B+ vertex; Need to fill new D0+ pT
            fitter.getTrack(1).getPxPyPzGlo(pvecbach); //momentum of pi+ at the B+ vertex; Need to fill new p+ pT

            pvecBCand = array{pvecbach[0] + pvecD0[0],
              pvecbach[1] + pvecD0[1],
              pvecbach[2] + pvecD0[2]};

            // invariant mass
            auto arrMom = array{pvecbach, pvecD0};
            double massBP = RecoDecay::M(arrMom, array{massPi, massD0}); //Plot B+ mass

          } //D0 and pi charge
        } //track loop
      }//D0 selection
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
