// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

// HF Configurable Classes
//
// Authors: Nima Zardoshti

#ifndef O2_ANALYSIS_HFCONFIGURABLES_H
#define O2_ANALYSIS_HFCONFIGURABLES_H

#include <TMath.h>

class HFTrackIndexSkimsCreatorConfigs
{
 public:
  HFTrackIndexSkimsCreatorConfigs() = default;
  ~HFTrackIndexSkimsCreatorConfigs() = default;

  // 2-prong cuts D0
  double mPtD0Min = 0.;
  double mInvMassD0Min = 1.46;
  double mInvMassD0Max = 2.26;
  double mCPAD0Min = 0.75;
  double mImpParProductD0Max = -0.00005;
  // 2-prong cuts Jpsi
  double mPtJpsiMin = 0.;
  double mInvMassJpsiMin = 2.75;
  double mInvMassJpsiMax = 3.45;
  double mCPAJpsiMin = -2;
  double mImpParProductJpsiMax = 1000.;
  // 3-prong cuts - D+
  double mPtDPlusMin = 0.;
  double mInvMassDPlusMin = 1.7;
  double mInvMassDPlusMax = 2.05;
  double mCPADPlusMin = 0.5;
  double mDecLenDPlusMin = 0.;
  // 3-prong cuts - Lc
  double mPtLcMin = 0.;
  double mInvMassLcMin = 2.1;
  double mInvMassLcMax = 2.5;
  double mCPALcMin = 0.5;
  double mDecLenLcMin = 0.;
  // 3-prong cuts - Ds
  double mPtDsMin = 0.;
  double mInvMassDsMin = 1.7;
  double mInvMassDsMax = 2.2;
  double mCPADsMin = 0.5;
  double mDecLenDsMin = 0.;

  double getPtD0Min() const
  {
    return mPtD0Min;
  }

  double getInvMassD0Min() const
  {
    return mInvMassD0Min;
  }

  double getInvMassD0Max() const
  {
    return mInvMassD0Max;
  }

  double getCPAD0Min() const
  {
    return mCPAD0Min;
  }

  double getImpParProductD0Max() const
  {
    return mImpParProductD0Max;
  }

  double getPtJpsiMin() const
  {
    return mPtJpsiMin;
  }

  double getInvMassJpsiMin() const
  {
    return mInvMassJpsiMin;
  }

  double getInvMassJpsiMax() const
  {
    return mInvMassJpsiMax;
  }

  double getCPAJpsiMin() const
  {
    return mCPAJpsiMin;
  }

  double getImpParProductJpsiMax() const
  {
    return mImpParProductJpsiMax;
  }

  double getPtDPlusMin() const
  {
    return mPtDPlusMin;
  }

  double getInvMassDPlusMin() const
  {
    return mInvMassDPlusMin;
  }

  double getInvMassDPlusMax() const
  {
    return mInvMassDPlusMax;
  }

  double getCPADPlusMin() const
  {
    return mCPADPlusMin;
  }

  double getDecLenDPlusMin() const
  {
    return mDecLenDPlusMin;
  }

  double getPtLcMin() const
  {
    return mPtLcMin;
  }

  double getInvMassLcMin() const
  {
    return mInvMassLcMin;
  }

  double getInvMassLcMax() const
  {
    return mInvMassLcMax;
  }

  double getCPALcMin() const
  {
    return mCPALcMin;
  }

  double getDecLenLcMin() const
  {
    return mDecLenLcMin;
  }

  double getPtDsMin() const
  {
    return mPtDsMin;
  }

  double getInvMassDsMin() const
  {
    return mInvMassDsMin;
  }

  double getInvMassDsMax() const
  {
    return mInvMassDsMax;
  }

  double getCPADsMin() const
  {
    return mCPADsMin;
  }

  double getDecLenDsMin() const
  {
    return mDecLenDsMin;
  }

 private:
  ClassDef(HFTrackIndexSkimsCreatorConfigs, 1);
};

#endif
