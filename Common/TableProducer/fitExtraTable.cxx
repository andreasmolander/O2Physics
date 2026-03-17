// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/FITExtra.h"
#include "Common/DataModel/FT0Corrected.h"

#include <CommonConstants/PhysicsConstants.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>

using namespace o2;
using namespace o2::framework;

#include "Framework/runDataProcessing.h"

struct fitExtraTable {
  // Producer
  Produces<o2::aod::FITExtra> table;

  // using CollisionEvSel = soa::Join<aod::Collisions, aod::EvSels>::iterator;
  static constexpr float invLightSpeedCm2NS = 1.f / o2::constants::physics::LightSpeedCm2NS;

  void init(InitContext const&)
  {
  }
  
  void process(soa::Join<aod::Collisions, aod::EvSels, aod::FT0sCorrected> const& collisions,
               aod::FT0s const&, aod::FV0As const&, aod::FDDs const&) {
    table.reserve(collisions.size());

    float pv = -200;
    int nContrib = -1;

    float ft0timeA = -200;
    float ft0timeC = -200;
    float ft0timeACorr = -200;
    float ft0timeCCorr = -200;
    float ft0time = -200;
    float ft0timeRes = -200;
    float ft0vtx = -200;
    float fv0time = -200;
    float fddtimeA = -200;
    float fddtimeC = -200;

    bool sel8, hasFT0, hasFV0, hasFDD;
    uint8_t ft0Triggers, fv0Triggers, fddTriggers;

    for (const auto& collision : collisions) {
      pv = collision.posZ();
      nContrib = collision.numContrib();
      sel8 = collision.sel8();

      ft0timeA = -200;
      ft0timeC = -200;
      ft0timeACorr = -200;
      ft0timeCCorr = -200;
      ft0time = -200;
      ft0timeRes = -200;
      ft0vtx = -200;
      fv0time = -200;
      fddtimeA = -200;
      fddtimeC = -200;

      hasFT0 = collision.has_foundFT0();
      hasFV0 = collision.has_foundFV0();
      hasFDD = collision.has_foundFDD();
      ft0Triggers = 0;
      fv0Triggers = 0;
      fddTriggers = 0;

      if (hasFT0) {
        auto ft0 = collision.foundFT0();
        ft0timeA = ft0.timeA();
        ft0timeC = ft0.timeC();
        ft0time = ft0.collTime();
        ft0vtx = ft0.posZ();
        ft0timeACorr = collision.t0ACorrected();
        ft0timeCCorr = collision.t0CCorrected();
        ft0timeRes = collision.t0resolution();
        ft0Triggers = ft0.triggerMask();
      }
      if (hasFV0) {
        auto fv0 = collision.foundFV0();
        fv0time = fv0.time();
        fv0Triggers = fv0.triggerMask();
      }
      if (hasFDD) {
        auto fdd = collision.foundFDD();
        fddtimeA = fdd.timeA();
        fddtimeC = fdd.timeC();
        fddTriggers = fdd.triggerMask();
      }
  
      // table(pv,
      //       nContrib,
      //       ft0timeA,
      //       ft0timeC,
      //       ft0timeACorr,
      //       ft0timeCCorr,
      //       ft0time,
      //       ft0timeRes,
      //       ft0vtx,
      //       fv0time,
      //       fddtimeA,
      //       fddtimeC,
      //       sel8,
      //       hasFT0,
      //       hasFV0,
      //       hasFDD,
      //       ft0Triggers,
      //       fv0Triggers,
      //       fddTriggers);
      table(sel8, hasFT0, hasFV0, hasFDD,
            ft0Triggers, fv0Triggers, fddTriggers,
            pv, nContrib,
            ft0timeA, ft0timeC, ft0timeACorr, ft0timeCCorr,
            ft0time, ft0timeRes, ft0vtx,
            fv0time, fddtimeA, fddtimeC);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<fitExtraTable>(cfgc)};
  return workflow;
}
