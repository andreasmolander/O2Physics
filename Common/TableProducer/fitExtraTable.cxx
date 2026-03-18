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
#include <string>

using namespace o2;
using namespace o2::framework;

#include "Framework/runDataProcessing.h"

using BCsWithTimestamps = soa::Join<aod::BCs, aod::Timestamps>;

struct fitExtraTable {
  // Producer
  Produces<o2::aod::FITExtra> table;

  // using CollisionEvSel = soa::Join<aod::Collisions, aod::EvSels>::iterator;
  static constexpr float invLightSpeedCm2NS = 1.f / o2::constants::physics::LightSpeedCm2NS;

  void init(InitContext const&)
  {
  }
  
  void process(soa::Join<aod::Collisions, aod::EvSels, aod::FT0sCorrected> const& collisions,
               aod::BCsWithTimestamps const&,
               aod::FT0s const&, aod::FV0As const&, aod::FDDs const&) {
    table.reserve(collisions.size());

    int runnumber = -1;
    uint64_t globalBC = 0;
    uint64_t ctpTriggerMask = 0;
    uint64_t ctpInputMask = 0;

    int32_t bcId = -1;
    float posX = -200;
    float posY = -200;
    float posZ = -200;
    uint16_t flags = 0;
    int nContrib = -1;
    float collisionTime = -200;
    float collisionTimeRes = -200;

    float ft0timeA = -200;
    float ft0timeC = -200;
    float ft0timeACorr = -200;
    float ft0timeCCorr = -200;
    float ft0time = -200;
    float ft0timeRes = -200;
    float ft0vtx = -200;
    std::vector<float> ft0chAmpl(o2::aod::fit::nChFT0, 0);
    float ft0totAmplA = 0;
    float ft0totAmplC = 0;
    float ft0totAmplACheck = 0;
    float ft0totAmplCCheck = 0;

    float fv0time = -200;
    std::vector<float> fv0chAmpl(o2::aod::fit::nChFV0, 0);
    float fv0totAmpl = 0;

    float fddtimeA = -200;
    float fddtimeC = -200;
    std::vector<float> fddchAmpl(o2::aod::fit::nChFDD, 0);
    float fddtotAmplA = 0;
    float fddtotAmplC = 0;

    bool sel8, hasFT0, hasFV0, hasFDD;
    uint8_t ft0Triggers, fv0Triggers, fddTriggers;

    for (const auto& collision : collisions) {
      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      runnumber = bc.runNumber();
      globalBC = bc.globalBC();
      ctpTriggerMask = bc.triggerMask();
      ctpInputMask = bc.inputMask();

      bcId = collision.bcId();
      posX = collision.posX();
      posY = collision.posY();
      posZ = collision.posZ();
      flags = collision.flags();
      nContrib = collision.numContrib();
      collisionTime = collision.collisionTime();
      collisionTimeRes = collision.collisionTimeRes();

      sel8 = collision.sel8();

      ft0timeA = -200;
      ft0timeC = -200;
      ft0timeACorr = -200;
      ft0timeCCorr = -200;
      ft0time = -200;
      ft0timeRes = -200;
      ft0vtx = -200;
      std::fill(ft0chAmpl.begin(), ft0chAmpl.end(), 0);
      ft0totAmplA = 0;
      ft0totAmplC = 0;
      ft0totAmplACheck = 0;
      ft0totAmplCCheck = 0;

      fv0time = -200;
      std::fill(fv0chAmpl.begin(), fv0chAmpl.end(), 0);
      fv0totAmpl = 0;

      fddtimeA = -200;
      fddtimeC = -200;
      std::fill(fddchAmpl.begin(), fddchAmpl.end(), 0);
      fddtotAmplA = 0;
      fddtotAmplC = 0;

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
        for (size_t i = 0; i < ft0.amplitudeA().size(); i++) {
          ft0chAmpl[ft0.channelA()[i]] = ft0.amplitudeA()[i];
          ft0totAmplA += ft0.amplitudeA()[i];
        }
        for (size_t i = 0; i < ft0.amplitudeC().size(); i++) {
          ft0chAmpl[ft0.channelC()[i] + o2::aod::fit::nChFT0A] = ft0.amplitudeC()[i]; // Channel IDs in the C-side array start from zero in AO2D (JIRA AFIT-129)
          ft0totAmplC += ft0.amplitudeC()[i];
        }
        ft0totAmplA = ft0.sumAmpA();
        ft0totAmplC = ft0.sumAmpC();
      }

      if (hasFV0) {
        auto fv0 = collision.foundFV0();
        fv0time = fv0.time();
        fv0Triggers = fv0.triggerMask();
        for (size_t i = 0; i < fv0.amplitude().size(); i++) {
          fv0chAmpl[fv0.channel()[i]] = fv0.amplitude()[i];
          fv0totAmpl += fv0.amplitude()[i];
        }
      }

      if (hasFDD) {
        auto fdd = collision.foundFDD();
        fddtimeA = fdd.timeA();
        fddtimeC = fdd.timeC();
        fddTriggers = fdd.triggerMask();
        for (size_t i = 0; i < 8; i++) {
          fddchAmpl[i + 8] = fdd.chargeA()[i];
          fddtotAmplA += fdd.chargeA()[i];
        }
        for (size_t i = 0; i < 8; i++) {
          fddchAmpl[i] = fdd.chargeC()[i];
          fddtotAmplC += fdd.chargeC()[i]; 
        }
      }
  
      table(sel8, hasFT0, hasFV0, hasFDD,
            ft0Triggers, fv0Triggers, fddTriggers,
            runnumber, globalBC, ctpTriggerMask, ctpInputMask,
            bcId, posX, posY, posZ, flags, nContrib, collisionTime, collisionTimeRes,
            ft0timeA, ft0timeC, ft0timeACorr, ft0timeCCorr,
            ft0time, ft0timeRes, ft0vtx,
            ft0chAmpl, ft0totAmplA, ft0totAmplC, ft0totAmplACheck, ft0totAmplCCheck,
            fv0time, fv0chAmpl, fv0totAmpl,
            fddtimeA, fddtimeC, fddchAmpl, fddtotAmplA, fddtotAmplC);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<fitExtraTable>(cfgc)};
  return workflow;
}
