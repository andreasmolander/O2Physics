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
#include <cstdint>
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

    // BC
    int runNumber = -1;
    uint64_t globalBC = 0;
    uint64_t ctpTriggerMask = 0;
    uint64_t ctpInputMask = 0;

    // Timestamps
    uint64_t timestamp = 0;

    // EvSel
    bool sel8 = false;
    bool hasFT0 = false;
    bool hasFV0 = false;
    bool hasFDD = false;

    // Collisions
    int32_t bcId = -1;
    float posX = -200;
    float posY = -200;
    float posZ = -200;
    uint16_t flags = 0;
    uint16_t numContrib = -1;
    float collisionTime = -200;
    float collisionTimeRes = -200;

    // FT0
    int32_t ft0BCId = -1;
    std::vector<float> ft0AmplitudeA;
    std::vector<int16_t> ft0ChannelA;
    std::vector<float> ft0AmplitudeC;
    std::vector<int16_t> ft0ChannelC;
    float ft0TimeA = -200;
    float ft0TimeC = -200;
    uint8_t ft0TriggerMask = 0;
    float ft0PosZ = -200;
    float ft0CollTime = -200;
    float ft0SumAmpA = 0;
    float ft0SumAmpC = 0;

    // FT0Corrected
    float t0ACorrected = -200;
    float t0CCorrected = -200;
    float t0AC = -200;
    float t0resolution = -200;

    // FT0 derived quantities
    std::vector<float> ft0ChAmpl(o2::aod::fit::nChFT0, 0);
    float ft0TotAmplA = 0;
    float ft0TotAmplC = 0;

    // FV0A
    int32_t fv0BCId = -1;
    std::vector<float> fv0Amplitude;
    std::vector<int16_t> fv0Channel;
    float fv0Time = -200;
    uint8_t fv0TriggerMask = 0;

    // FV0 derived quantities
    std::vector<float> fv0ChAmpl(o2::aod::fit::nChFV0, 0);
    float fv0TotAmpl = 0;

    // FDD
    int32_t fddBCId = -1;
    int16_t fddChargeA[8] = {0};
    int16_t fddChargeC[8] = {0};
    float fddTimeA = -200;
    float fddTimeC = -200;
    uint8_t fddTriggerMask = 0;

    // FDD derived quantities
    std::vector<float> fddChAmpl(o2::aod::fit::nChFDD, 0);
    float fddTotAmplA = 0;
    float fddTotAmplC = 0;

    for (const auto& collision : collisions) {
      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      
      // BC
      runNumber = bc.runNumber();
      globalBC = bc.globalBC();
      ctpTriggerMask = bc.triggerMask();
      ctpInputMask = bc.inputMask();
      
      // Timestamps
      timestamp = bc.timestamp();

      // EvSel
      sel8 = collision.sel8();
      hasFT0 = collision.has_foundFT0();
      hasFV0 = collision.has_foundFV0();
      hasFDD = collision.has_foundFDD();

      // Collisions
      bcId = collision.bcId();
      posX = collision.posX();
      posY = collision.posY();
      posZ = collision.posZ();
      flags = collision.flags();
      numContrib = collision.numContrib();
      collisionTime = collision.collisionTime();
      collisionTimeRes = collision.collisionTimeRes();

      // FT0
      ft0BCId = -1;
      ft0AmplitudeA.clear();
      ft0ChannelA.clear();
      ft0AmplitudeC.clear();
      ft0ChannelC.clear();
      ft0TimeA = -200;
      ft0TimeC = -200;
      ft0TriggerMask = 0;
      ft0PosZ = -200;
      ft0CollTime = -200;
      ft0SumAmpA = 0;
      ft0SumAmpC = 0;

      // FT0Corrected
      t0ACorrected = collision.t0ACorrected();
      t0CCorrected = collision.t0CCorrected();
      t0AC = collision.t0AC();
      t0resolution = collision.t0resolution();

      // FT0 derived quantities
      ft0ChAmpl.assign(o2::aod::fit::nChFT0, 0);
      ft0TotAmplA = 0;
      ft0TotAmplC = 0;

      // FV0A
      fv0BCId = -1;
      fv0Amplitude.clear();
      fv0Channel.clear();
      fv0Time = -200;
      fv0TriggerMask = 0;

      // FV0 derived quantities
      fv0ChAmpl.assign(o2::aod::fit::nChFV0, 0);
      fv0TotAmpl = 0;

      // FDD
      fddBCId = -1;
      std::fill(std::begin(fddChargeA), std::end(fddChargeA), 0);
      std::fill(std::begin(fddChargeC), std::end(fddChargeC), 0);
      fddTimeA = -200;
      fddTimeC = -200;
      fddTriggerMask = 0;

      // FDD derived quantities
      fddChAmpl.assign(o2::aod::fit::nChFDD, 0);
      fddTotAmplA = 0;
      fddTotAmplC = 0;

      if (hasFT0) {
        // FT0
        auto ft0 = collision.foundFT0();
        ft0BCId = ft0.bcId();
        ft0TimeA = ft0.timeA();
        ft0TimeC = ft0.timeC();
        ft0TriggerMask = ft0.triggerMask();
        ft0PosZ = ft0.posZ();
        ft0CollTime = ft0.collTime();
        ft0SumAmpA = ft0.sumAmpA();
        ft0SumAmpC = ft0.sumAmpC();

        for (size_t i = 0; i < ft0.amplitudeA().size(); i++) {
          // FT0
          ft0AmplitudeA.push_back(ft0.amplitudeA()[i]);
          ft0ChannelA.push_back(ft0.channelA()[i]);

          // FT0 derived quantities
          ft0ChAmpl[ft0.channelA()[i]] = ft0.amplitudeA()[i];
          ft0TotAmplA += ft0.amplitudeA()[i];
        }
        for (size_t i = 0; i < ft0.amplitudeC().size(); i++) {
          // FT0
          ft0AmplitudeC.push_back(ft0.amplitudeC()[i]);
          ft0ChannelC.push_back(ft0.channelC()[i]);

          // FT0 derived quantities
          ft0ChAmpl[ft0.channelC()[i] + o2::aod::fit::nChFT0A] = ft0.amplitudeC()[i]; // Channel IDs in the C-side array start from zero in AO2D (JIRA AFIT-129)
          ft0TotAmplC += ft0.amplitudeC()[i];
        }
      }

      if (hasFV0) {
        // FV0A
        auto fv0 = collision.foundFV0();
        fv0BCId = fv0.bcId();
        fv0Time = fv0.time();
        fv0TriggerMask = fv0.triggerMask();
        
        for (size_t i = 0; i < fv0.amplitude().size(); i++) {
          // FV0A
          fv0Amplitude.push_back(fv0.amplitude()[i]);
          fv0Channel.push_back(fv0.channel()[i]);

          // FV0 derived quantities
          fv0ChAmpl[fv0.channel()[i]] = fv0.amplitude()[i];
          fv0TotAmpl += fv0.amplitude()[i];
        }
      }

      if (hasFDD) {
        // FDD
        auto fdd = collision.foundFDD();
        fddBCId = fdd.bcId();
        fddTimeA = fdd.timeA();
        fddTimeC = fdd.timeC();
        fddTriggerMask = fdd.triggerMask();

        for (size_t i = 0; i < 8; i++) {
          // FDD
          fddChargeA[i] = fdd.chargeA()[i];
          fddChargeC[i] = fdd.chargeC()[i];

          // FDD derived quantities
          fddChAmpl[i + 8] = fdd.chargeA()[i];
          fddChAmpl[i] = fdd.chargeC()[i];
          fddTotAmplA += fdd.chargeA()[i];
          fddTotAmplC += fdd.chargeC()[i];
        }
      }

      table(runNumber, globalBC, ctpTriggerMask, ctpInputMask,
            timestamp,
            sel8, hasFT0, hasFV0, hasFDD,
            bcId, posX, posY, posZ, flags, numContrib, collisionTime, collisionTimeRes,
            ft0BCId, ft0AmplitudeA, ft0ChannelA, ft0AmplitudeC, ft0ChannelC,
            ft0TimeA, ft0TimeC, ft0TriggerMask, ft0PosZ, ft0CollTime, ft0SumAmpA, ft0SumAmpC,
            t0ACorrected, t0CCorrected, t0AC, t0resolution,
            ft0ChAmpl, ft0TotAmplA, ft0TotAmplC,
            fv0BCId, fv0Amplitude, fv0Channel, fv0Time, fv0TriggerMask,
            fv0ChAmpl, fv0TotAmpl,
            fddBCId, fddChargeA, fddChargeC, fddTimeA, fddTimeC, fddTriggerMask,
            fddChAmpl, fddTotAmplA, fddTotAmplC);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<fitExtraTable>(cfgc)};
  return workflow;
}
