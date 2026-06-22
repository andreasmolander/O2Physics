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

/// \file   fitAll.cxx
/// \brief  FITAll table definition and producer. Standalone derived data used for FIT studies.
///
/// \author Andreas Molander andreas.molander@cern.ch

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/FT0Corrected.h"
#include "Common/DataModel/Multiplicity.h"

#include <FDDBase/Constants.h>
#include <FV0Base/Constants.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

using namespace o2;
using namespace o2::framework;

namespace o2::aod
{
namespace fit
{

// Constants
static constexpr int NchFt0 = 208;                              ///< Number of FT0 channels
static constexpr int NchFt0A = 96;                              ///< Number of FT0A channels
static constexpr int NchFv0 = o2::fv0::Constants::nFv0Channels; ///< Number of FV0 channels
static constexpr int NchFdd = o2::fdd::Nchannels;               ///< Number of FDD channels

// Collisions
DECLARE_SOA_COLUMN(PosX, posX, float);
DECLARE_SOA_COLUMN(PosY, posY, float);
DECLARE_SOA_COLUMN(PosZ, posZ, float);
DECLARE_SOA_COLUMN(Flags, flags, uint16_t);
DECLARE_SOA_COLUMN(NumContrib, numContrib, uint16_t);
DECLARE_SOA_COLUMN(CollisionTime, collisionTime, float);
DECLARE_SOA_COLUMN(CollisionTimeRes, collisionTimeRes, float);

// Collisions BCs with timestamps
// TODO: remove these, for temporary debugging. EvSels BCs with timestamps
// are the ones to be used
DECLARE_SOA_COLUMN(CollRunNumber, collRunNumber, int);
DECLARE_SOA_COLUMN(CollGlobalBc, collGlobalBc, uint64_t);
DECLARE_SOA_COLUMN(CollCtpTriggerMask, collCtpTriggerMask, uint64_t);
DECLARE_SOA_COLUMN(CollCtpInputMask, collCtpInputMask, uint64_t);
DECLARE_SOA_COLUMN(CollTimestamp, collTimestamp, uint64_t);

// EvSels
DECLARE_SOA_BITMAP_COLUMN(TriggerAlias, triggerAlias, 32);
DECLARE_SOA_BITMAP_COLUMN(SelectionFlags, selectionFlags, 64);
DECLARE_SOA_BITMAP_COLUMN(RctFlags, rctFlags, 32);
DECLARE_SOA_COLUMN(Sel8, sel8, bool);
DECLARE_SOA_COLUMN(HasFoundBc, hasFoundBc, bool);  // TODO: remove? At the moment always true
DECLARE_SOA_COLUMN(HasFoundFt0, hasFoundFt0, bool);
DECLARE_SOA_COLUMN(HasFoundFv0, hasFoundFv0, bool);
DECLARE_SOA_COLUMN(HasFoundFdd, hasFoundFdd, bool);
DECLARE_SOA_COLUMN(HasFoundZdc, hasFoundZdc, bool);

// EvSels BCs with timestamps
DECLARE_SOA_COLUMN(EvSelRunNumber, evSelRunNumber, int);
DECLARE_SOA_COLUMN(EvSelGlobalBc, evSelGlobalBc, uint64_t);
DECLARE_SOA_COLUMN(EvSelCtpTriggerMask, evSelCtpTriggerMask, uint64_t);
DECLARE_SOA_COLUMN(EvSelCtpInputMask, evSelCtpInputMask, uint64_t);
DECLARE_SOA_COLUMN(EvSelTimestamp, evSelTimestamp, uint64_t);

// Mults
DECLARE_SOA_COLUMN(MultFt0A, multFt0A, float);
DECLARE_SOA_COLUMN(MultFt0C, multFt0C, float);
DECLARE_SOA_COLUMN(MultFv0A, multFv0A, float);
DECLARE_SOA_COLUMN(MultFddA, multFddA, float);
DECLARE_SOA_COLUMN(MultFddC, multFddC, float);
DECLARE_SOA_COLUMN(MultZnA, multZnA, float);
DECLARE_SOA_COLUMN(MultZnC, multZnC, float);
DECLARE_SOA_COLUMN(MultZem1, multZem1, float);
DECLARE_SOA_COLUMN(MultZem2, multZem2, float);
DECLARE_SOA_COLUMN(MultZpA, multZpA, float);
DECLARE_SOA_COLUMN(MultZpC, multZpC, float);

// FT0s
DECLARE_SOA_COLUMN(Ft0GlobalBc, ft0GlobalBc, uint64_t);
DECLARE_SOA_COLUMN(Ft0AmplitudeA, ft0AmplitudeA, std::vector<float>);
DECLARE_SOA_COLUMN(Ft0ChannelA, ft0ChannelA, std::vector<uint8_t>);
DECLARE_SOA_COLUMN(Ft0AmplitudeC, ft0AmplitudeC, std::vector<float>);
DECLARE_SOA_COLUMN(Ft0ChannelC, ft0ChannelC, std::vector<uint8_t>);
DECLARE_SOA_COLUMN(Ft0TimeA, ft0TimeA, float);
DECLARE_SOA_COLUMN(Ft0TimeC, ft0TimeC, float);
DECLARE_SOA_COLUMN(Ft0TriggerMask, ft0TriggerMask, uint8_t);
DECLARE_SOA_COLUMN(Ft0PosZ, ft0PosZ, float);
DECLARE_SOA_COLUMN(Ft0CollTime, ft0CollTime, float);
DECLARE_SOA_COLUMN(Ft0SumAmpA, ft0SumAmpA, float);
DECLARE_SOA_COLUMN(Ft0SumAmpC, ft0SumAmpC, float);

// FT0sCorrected
DECLARE_SOA_COLUMN(T0aCorrected, t0aCorrected, float);
DECLARE_SOA_COLUMN(T0cCorrected, t0cCorrected, float);
DECLARE_SOA_COLUMN(T0ac, t0ac, float);
DECLARE_SOA_COLUMN(T0resolution, t0resolution, float);

// FT0 derived quantities
DECLARE_SOA_COLUMN(Ft0ChAmpl, ft0ChAmpl, std::vector<float>); ///< FT0 channel amplitudes, vector idx = ch ID
DECLARE_SOA_COLUMN(Ft0TotAmplA, ft0TotAmplA, float);          ///< FT0-A total amplitude computed from channel amplitudes (for cross check)
DECLARE_SOA_COLUMN(Ft0TotAmplC, ft0TotAmplC, float);          ///< FT0-C total amplitude computed from channel amplitudes (for cross check)

// FV0As
DECLARE_SOA_COLUMN(Fv0GlobalBc, fv0GlobalBc, uint64_t);
DECLARE_SOA_COLUMN(Fv0Amplitude, fv0Amplitude, std::vector<float>);
DECLARE_SOA_COLUMN(Fv0Channel, fv0Channel, std::vector<uint8_t>);
DECLARE_SOA_COLUMN(Fv0Time, fv0Time, float);
DECLARE_SOA_COLUMN(Fv0TriggerMask, fv0TriggerMask, uint8_t);

// FV0 derived quantities
DECLARE_SOA_COLUMN(Fv0ChAmpl, fv0ChAmpl, std::vector<float>); ///< FV0 channel amplitudes, vector idx = ch ID
DECLARE_SOA_COLUMN(Fv0TotAmpl, fv0TotAmpl, float);            ///< FV0 total amplitude computed from channel amplitudes (for cross check)

// FDDs
DECLARE_SOA_COLUMN(FddGlobalBc, fddGlobalBc, uint64_t);
DECLARE_SOA_COLUMN(FddChargeA, fddChargeA, int16_t[8]);
DECLARE_SOA_COLUMN(FddChargeC, fddChargeC, int16_t[8]);
DECLARE_SOA_COLUMN(FddTimeA, fddTimeA, float);
DECLARE_SOA_COLUMN(FddTimeC, fddTimeC, float);
DECLARE_SOA_COLUMN(FddTriggerMask, fddTriggerMask, uint8_t);

// FDD derived quantities
DECLARE_SOA_COLUMN(FddChAmpl, fddChAmpl, std::vector<float>); ///< FDD channel amplitudes, vector idx = ch ID
DECLARE_SOA_COLUMN(FddTotAmplA, fddTotAmplA, float);          ///< FDD-A total amplitude computed from channel amplitudes (for cross check)
DECLARE_SOA_COLUMN(FddTotAmplC, fddTotAmplC, float);          ///< FDD-C total amplitude computed from channel amplitudes (for cross check)

// ZDCs
DECLARE_SOA_COLUMN(EnergyCommonZnA, energyCommonZnA, float);
DECLARE_SOA_COLUMN(EnergyCommonZnC, energyCommonZnC, float);
DECLARE_SOA_COLUMN(TimeZem1, timeZem1, float);
DECLARE_SOA_COLUMN(TimeZem2, timeZem2, float);

} // namespace fit

DECLARE_SOA_TABLE(FITsAll, "AOD", "FITALL", ///< Standalone derived data used for FIT studies.
                  fit::PosX, fit::PosY, fit::PosZ, fit::Flags, fit::NumContrib, fit::CollisionTime, fit::CollisionTimeRes,
                  fit::CollRunNumber, fit::CollGlobalBc, fit::CollCtpTriggerMask, fit::CollCtpInputMask, fit::CollTimestamp,
                  fit::TriggerAlias, fit::SelectionFlags, fit::RctFlags, fit::Sel8, fit::HasFoundBc, fit::HasFoundFt0, fit::HasFoundFv0, fit::HasFoundFdd, fit::HasFoundZdc,
                  fit::EvSelRunNumber, fit::EvSelGlobalBc, fit::EvSelCtpTriggerMask, fit::EvSelCtpInputMask, fit::EvSelTimestamp,
                  fit::MultFt0A, fit::MultFt0C, fit::MultFv0A, fit::MultFddA, fit::MultFddC, fit::MultZnA, fit::MultZnC, fit::MultZem1, fit::MultZem2, fit::MultZpA, fit::MultZpC,
                  fit::Ft0GlobalBc, fit::Ft0AmplitudeA, fit::Ft0ChannelA, fit::Ft0AmplitudeC, fit::Ft0ChannelC,
                  fit::Ft0TimeA, fit::Ft0TimeC, fit::Ft0TriggerMask, fit::Ft0PosZ, fit::Ft0CollTime, fit::Ft0SumAmpA, fit::Ft0SumAmpC,
                  fit::T0aCorrected, fit::T0cCorrected, fit::T0ac, fit::T0resolution,
                  fit::Ft0ChAmpl, fit::Ft0TotAmplA, fit::Ft0TotAmplC,
                  fit::Fv0GlobalBc, fit::Fv0Amplitude, fit::Fv0Channel, fit::Fv0Time, fit::Fv0TriggerMask,
                  fit::Fv0ChAmpl, fit::Fv0TotAmpl,
                  fit::FddGlobalBc, fit::FddChargeA, fit::FddChargeC, fit::FddTimeA, fit::FddTimeC, fit::FddTriggerMask,
                  fit::FddChAmpl, fit::FddTotAmplA, fit::FddTotAmplC,
                  fit::EnergyCommonZnA, fit::EnergyCommonZnC, fit::TimeZem1, fit::TimeZem2);

using FITAll = FITsAll::iterator;

} // namespace o2::aod

struct FitAll {
  
  enum statBins {
    kNcoll,
    kNhasFoundBc,
    kHasFt0,
    kHasFv0,
    kHasFdd,
    kCollBcIsEvSelBc,
    kCollBcNotEvSelBc,
    kCollBcIsFt0Bc,
    kCollBcNotFt0Bc,
    kEvSelBcIsFt0Bc,
    kEvSelBcNotFt0Bc,
    kEvSelBcIsFv0Bc,
    kEvSelBcNotFv0Bc,
    kEvSelBcIsFddBc,
    kEvSelBcNotFddBc,
    kNstatBins
  };

  // Configurables
  Configurable<bool> doQa{"doQA", false, "Add QA histograms"};

  // Producer
  Produces<o2::aod::FITsAll> table;

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::QAObject};

  void init(InitContext const&)
  {
    if (doQa) {
      // TODO: review, this is temporary debugging
      histos.add("stats", "stats", kTH1I, {{kNstatBins, 0, kNstatBins, ""}});
      auto h = histos.get<TH1>(HIST("stats"));
      h->GetXaxis()->SetBinLabel(kNcoll + 1, "Coll");
      h->GetXaxis()->SetBinLabel(kNhasFoundBc + 1, "HasFoundBC");
      h->GetXaxis()->SetBinLabel(kHasFt0 + 1, "HasFT0");
      h->GetXaxis()->SetBinLabel(kHasFv0 + 1, "HasFV0");
      h->GetXaxis()->SetBinLabel(kHasFdd + 1, "HasFDD");
      h->GetXaxis()->SetBinLabel(kCollBcIsEvSelBc + 1, "CollBC == EvSelBC");
      h->GetXaxis()->SetBinLabel(kCollBcNotEvSelBc + 1, "CollBC != EvSelBC");
      h->GetXaxis()->SetBinLabel(kCollBcIsFt0Bc + 1, "CollBC == FT0BC");
      h->GetXaxis()->SetBinLabel(kCollBcNotFt0Bc + 1, "CollBC != FT0BC");
      h->GetXaxis()->SetBinLabel(kEvSelBcIsFt0Bc + 1, "EvSelBC == FT0BC");
      h->GetXaxis()->SetBinLabel(kEvSelBcNotFt0Bc + 1, "EvSelBC != FT0BC");
      h->GetXaxis()->SetBinLabel(kEvSelBcIsFv0Bc + 1, "EvSelBC == FV0BC");
      h->GetXaxis()->SetBinLabel(kEvSelBcNotFv0Bc + 1, "EvSelBC != FV0BC");
      h->GetXaxis()->SetBinLabel(kEvSelBcIsFddBc + 1, "EvSelBC == FDDBC");
      h->GetXaxis()->SetBinLabel(kEvSelBcNotFddBc + 1, "EvSelBC != FDDBC");
    }
  }

  void process(soa::Join<aod::Collisions, aod::EvSels, aod::MultsRun3, aod::FT0sCorrected> const& collisions,
               aod::BCsWithTimestamps const&,
               aod::FT0s const&, aod::FV0As const&, aod::FDDs const&, aod::Zdcs const&)
  {
    table.reserve(collisions.size());

    for (const auto& collision : collisions) {
      if (doQa) {
        histos.fill(HIST("stats"), kNcoll);
      }

      // Collision
      float posX = collision.posX();
      float posY = collision.posY();
      float posZ = collision.posZ();
      uint16_t flags = collision.flags();
      uint16_t numContrib = collision.numContrib();
      float collisionTime = collision.collisionTime();
      float collisionTimeRes = collision.collisionTimeRes();

      // Collision BC with timestamp
      auto collBc = collision.bc_as<aod::BCsWithTimestamps>();
      int collRunNumber = collBc.runNumber();
      uint64_t collGlobalBc = collBc.globalBC();
      uint64_t collCtpTriggerMask = collBc.triggerMask();
      uint64_t collCtpInputMask = collBc.inputMask();
      uint64_t collTimestamp = collBc.timestamp();

      // EvSel
      uint32_t triggerAlias = collision.alias_raw();
      uint64_t selectionFlags = collision.selection_raw();
      uint32_t rctFlags = collision.rct_raw();
      bool sel8 = collision.sel8();
      bool hasFoundBc = collision.has_foundBC();
      bool hasFoundFt0 = collision.has_foundFT0();
      bool hasFoundFv0 = collision.has_foundFV0();
      bool hasFoundFdd = collision.has_foundFDD();
      bool hasFoundZdc = collision.has_foundZDC();

      // EvSel BC with timestamp
      int evSelRunNumber = -1;
      uint64_t evSelGlobalBc = 0;
      uint64_t evSelCtpTriggerMask = 0;
      uint64_t evSelCtpInputMask = 0;
      uint64_t evSelTimestamp = 0;
      if (hasFoundBc) {
        auto evSelBc = collision.foundBC_as<aod::BCsWithTimestamps>();
        evSelRunNumber = evSelBc.runNumber();
        evSelGlobalBc = evSelBc.globalBC();
        evSelCtpTriggerMask = evSelBc.triggerMask();
        evSelCtpInputMask = evSelBc.inputMask();
        evSelTimestamp = evSelBc.timestamp();

        if (doQa) {
          histos.fill(HIST("stats"), kNhasFoundBc);
          if (collGlobalBc == evSelGlobalBc) {
            histos.fill(HIST("stats"), kCollBcIsEvSelBc);
          } else {
            histos.fill(HIST("stats"), kCollBcNotEvSelBc);
          }
        }
      }

      // Multiplicities
      float multFT0A = collision.multFT0A();
      float multFT0C = collision.multFT0C();
      float multFV0A = collision.multFV0A();
      float multFDDA = collision.multFDDA();
      float multFDDC = collision.multFDDC();
      float multZNA = collision.multZNA();
      float multZNC = collision.multZNC();
      float multZEM1 = collision.multZEM1();
      float multZEM2 = collision.multZEM2();
      float multZPA = collision.multZPA();
      float multZPC = collision.multZPC();

      // FT0Corrected
      float t0ACorrected = collision.t0ACorrected();
      float t0CCorrected = collision.t0CCorrected();
      float t0AC = collision.t0AC();
      float t0resolution = collision.t0resolution();

      // FT0
      uint64_t ft0GlobalBc = 0;
      std::vector<float> ft0AmplitudeA;
      std::vector<uint8_t> ft0ChannelA;
      std::vector<float> ft0AmplitudeC;
      std::vector<uint8_t> ft0ChannelC;
      float ft0TimeA = -200;
      float ft0TimeC = -200;
      uint8_t ft0TriggerMask = 0;
      float ft0PosZ = -200;
      float ft0CollTime = -200;
      float ft0SumAmpA = 0;
      float ft0SumAmpC = 0;

      // FT0 derived quantities
      std::vector<float> ft0ChAmpl(o2::aod::fit::NchFt0, 0);
      float ft0TotAmplA = 0;
      float ft0TotAmplC = 0;

      if (hasFoundFt0) {
        // FT0
        auto ft0 = collision.foundFT0();
        auto ft0bc = ft0.bc_as<aod::BCsWithTimestamps>();
        ft0GlobalBc = ft0bc.globalBC();
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
          ft0ChAmpl[ft0.channelC()[i] + o2::aod::fit::NchFt0A] = ft0.amplitudeC()[i]; // Channel IDs in the C-side array start from zero in AO2D (JIRA AFIT-129)
          ft0TotAmplC += ft0.amplitudeC()[i];
        }

        if (doQa) {
          histos.fill(HIST("stats"), kHasFt0);
          if (collGlobalBc == ft0GlobalBc) {
            histos.fill(HIST("stats"), kCollBcIsFt0Bc);
          } else {
            histos.fill(HIST("stats"), kCollBcNotFt0Bc);
          }
          if (hasFoundBc) {
            if (evSelGlobalBc == ft0GlobalBc) {
              histos.fill(HIST("stats"), kEvSelBcIsFt0Bc);
            } else {
              histos.fill(HIST("stats"), kEvSelBcNotFt0Bc);
            }
          }
        }
      }

      // FV0A
      uint64_t fv0GlobalBc = 0;
      std::vector<float> fv0Amplitude;
      std::vector<uint8_t> fv0Channel;
      float fv0Time = -200;
      uint8_t fv0TriggerMask = 0;

      // FV0 derived quantities
      std::vector<float> fv0ChAmpl(o2::aod::fit::NchFv0, 0);
      float fv0TotAmpl = 0;

      if (hasFoundFv0) {
        // FV0A
        auto fv0 = collision.foundFV0();
        auto fv0bc = fv0.bc_as<aod::BCsWithTimestamps>();
        fv0GlobalBc = fv0bc.globalBC();
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

        if (doQa) {
          histos.fill(HIST("stats"), kHasFv0);
          if (hasFoundBc) {
            if (evSelGlobalBc == fv0GlobalBc) {
              histos.fill(HIST("stats"), kEvSelBcIsFv0Bc);
            } else {
              histos.fill(HIST("stats"), kEvSelBcNotFv0Bc);
            }
          }
        }
      }

      // FDD
      uint64_t fddGlobalBc = 0;
      int16_t fddChargeA[8] = {0};
      int16_t fddChargeC[8] = {0};
      float fddTimeA = -200;
      float fddTimeC = -200;
      uint8_t fddTriggerMask = 0;

      // FDD derived quantities
      std::vector<float> fddChAmpl(o2::aod::fit::NchFdd, 0);
      float fddTotAmplA = 0;
      float fddTotAmplC = 0;

      if (hasFoundFdd) {
        // FDD
        auto fdd = collision.foundFDD();
        auto fddbc = fdd.bc_as<aod::BCsWithTimestamps>();
        fddGlobalBc = fddbc.globalBC();
        fddTimeA = fdd.timeA();
        fddTimeC = fdd.timeC();
        fddTriggerMask = fdd.triggerMask();

        for (size_t i = 0; i < o2::aod::fit::NchFdd / 2; i++) {
          // FDD
          fddChargeA[i] = fdd.chargeA()[i];
          fddChargeC[i] = fdd.chargeC()[i];

          // FDD derived quantities
          fddChAmpl[i + o2::aod::fit::NchFdd / 2] = fdd.chargeA()[i];
          fddChAmpl[i] = fdd.chargeC()[i];
          fddTotAmplA += fdd.chargeA()[i];
          fddTotAmplC += fdd.chargeC()[i];
        }

        if (doQa) {
          histos.fill(HIST("stats"), kHasFdd);
          if (hasFoundBc) {
            if (evSelGlobalBc == fddGlobalBc) {
              histos.fill(HIST("stats"), kEvSelBcIsFddBc);
            } else {
              histos.fill(HIST("stats"), kEvSelBcNotFddBc);
            }
          }
        }
      }

      // ZDC
      float energyCommonZNA = -std::numeric_limits<float>::infinity();
      float energyCommonZNC = -std::numeric_limits<float>::infinity();
      float timeZEM1 = -std::numeric_limits<float>::infinity();
      float timeZEM2 = -std::numeric_limits<float>::infinity();

      if (hasFoundZdc) {
        // ZDC
        auto zdc = collision.foundZDC();
        energyCommonZNA = zdc.energyCommonZNA();
        energyCommonZNC = zdc.energyCommonZNC();
        timeZEM1 = zdc.timeZEM1();
        timeZEM2 = zdc.timeZEM2();
      }

      table(posX, posY, posZ, flags, numContrib, collisionTime, collisionTimeRes,
            collRunNumber, collGlobalBc, collCtpTriggerMask, collCtpInputMask, collTimestamp,
            triggerAlias, selectionFlags, rctFlags, sel8, hasFoundBc, hasFoundFt0, hasFoundFv0, hasFoundFdd, hasFoundZdc,
            evSelRunNumber, evSelGlobalBc, evSelCtpTriggerMask, evSelCtpInputMask, evSelTimestamp,
            multFT0A, multFT0C, multFV0A, multFDDA, multFDDC, multZNA, multZNC, multZEM1, multZEM2, multZPA, multZPC,
            ft0GlobalBc, ft0AmplitudeA, ft0ChannelA, ft0AmplitudeC, ft0ChannelC,
            ft0TimeA, ft0TimeC, ft0TriggerMask, ft0PosZ, ft0CollTime, ft0SumAmpA, ft0SumAmpC,
            t0ACorrected, t0CCorrected, t0AC, t0resolution,
            ft0ChAmpl, ft0TotAmplA, ft0TotAmplC,
            fv0GlobalBc, fv0Amplitude, fv0Channel, fv0Time, fv0TriggerMask,
            fv0ChAmpl, fv0TotAmpl,
            fddGlobalBc, fddChargeA, fddChargeC, fddTimeA, fddTimeC, fddTriggerMask,
            fddChAmpl, fddTotAmplA, fddTotAmplC,
            energyCommonZNA, energyCommonZNC, timeZEM1, timeZEM2);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<FitAll>(cfgc)};
  return workflow;
}
