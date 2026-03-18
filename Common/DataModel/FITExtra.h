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

#ifndef COMMON_DATAMODEL_FITEXTRA_H_
#define COMMON_DATAMODEL_FITEXTRA_H_

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include <cstdint>

namespace o2::aod
{
namespace fit
{

// Constants
// TODO: add reference channels?
static constexpr int nChFT0 = 208; ///< Number of FT0 channels
static constexpr int nChFT0A = 96; ///< Number of FT0A channels (-> number of FT0C channels = nChFT0 - nChFT0A)
static constexpr int nChFV0 = 48;  ///< Number of FV0 channels
static constexpr int nChFDD = 16;  ///< Number of FDD channels
static constexpr int nADC = 4096;  ///< Number of ADC channels

// Quantities copied straight from AOD
// TODO: do we need them here?

// BCs
DECLARE_SOA_COLUMN(RunNumber, runNumber, int);
DECLARE_SOA_COLUMN(GlobalBC, globalBC, uint64_t);
DECLARE_SOA_COLUMN(CTPtriggerMask, ctpTriggerMask, uint64_t);
DECLARE_SOA_COLUMN(CTPinputMask, ctpInputMask, uint64_t);

// Timestamps
DECLARE_SOA_COLUMN(Timestamp, timestamp, uint64_t);

// EvSels
DECLARE_SOA_COLUMN(Sel8, sel8, bool);
DECLARE_SOA_COLUMN(HasFT0, hasFT0, bool);
DECLARE_SOA_COLUMN(HasFV0, hasFV0, bool);
DECLARE_SOA_COLUMN(HasFDD, hasFDD, bool);

// Collisions
DECLARE_SOA_COLUMN(BCId, bcId, int32_t);
DECLARE_SOA_COLUMN(PosX, posX, float);
DECLARE_SOA_COLUMN(PosY, posY, float);
DECLARE_SOA_COLUMN(PosZ, posZ, float);
DECLARE_SOA_COLUMN(Flags, flags, uint16_t);
DECLARE_SOA_COLUMN(NumContrib, numContrib, uint16_t);
DECLARE_SOA_COLUMN(CollisionTime, collisionTime, float);
DECLARE_SOA_COLUMN(CollisionTimeRes, collisionTimeRes, float);

// FT0s
DECLARE_SOA_COLUMN(FT0BCId, ft0BCId, int32_t);
DECLARE_SOA_COLUMN(FT0AmplitudeA, ft0AmplitudeA, std::vector<float>);
DECLARE_SOA_COLUMN(FT0ChannelA, ft0ChannelA, std::vector<uint8_t>);
DECLARE_SOA_COLUMN(FT0AmplitudeC, ft0AmplitudeC, std::vector<float>);
DECLARE_SOA_COLUMN(FT0ChannelC, ft0ChannelC, std::vector<uint8_t>);
DECLARE_SOA_COLUMN(FT0TimeA, ft0TimeA, float);
DECLARE_SOA_COLUMN(FT0TimeC, ft0TimeC, float);
DECLARE_SOA_COLUMN(FT0TriggerMask, ft0TriggerMask, uint8_t);
DECLARE_SOA_COLUMN(FT0PosZ, ft0PosZ, float);
DECLARE_SOA_COLUMN(FT0CollTime, ft0CollTime, float);
DECLARE_SOA_COLUMN(FT0SumAmpA, ft0SumAmpA, float);
DECLARE_SOA_COLUMN(FT0SumAmpC, ft0SumAmpC, float);

// FT0sCorrected
DECLARE_SOA_COLUMN(T0ACorrected, t0ACorrected, float);
DECLARE_SOA_COLUMN(T0CCorrected, t0CCorrected, float);
DECLARE_SOA_COLUMN(T0AC, t0AC, float);
DECLARE_SOA_COLUMN(T0Resolution, t0resolution, float);

// FT0 derived quantities
DECLARE_SOA_COLUMN(FT0ChAmpl, ft0ChAmpl, std::vector<float>);  //! FT0 channel amplitudes, vector idx = ch ID
DECLARE_SOA_COLUMN(FT0TotAmplA, ft0TotAmplA, float); //! FT0-A total amplitude computed from channel amplitudes (for cross check)
DECLARE_SOA_COLUMN(FT0TotAmplC, ft0TotAmplC, float); //! FT0-C total amplitude computed from channel amplitudes(for cross check)

// FV0As
DECLARE_SOA_COLUMN(FV0BCId, fv0BCId, int32_t);
DECLARE_SOA_COLUMN(FV0Amplitude, fv0Amplitude, std::vector<float>);
DECLARE_SOA_COLUMN(FV0Channel, fv0Channel, std::vector<uint8_t>);
DECLARE_SOA_COLUMN(FV0Time, fv0Time, float);
DECLARE_SOA_COLUMN(FV0TriggerMask, fv0TriggerMask, uint8_t);

// FV0 derived quantities
DECLARE_SOA_COLUMN(FV0ChAmpl, fv0ChAmpl, std::vector<float>); //! FV0 channel amplitudes, vector idx = ch ID
DECLARE_SOA_COLUMN(FV0TotAmpl, fv0TotAmpl, float);

// FDD
DECLARE_SOA_COLUMN(FDDBCId, fddBCId, int32_t);
DECLARE_SOA_COLUMN(FDDChargeA, fddChargeA, int16_t[8]);
DECLARE_SOA_COLUMN(FDDChargeC, fddChargeC, int16_t[8]);
DECLARE_SOA_COLUMN(FDDTimeA, fddTimeA, float);
DECLARE_SOA_COLUMN(FDDTimeC, fddTimeC, float);
DECLARE_SOA_COLUMN(FDDTriggerMask, fddTriggerMask, uint8_t);

// FDD derived quantities
DECLARE_SOA_COLUMN(FDDChAmpl, fddChAmpl, std::vector<float>); //! FDD channel amplitudes, vector idx = ch ID
DECLARE_SOA_COLUMN(FDDTotAmplA, fddTotAmplA, float);
DECLARE_SOA_COLUMN(FDDTotAmplC, fddTotAmplC, float);

} // namespace fit

DECLARE_SOA_TABLE(FITExtras, "AOD", "FITEXTRA", //! Table with extra FIT information
                  fit::RunNumber, fit::GlobalBC, fit::CTPtriggerMask, fit::CTPinputMask,
                  fit::Timestamp,
                  fit::Sel8, fit::HasFT0, fit::HasFV0, fit::HasFDD,
                  fit::BCId, fit::PosX, fit::PosY, fit::PosZ, fit::Flags, fit::NumContrib, fit::CollisionTime, fit::CollisionTimeRes,
                  fit::FT0BCId, fit::FT0AmplitudeA, fit::FT0ChannelA, fit::FT0AmplitudeC, fit::FT0ChannelC,
                  fit::FT0TimeA, fit::FT0TimeC, fit::FT0TriggerMask, fit::FT0PosZ, fit::FT0CollTime, fit::FT0SumAmpA, fit::FT0SumAmpC,
                  fit::T0ACorrected, fit::T0CCorrected, fit::T0AC, fit::T0Resolution,
                  fit::FT0ChAmpl, fit::FT0TotAmplA, fit::FT0TotAmplC,
                  fit::FV0BCId, fit::FV0Amplitude, fit::FV0Channel, fit::FV0Time, fit::FV0TriggerMask,
                  fit::FV0ChAmpl, fit::FV0TotAmpl,
                  fit::FDDBCId, fit::FDDChargeA, fit::FDDChargeC, fit::FDDTimeA, fit::FDDTimeC, fit::FDDTriggerMask,
                  fit::FDDChAmpl, fit::FDDTotAmplA, fit::FDDTotAmplC);

using FITExtra = FITExtras::iterator;

} // namespace o2::aod

#endif // COMMON_DATAMODEL_FITEXTRA_H_