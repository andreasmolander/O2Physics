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

namespace o2::aod
{
namespace fit
{
// Quantities copied straight from AOD
// TODO: do we need them here?
DECLARE_SOA_COLUMN(PV, pv, float);                     //! Primary vertex position in cm (o2::aod::‌collision::PosZ)
DECLARE_SOA_COLUMN(NContrib, nContrib, int);           //! Number of contributors to primary vertex (o2::aod::‌collision::NumContrib)
DECLARE_SOA_COLUMN(FT0TimeA, ft0timeA, float);         //! FT0-A average time in ns (o2::aod::ft0::TimeA)
DECLARE_SOA_COLUMN(FT0TimeC, ft0timeC, float);         //! FT0-C average time in ns (o2::aod::ft0::TimeC)
DECLARE_SOA_COLUMN(FT0TimeACorr, ft0timeACorr, float); //! FT0-A average time in ns corrected PV (o2::aod::ft0::T0ACorrected)
DECLARE_SOA_COLUMN(FT0TimeCCorr, ft0timeCCorr, float); //! FT0-C average time in ns corrected PV (o2::aod::ft0::T0CCorrected)
DECLARE_SOA_COLUMN(FT0Time, ft0time, float);           //! FT0 collision time in ns (o2::aod::ft0::CollTime)
DECLARE_SOA_COLUMN(FT0TimeRes, ft0timeRes, float);     //! FT0 collision time resolution in ns (o2::aod::ft0::T0Resolution)
DECLARE_SOA_COLUMN(FT0Vtx, ft0vtx, float);             //! FT0 vertex in cm (o2::aod::ft0::PosZ)
DECLARE_SOA_COLUMN(FV0Time, fv0time, float);           //! FV0 average time in ns (o2::aod::fv0a::Time)
DECLARE_SOA_COLUMN(FDDTimeA, fddtimeA, float);         //! FDD-A average time in ns (o2::aod::fdd::TimeA)
DECLARE_SOA_COLUMN(FDDTimeC, fddtimeC, float);         //! FDD-C average time in ns (o2::aod::fdd::TimeC)

// Derived quantities

// Event selection conditions straigt from AOD
// TODO: do we need them here?
DECLARE_SOA_COLUMN(Sel8, sel8, bool);                  //! (o2::aod::evsel::Sel8)
DECLARE_SOA_COLUMN(HasFT0, hasFT0, bool);              //! (o2::aod::collision::has_foundFT0())
DECLARE_SOA_COLUMN(HasFV0, hasFV0, bool);              //! (o2::aod::collision::has_foundFV0())
DECLARE_SOA_COLUMN(HasFDD, hasFDD, bool);              //! (o2::aod::collision::has_foundFDD())
DECLARE_SOA_COLUMN(FT0Triggers, ft0Triggers, uint8_t); //! FT0 trigger mask (o2::aod::ft0::TriggerMask)
DECLARE_SOA_COLUMN(FV0Triggers, fv0Triggers, uint8_t); //! FV0 trigger mask (o2::aod::fv0a::TriggerMask)
DECLARE_SOA_COLUMN(FDDTriggers, fddTriggers, uint8_t); //! FDD trigger mask (o2::aod::fdd::TriggerMask)
} // namespace fit

DECLARE_SOA_TABLE(FITExtras, "AOD", "FITEXTRA", //! Table with extra FIT information
                  fit::Sel8, fit::HasFT0, fit::HasFV0, fit::HasFDD,
                  fit::FT0Triggers, fit::FV0Triggers, fit::FDDTriggers,
                  fit::PV, fit::NContrib,
                  fit::FT0TimeA, fit::FT0TimeC, fit::FT0TimeACorr, fit::FT0TimeCCorr,
                  fit::FT0Time, fit::FT0TimeRes, fit::FT0Vtx,
                  fit::FV0Time, fit::FDDTimeA, fit::FDDTimeC);

using FITExtra = FITExtras::iterator;

} // namespace o2::aod

#endif // COMMON_DATAMODEL_FITEXTRA_H_