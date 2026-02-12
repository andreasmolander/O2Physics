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

/// \file   qaFIT.cxx
/// \author Andreas Molander andreas.molander@cern.ch
/// \brief  FIT QA

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/FT0Corrected.h"

#include <DataFormatsFIT/Triggers.h>
#include <DataFormatsFDD/RecPoint.h>
#include <DataFormatsFT0/RecPoints.h>
#include <DataFormatsFV0/RecPoints.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/InitContext.h>
#include <Framework/Output.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <TH1F.h>
#include <TH2.h>
#include <TH2F.h>
#include <TH3.h>
#include <TH3F.h>

#include <array>
#include <bitset>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;

struct fitQa {

  /* 
  General comments:
  - One has to be a bit careful with axis ranges and binning
    as well as the amount of event selection conditions to consider,
    especially for the 3D histograms. ROOT object buffers can easily
    exceed the 1GB limit.

    TODO:
    - finer resolution for time res plots
  */

  /* Constants - TODO: get from somehwere else, don't hardcode here? */

  static constexpr int nChFT0 = 208; ///< Number of FT0 channels
  static constexpr int nChFT0A = 96; ///< Number of FT0A channels (-> number of FT0C channels = nChFT0 - nChFT0A)
  static constexpr int nChFV0 = 48;  ///< Number of FV0 channels
  static constexpr int nChFDD = 16;  ///< Number of FDD channels
  static constexpr int nADC = 4096;  ///< Number of ADC channels

  /* Common histogram properties - TODO: make configurable */

  static constexpr int nBinsNContrib = 10000;
  static constexpr float nContribMin = -0.;
  static constexpr float nContribMax = 10000.;

  // Time axis limits for FT0 and FV0
  // Times are in ns, these give a resolution of 20 ps per bin
  static constexpr int nBinsT = 500;
  static constexpr float tMin = -5.;
  static constexpr float tMax = 5.;

  // Time axis limits for FDD (it has larger spread)
  // Times are in ns, these give a resolution of 20 ps per bin
  static constexpr int nBinsTFDD = 2000;
  static constexpr float tMinFDD = -20.;
  static constexpr float tMaxFDD = 20.;

  // Axis limits for histograms from which we deduce collision time for FT0 and FV0
  // Times are in ns, these give a resolution of 1 ps per bin
  static constexpr int nBinsTRes = 2000;
  static constexpr float tResMin = -1.;
  static constexpr float tResMax = 1.;

  // Axis limits for histograms from which we deduce collision time for FDD (it has larger spread)
  // Times are in ns, these give a resolution of 2 ps per bin
  static constexpr int nBinsTResFDD = 4000;
  static constexpr float tResMinFDD = -4.;
  static constexpr float tResMaxFDD = 4.;

  // Vertex axis limits
  // Vertex is in cm and these give a resolution of 0.2 cm per bin
  static constexpr int nBinsVtx = 300;
  static constexpr float vtxMin = -30.;
  static constexpr float vtxMax = 30.;

  // Vertex axis limits for FDD (it has larger spread)
  // Vertex is in cm and these give a resolution of 0.2 cm per bin
  static constexpr int nBinsVtxFDD = 600;
  static constexpr float vtxMinFDD = -60.;
  static constexpr float vtxMaxFDD = 60.;

  static constexpr int nBinsTotAmpl = 30000;
  static constexpr float totAmplMin = 0.;
  static constexpr float totAmplMax = 300000.;

  /* Helper functions */
  static float cm2ns(float cm) { return cm / o2::constants::physics::LightSpeedCm2NS; }
  static float ns2cm(float ns) { return ns * o2::constants::physics::LightSpeedCm2NS; }
  static float collTime(float tA, float tC) { return 0.5f * (tA + tC); }
  static float vtxTime(float tA, float tC) { return 0.5f * (tC - tA); }

  /* Variables storing AO2D quantities and conditions, these are set in each processing step.
     NOTE: if quantities and conditions are added, they need to be reset in
     resetVariables() */
  
  // Quantities
  float pv = -200.f;                                             ///< Primary vertex position in cm (o2::aod::‌collision::PosZ)
  int nContrib = 0;                                              ///< Number of contributors to primary vertex (o2::aod::‌collision::NumContrib)
  float ft0timeA = o2::ft0::RecPoints::sDummyCollissionTime;     ///< FT0-A average time in ns (o2::aod::ft0::TimeA)
  float ft0timeC = o2::ft0::RecPoints::sDummyCollissionTime;     ///< FT0-A average time in ns (o2::aod::ft0::TimeA)
  float ft0timeACorr = o2::ft0::RecPoints::sDummyCollissionTime; ///< FT0-A average time in ns corrected PV (o2::aod::ft0::T0ACorrected)
  float ft0timeCCorr = o2::ft0::RecPoints::sDummyCollissionTime; ///< FT0-A average time in ns corrected PV (o2::aod::ft0::T0ACorrected)
  float ft0time = o2::ft0::RecPoints::sDummyCollissionTime;      ///< FT0 collision time in ns (o2::aod::ft0::CollTime)
  float ft0timeRes = o2::ft0::RecPoints::sDummyCollissionTime;   ///< FT0 collision time resolution in ns (o2::aod::ft0::T0Resolution)
  float ft0vtx = -200.f;                                         ///< FT0 vertex in ns (o2::aod::ft0::PosZ)
  float fv0time = o2::fv0::RecPoints::sDummyCollissionTime;      ///< FV0 average time in ns (o2::aod::fv0a::Time)
  float fddtimeA = o2::fdd::RecPoint::sDummyCollissionTime;      ///< FDD-A average time in ns (o2::aod::fdd::TimeA)
  float fddtimeC = o2::fdd::RecPoint::sDummyCollissionTime;      ///< FDD-C average time in ns (o2::aod::fdd::TimeC)
  
  std::array<float, nChFT0> ft0ChAmpl{}; ///< FT0 channel amplitudes (o2::aod::ft0::AmplitudeA and o2::aod::ft0::AmplitudeC)
  std::array<float, nChFV0> fv0ChAmpl{}; ///< FV0 channel amplitudes (o2::aod::fv0::Amplitude)
  std::array<float, nChFDD> fddChAmpl{}; ///< FDD channel amplitudes (o2::aod::fdd::ChargeA and o2::aod::fdd::ChargeC)

  // Derived quantities (not directly from AO2D)
  float ft0TotAmpl = 0;
  float ft0TotAmplA = 0;
  float ft0TotAmplC = 0;
  float fv0TotAmpl = 0;
  float fddTotAmpl = 0;
  float fddTotAmplA = 0;
  float fddTotAmplC = 0;

  // Event selection conditions
  bool isSel8 = false;            ///< (o2::aod::evsel::Sel8)
  bool hasFT0 = false;            ///< (o2::aod::collision::has_foundFT0())
  bool hasFV0 = false;            ///< (o2::aod::collision::has_foundFV0())
  bool hasFDD = false;            ///< (o2::aod::collision::has_foundFDD())
  std::bitset<8> ft0Triggers;     ///< FT0 trigger bits (o2::aod::ft0::TriggerMask)
  std::bitset<8> fv0Triggers;     ///< FV0 trigger bits (o2::aod::fv0::TriggerMask)
  std::bitset<8> fddTriggers;     ///< FDD trigger bits (o2::aod::fdd::TriggerMask)

  /// Helper struct for event selection conditions
  struct Condition {
    std::string name;           ///< Condition name - TODO: maybe not needed
    std::string title;          ///< Condition title
    std::function<bool()> eval; ///< How to evaluate the condition
  };

  std::vector<Condition> conditions{};  //< Event selection conditions to be considered
  std::vector<Condition> ft0Conditions{}; //< FT0 related event selection conditions to be considered 
  std::vector<Condition> fv0Conditions{}; //< FT0 related event selection conditions to be considered 
  std::vector<Condition> fddConditions{}; //< FT0 related event selection conditions to be considered 
  
  /* Output objects
     NOTE1: 1D Quantities are stored in 2D hists, with event selection conditions on the Y-axis.
     I.e. there's one version of the 1D histogram per Y bin (Condition).
     Similar for 2D Quantities, they are stored in 3D histograms. 
    
    NOTE2: There's a limit of number of task struct members, so some OutputObj's are disabled for now
  */

  /* 1D */

  // General
  OutputObj<TH1F> ooStats{"Stats"}; ///< Event selection statistics

  // Collision
  OutputObj<TH2F> ooPV{"collPV"};             ///< Primary vertex (cm)
  OutputObj<TH2F> ooPVns{"collPVns"};         ///< Primary vertex (ns)
  OutputObj<TH2F> ooNcontrib{"collNcontrib"}; ///< Number of contributors
  OutputObj<TH2F> ooNcontribFT0{"collNcontribFT0"}; ///< Number of contributors, FT0 triggers
  OutputObj<TH2F> ooNcontribFV0{"collNcontribFV0"}; ///< Number of contributors, FV0 triggers
  OutputObj<TH2F> ooNcontribFDD{"collNcontribFDD"}; ///< Number of contributors, FDD triggers
  
  // FT0
  OutputObj<TH2F> ooFT0TimeA{"FT0TimeA"};         ///< FT0A average time (ns)
  OutputObj<TH2F> ooFT0TimeC{"FT0TimeC"};         ///< FT0C average time (ns)
  OutputObj<TH2F> ooFT0TimeACorr{"FT0TimeACorr"}; ///< PV corrected FT0A average time (ns)
  OutputObj<TH2F> ooFT0TimeCCorr{"FT0TimeCCorr"}; ///< PV corrected FT0C average time (ns)
  OutputObj<TH2F> ooFT0Time{"FT0Time"};           ///< FT0 collision time (ns)
  OutputObj<TH2F> ooFT0TimeRes{"FT0TimeRes"};     ///< FT0 collision time resolution (ns)
  OutputObj<TH2F> ooFT0Vtx{"FT0Vtx"};             ///< FT0 vertex (cm)
  OutputObj<TH2F> ooFT0VtxNS{"FT0VtxNS"};         ///< FT0 vertex (ns)
  OutputObj<TH2F> ooFT0TotAmpl{"FT0TotAmpl"};     ///< FT0 total amplitude (ADC)
  OutputObj<TH2F> ooFT0TotAmplA{"FT0TotAmplA"};   ///< FT0A total amplitude (ADC)
  OutputObj<TH2F> ooFT0TotAmplC{"FT0TotAmplC"};   ///< FT0C total amplitude (ADC)

  // FV0
  OutputObj<TH2F> ooFV0Time{"FV0Time"};                     ///< FV0 average time (ns)
  OutputObj<TH2F> ooFV0TimeCorr{"FV0TimeCorr"};             ///< PV corrected FV0 average time (ns)
  OutputObj<TH2F> ooPVFV0FT0CVtxDiffNS{"PVFV0FT0CVtxDiff"}; ///< PV - FV0-FT0C vertex (ns)
  OutputObj<TH2F> ooFV0TotAmpl{"FV0TotAmpl"};               ///< FV0 total amplitude (ADC)

  // FDD
  OutputObj<TH2F> ooFDDTimeA{"FDDTimeA"};             ///< FDDA average time (ns)
  OutputObj<TH2F> ooFDDTimeC{"FDDTimeC"};             ///< FDDC average time (ns)
  OutputObj<TH2F> ooFDDTimeACorr{"FDDTimeACorr"};     ///< PV corrected FDDA average time (ns)
  OutputObj<TH2F> ooFDDTimeCCorr{"FDDTimeCCorr"};     ///< PV corrected FDDC average time (ns)
  OutputObj<TH2F> ooFDDTime{"FDDTime"};               ///< FDD collision time (ns)
  OutputObj<TH2F> ooFDDVtx{"FDDVtx"};                 ///< FDD vertex (cm)
  OutputObj<TH2F> ooFDDVtxNS{"FDDVtxNS"};             ///< FDD vertex (ns)
  OutputObj<TH2F> ooPVFDDVtxDiffNS{"PVFDDVtxDiffNS"}; ///< PV - FDD vertex (ns)
  OutputObj<TH2F> ooFDDTotAmpl{"FDDTotAmpl"};         ///< FDD total amplitude (ADC)
  OutputObj<TH2F> ooFDDTotAmplA{"FDDTotAmplA"};       ///< FDD-A total amplitude (ADC)
  OutputObj<TH2F> ooFDDTotAmplC{"FDDTotAmplC"};       ///< FDD-C total amplitude (ADC)

  // FT0, FV0
  OutputObj<TH2F> ooFT0TimeFV0TimeDiff{"FT0TimeFV0TimeDiff"};   ///< FT0 collision time - FV0 average time (ns)
  OutputObj<TH2F> ooFT0TimeAFV0TimeDiff{"FT0TimeAFV0TimeDiff"}; ///< FT0A average time - FV0 average time (ns)

  // FT0, FDD
  OutputObj<TH2F> ooFT0TimeFDDTimeDiff{"FT0TimeFDDTimeDiff"}; ///< FT0 collision time - FDD collision time (ns)
  OutputObj<TH2F> ooFT0VtxFDDVtxDiffNS{"FT0VtxFDDVtxDiffNS"}; ///< FT0 vertex - FDD vertex (ns)

  // FV0, FDD
  OutputObj<TH2F> ooFV0TimeFDDTimeDiff{"FV0TimeFDDTimeDiff"}; ///< FV0 average time - FDD collision time (ns)

  /* 2D */
  
  // FT0
  OutputObj<TH3F> ooFT0TimeVsFT0Vtx{"FT0TimeVsFT0Vtx"};               ///< FT0 collision time vs FT0 vertex
  OutputObj<TH3F> ooPVvsFT0Vtx{"PVvsFT0Vtx"};                         ///< PV vs FT0 vertex
  OutputObj<TH3F> ooFT0TimeResVsNContrib{"FT0TimeResVsNContrib"};     ///< FT0 time resolution vs number of contributors
  OutputObj<TH3F> ooFT0TimeResVsFT0TotAmpl{"FT0TimeResVsFT0TotAmpl"}; ///< FT0 time resolution vs FT0 total amplitude

  // FV0
  OutputObj<TH3F> ooPVFV0FT0CVtxDiffNSVsNContrib{"PVFV0FT0CVtxDiffNSvsNContrib"};     ///< PV - FV0-FT0C vertex vs number of contributors
  OutputObj<TH3F> ooPVFV0FT0CVtxDiffNSVsFV0TotAmpl{"PVFV0FT0CVtxDiffNSvsFV0TotAmpl"}; ///< PV - FV0-FT0C vertex vs FV0 total amplitude

  // FDD
  OutputObj<TH3F> ooFDDTimeVsFDDVtx{"FDDTimeVsFDDVtx"};                       ///< FDD collision time vs FDD vertex
  OutputObj<TH3F> ooPVvsFDDVtx{"PVvsFDDVtx"};                                 ///< PV vs FDD vertex
  OutputObj<TH3F> ooPVFDDVtxDiffNSVsNContrib{"PVFDDVtxDiffNSvsNContrib"};     ///< PV - FDD vertex vs number of contributors
  OutputObj<TH3F> ooPVFDDVtxDiffNSVsFDDTotAmpl{"PVFDDVtxDiffNSvsFDDTotAmpl"}; ///< PV - FDD vertex vs FDD total amplitude

  // FT0, FV0
  OutputObj<TH3F> ooFT0TimeAFV0TimeDiffVsNContrib{"FT0TimeAFV0TimeDiffVsNContrib"};     ///< FT0A average time - FV0 average time vs number of contributors
  OutputObj<TH3F> ooFT0TimeAFV0TimeDiffVsFV0TotAmpl{"FT0TimeAFV0TimeDiffVsFV0TotAmpl"}; ///< FT0A average time - FV0 average time vs FV0 total amplitude
  OutputObj<TH3F> ooFT0TotAmplVsFV0TotAmpl{"FT0TotAmplVsFV0TotAmpl"};                   ///< FT0 total amplitude vs FV0 total amplitude
  OutputObj<TH3F> ooFT0TotAmplAVsFV0TotAmpl{"FT0TotAmplAVsFV0TotAmpl"};                 ///< FT0A total amplitude vs FV0 total amplitude
  OutputObj<TH3F> ooFT0TotAmplCVsFV0TotAmpl{"FT0TotAmplCVsFV0TotAmpl"};                 ///< FT0C total amplitude vs FV0 total amplitude

  // FT0, FDD
  OutputObj<TH3F> ooFT0TotAmplVsFDDTotAmpl{"FT0TotAmplVsFDDTotAmpl"};   ///< FT0 total amplitude vs FDD total amplitude
  // OutputObj<TH3F> ooFT0TotAmplVsFDDTotAmplA{"FT0TotAmplVsFDDTotAmplA"}; ///< FT0 total amplitude vs FDD-A total amplitude
  // OutputObj<TH3F> ooFT0TotAmplVsFDDTotAmplC{"FT0TotAmplVsFDDTotAmplC"}; ///< FT0 total amplitude vs FDD-C total amplitude

  // OutputObj<TH3F> ooFT0TotAmplAVsFDDTotAmpl{"FT0TotAmplAVsFDDTotAmpl"};   ///< FT0A total amplitude vs FDD total amplitude
  // OutputObj<TH3F> ooFT0TotAmplAVsFDDTotAmplA{"FT0TotAmplAVsFDDTotAmplA"}; ///< FT0A total amplitude vs FDD-A total amplitude
  // OutputObj<TH3F> ooFT0TotAmplAVsFDDTotAmplC{"FT0TotAmplAVsFDDTotAmplC"}; ///< FT0A total amplitude vs FDD-C total amplitude
  
  // OutputObj<TH3F> ooFT0TotAmplCVsFDDTotAmpl{"FT0TotAmplCVsFDDTotAmpl"};   ///< FT0C total amplitude vs FDD total amplitude
  // OutputObj<TH3F> ooFT0TotAmplCVsFDDTotAmplA{"FT0TotAmplCVsFDDTotAmplA"}; ///< FT0C total amplitude vs FDD-A total amplitude
  // OutputObj<TH3F> ooFT0TotAmplCVsFDDTotAmplC{"FT0TotAmplCVsFDDTotAmplC"}; ///< FT0C total amplitude vs FDD-C total amplitude

  // // FV0, FDD
  OutputObj<TH3F> ooFV0TotAmplVsFDDTotAmpl{"FV0TotAmplVsFDDTotAmpl"};   ///< FV0 total amplitude vs FDD total amplitude
  // OutputObj<TH3F> ooFV0TotAmplVsFDDTotAmplA{"FV0TotAmplVsFDDTotAmplA"}; ///< FV0 total amplitude vs FDD-A total amplitude
  // OutputObj<TH3F> ooFV0TotAmplVsFDDTotAmplC{"FV0TotAmplVsFDDTotAmplC"}; ///< FV0 total amplitude vs FDD-C total amplitude

  /* 2D Quantities per X bin */

  OutputObj<TH3F> ooFT0AmplPerCh{"FT0AmplPerCh"}; ///< FT0 amplitude per channel
  OutputObj<TH3F> ooFV0AmplPerCh{"FV0AmplPerCh"}; ///< FV0 amplitude per channel
  OutputObj<TH3F> ooFDDAmplPerCh{"FDDAmplPerCh"}; ///< FDD amplitude per channel 

  /* Mapping of OutputObj pointers to functions that compute the quantity values */

  /// 1D quantities (one value per collision)
  std::unordered_map<OutputObj<TH2F>*, std::function<float()>> objs;
  std::unordered_map<OutputObj<TH2F>*, std::function<float()>> objsft0;
  std::unordered_map<OutputObj<TH2F>*, std::function<float()>> objsfv0;
  std::unordered_map<OutputObj<TH2F>*, std::function<float()>> objsfdd;
  /// 2D quantities (one value per collision, e.g. collision time vs vertex)
  std::unordered_map<OutputObj<TH3F>*, std::function<std::pair<float, float>()>> objs2D;
  /// 2D quantities per X bin (one value per X bin per collision, e.g. amplitude per channel)
  std::unordered_map<OutputObj<TH3F>*, std::function<float(int)>> objs2DPerXBin;

  /// Resetting all AO2D quantities and conditions to default values
  void resetVariables() {
    pv = -200.f;
    nContrib = 0;
    ft0timeA = o2::ft0::RecPoints::sDummyCollissionTime;
    ft0timeC = o2::ft0::RecPoints::sDummyCollissionTime;
    ft0timeACorr = o2::ft0::RecPoints::sDummyCollissionTime;
    ft0timeCCorr = o2::ft0::RecPoints::sDummyCollissionTime;
    ft0time = o2::ft0::RecPoints::sDummyCollissionTime;
    ft0timeRes = o2::ft0::RecPoints::sDummyCollissionTime;
    ft0vtx = -200.f;
    fv0time = o2::fv0::RecPoints::sDummyCollissionTime;
    fddtimeA = o2::fdd::RecPoint::sDummyCollissionTime;
    fddtimeC = o2::fdd::RecPoint::sDummyCollissionTime;

    // Same default values as in ChannelDataFloat
    ft0ChAmpl.fill(-20000.f);
    fv0ChAmpl.fill(-20000.f);
    fddChAmpl.fill(-20000.f);

    ft0TotAmpl = 0;
    ft0TotAmplA = 0;
    ft0TotAmplC = 0;
    fv0TotAmpl = 0;
    fddTotAmpl = 0;
    fddTotAmplA = 0;
    fddTotAmplC = 0;

    isSel8 = false;
    hasFT0 = false;
    hasFV0 = false;
    hasFDD = false;
    ft0Triggers.reset();
    fv0Triggers.reset();
    fddTriggers.reset();
  }

  void init(InitContext&) {
    /* Init event selection conditions - TODO make configurable */

    conditions.push_back({"All", "all collisions", [&]() { return true; }});
    // conditions.push_back({"sel8", "sel8", [&]() { return isSel8; }});
    // conditions.push_back({"HasFT0", "has FT0", [&]() { return hasFT0; }});
    // conditions.push_back({"HasFV0", "has FV0", [&]() { return hasFV0; }});
    // conditions.push_back({"HasFDD", "has FDD", [&]() { return hasFDD; }});
    conditions.push_back({"FT0VTX", "FT0 vertex", [&]() { return ft0Triggers.test(o2::fit::Triggers::bitVertex); }});
    // conditions.push_back({"FT0CE", "FT0 CE", [&]() { return isFT0CE; }});
    // conditions.push_back({"FT0SCE", "FT0 SCE", [&]() { return isFT0SCE; }});
    // conditions.push_back({"FV0ORA", "FV0 ORA", [&]() { return isFV0ORA; }});
    // conditions.push_back({"FV0CH", "FV0 CH", [&]() { return isFV0CH; }});
    // conditions.push_back({"FV0IN", "FV0 IN", [&]() { return isFV0IN; }});
    // conditions.push_back({"FDDVTX", "FDD vertex", [&]() { return isFDDVTX; }});
    conditions.push_back({"FT0VTXandFV0ORA", "FT0 vertex AND FV0 ORA", [&]() { return ft0Triggers.test(o2::fit::Triggers::bitVertex) && fv0Triggers.test(o2::fit::Triggers::bitA); }});
    // conditions.push_back({"FT0VTXandFV0CH", "FT0 vertex AND FV0 CH", [&]() { return isFT0VTX && isFV0CH; }});
    // conditions.push_back({"FT0VTXandFV0IN", "FT0 vertex AND FV0 IN", [&]() { return isFT0VTX && isFV0IN; }});
    conditions.push_back({"FT0VTXandFDDVTX", "FT0 vertex AND FDD vertex", [&]() { return ft0Triggers.test(o2::fit::Triggers::bitVertex) && fddTriggers.test(o2::fit::Triggers::bitVertex); }});

    ft0Conditions.push_back({"All", "all collisions", [&]() { return true; }});
    ft0Conditions.push_back({"FT0ORA", "FT0 ORA", [&]() { return ft0Triggers.test(o2::fit::Triggers::bitA); }});
    ft0Conditions.push_back({"FT0ORC", "FT0 ORC", [&]() { return ft0Triggers.test(o2::fit::Triggers::bitC); }});
    ft0Conditions.push_back({"FT0VTX", "FT0 VTX", [&]() { return ft0Triggers.test(o2::fit::Triggers::bitVertex); }});
    ft0Conditions.push_back({"FT0CE", "FT0 CE", [&]() { return ft0Triggers.test(o2::fit::Triggers::bitCen); }});
    ft0Conditions.push_back({"FT0SCE", "FT0 SCE", [&]() { return ft0Triggers.test(o2::fit::Triggers::bitSCen); }});
    ft0Conditions.push_back({"FT0ACTIVEA", "FT0 ACTIVE A", [&]() { return ft0Triggers.test(o2::fit::Triggers::bitLaser); }});
    ft0Conditions.push_back({"FT0ACTIVEC", "FT0 ACTIVE C", [&]() { return ft0Triggers.test(o2::fit::Triggers::bitOutputsAreBlocked); }});
    ft0Conditions.push_back({"FT0FLANGE", "FT0 FLANGE", [&]() { return ft0Triggers.test(o2::fit::Triggers::bitDataIsValid); }});
    ft0Conditions.push_back({"FT0MINBIAS", "FT0 MIN BIAS", [&]() { return ft0Triggers.test(o2::fit::Triggers::bitVertex) && (ft0Triggers.test(o2::fit::Triggers::bitCen) || ft0Triggers.test(o2::fit::Triggers::bitSCen)); }});

    fv0Conditions.push_back({"All", "all collisions", [&]() { return true; }});
    fv0Conditions.push_back({"FV0ORA", "FV0 ORA", [&]() { return fv0Triggers.test(o2::fit::Triggers::bitA); }});
    fv0Conditions.push_back({"FV0NCHAN", "FV0 NCHAN", [&]() { return fv0Triggers.test(o2::fit::Triggers::bitTrgNchan); }});
    fv0Conditions.push_back({"FV0CH", "FV0 CH", [&]() { return fv0Triggers.test(o2::fit::Triggers::bitTrgCharge); }});
    fv0Conditions.push_back({"FV0OUT", "FV0 OUT", [&]() { return fv0Triggers.test(o2::fit::Triggers::bitAOut); }});
    fv0Conditions.push_back({"FV0CH", "FV0 IN", [&]() { return fv0Triggers.test(o2::fit::Triggers::bitAIn); }});

    fddConditions.push_back({"All", "all collisions", [&]() { return true; }});
    fddConditions.push_back({"FDDORA", "FDD ORA", [&]() { return fddTriggers.test(o2::fit::Triggers::bitA); }});
    fddConditions.push_back({"FDDORC", "FDD ORC", [&]() { return fddTriggers.test(o2::fit::Triggers::bitC); }});
    fddConditions.push_back({"FDDVTX", "FDD vertex", [&]() { return fddTriggers.test(o2::fit::Triggers::bitVertex); }});
    fddConditions.push_back({"FDDCE", "FDD CE", [&]() { return fddTriggers.test(o2::fit::Triggers::bitCen); }});
    fddConditions.push_back({"FDDSCE", "FDD SCE", [&]() { return fddTriggers.test(o2::fit::Triggers::bitSCen); }});

    /* Init OutputObj's */

    ooStats.setObject(new TH1F("Stats", "Event selection statistics;;Collisions", conditions.size(), 0, conditions.size()));
    for (size_t c = 0; c < conditions.size(); c++) {
      ooStats.object->GetXaxis()->SetBinLabel(c + 1, conditions[c].title.c_str());
    }

    // Base quantities (AO2D data)

    ooPV.setObject(new TH2F(ooPV.label.c_str(), "Primary vertex;Primary vertex z position (cm)", nBinsVtx, vtxMin, vtxMax, conditions.size(), 0, conditions.size()));
    objs[&ooPV] = [&]() { return pv; };

    ooNcontrib.setObject(new TH2F(ooNcontrib.label.c_str(), "Number of contributors to primary vertex;Number of contributors to primary vertex", nBinsNContrib, nContribMin, nContribMax, conditions.size(), 0, conditions.size()));
    objs[&ooNcontrib] = [&]() { return nContrib; };

    ooNcontribFT0.setObject(new TH2F(ooNcontribFT0.label.c_str(), "Number of contributors to primary vertex;Number of contributors to primary vertex", nBinsNContrib, nContribMin, nContribMax, ft0Conditions.size(), 0, ft0Conditions.size()));
    objsft0[&ooNcontribFT0] = [&]() { return nContrib; };

    ooNcontribFV0.setObject(new TH2F(ooNcontribFV0.label.c_str(), "Number of contributors to primary vertex;Number of contributors to primary vertex", nBinsNContrib, nContribMin, nContribMax, fv0Conditions.size(), 0, fv0Conditions.size()));
    objsfv0[&ooNcontribFV0] = [&]() { return nContrib; };

    ooNcontribFDD.setObject(new TH2F(ooNcontribFDD.label.c_str(), "Number of contributors to primary vertex;Number of contributors to primary vertex", nBinsNContrib, nContribMin, nContribMax, fddConditions.size(), 0, fddConditions.size()));
    objsfdd[&ooNcontribFDD] = [&]() { return nContrib; };

    ooFT0TimeA.setObject(new TH2F(ooFT0TimeA.label.c_str(), "FT0A time;$\\langle t_{\\text{FT0A}} \\rangle \\text{ (ns)}$", nBinsT, tMin, tMax, conditions.size(), 0, conditions.size()));
    objs[&ooFT0TimeA] = [&]() { return ft0timeA; };

    ooFT0TimeC.setObject(new TH2F(ooFT0TimeC.label.c_str(), "FT0C time;$\\langle t_{\\text{FT0C}} \\rangle \\text{ (ns)}$", nBinsT, tMin, tMax, conditions.size(), 0, conditions.size()));
    objs[&ooFT0TimeC] = [&]() { return ft0timeC; };

    ooFT0TimeACorr.setObject(new TH2F(ooFT0TimeACorr.label.c_str(), "PV corrected FT0A time;$\\langle t_{\\text{FT0A}} \\rangle + \\text{PV} \\text{ (ns)}$", nBinsT, tMin, tMax, conditions.size(), 0, conditions.size()));
    objs[&ooFT0TimeACorr] = [&]() { return ft0timeACorr; };

    ooFT0TimeCCorr.setObject(new TH2F(ooFT0TimeCCorr.label.c_str(), "PV corrected FT0C time;$\\langle t_{\\text{FT0C}} \\rangle - \\text{PV} \\text{ (ns)}$", nBinsT, tMin, tMax, conditions.size(), 0, conditions.size()));
    objs[&ooFT0TimeCCorr] = [&]() { return ft0timeCCorr; };

    ooFT0Time.setObject(new TH2F(ooFT0Time.label.c_str(), "FT0 time;$(\\langle t_{\\text{FT0A}} \\rangle + \\langle t_{\\text{FT0C}} \\rangle)/2 \\text{ (ns)}$", nBinsT, tMin, tMax, conditions.size(), 0, conditions.size()));
    objs[&ooFT0Time] = [&]() { return ft0time; };

    ooFT0TimeRes.setObject(new TH2F(ooFT0TimeRes.label.c_str(), "FT0 time resolution;$\\text{PV} - (\\langle t_{\\text{FT0C}} \\rangle - \\langle t_{\\text{FT0A}} \\rangle)/2 \\text{ (ns)}$", nBinsTRes, tResMin, tResMax, conditions.size(), 0, conditions.size()));
    objs[&ooFT0TimeRes] = [&]() { return ft0timeRes; };

    ooFT0Vtx.setObject(new TH2F(ooFT0Vtx.label.c_str(), "FT0 vertex;$(\\langle t_{\\text{FT0C}} \\rangle - \\langle t_{\\text{FT0A}} \\rangle)/2 \\text{ (cm)}$", nBinsVtx, vtxMin, vtxMax, conditions.size(), 0, conditions.size()));
    objs[&ooFT0Vtx] = [&]() { return ft0vtx; };

    ooFV0Time.setObject(new TH2F(ooFV0Time.label.c_str(), "FV0 time;$\\langle t_{\\text{FV0}} \\rangle \\text{ (ns)}$", nBinsT, tMin, tMax, conditions.size(), 0, conditions.size()));
    objs[&ooFV0Time] = [&]() { return fv0time; };

    ooFDDTimeA.setObject(new TH2F(ooFDDTimeA.label.c_str(), "FDDA time;$\\langle t_{\\text{FDDA}} \\rangle \\text{ (ns)}$", nBinsTFDD, tMinFDD, tMaxFDD, conditions.size(), 0, conditions.size()));
    objs[&ooFDDTimeA] = [&]() { return fddtimeA; };

    ooFDDTimeC.setObject(new TH2F(ooFDDTimeC.label.c_str(), "FDDC time;$\\langle t_{\\text{FDDC}} \\rangle \\text{ (ns)}$", nBinsTFDD, tMinFDD, tMaxFDD, conditions.size(), 0, conditions.size()));
    objs[&ooFDDTimeC] = [&]() { return fddtimeC; };

    // TODO: fill multiplicity from mults tables and compare
    ooFT0TotAmpl.setObject(new TH2F(ooFT0TotAmpl.label.c_str(), "FT0 total amplitude;FT0 amplitude (ADC)", nBinsTotAmpl, totAmplMin, totAmplMax, conditions.size(), 0, conditions.size()));
    objs[&ooFT0TotAmpl] = [&]() { return ft0TotAmpl; };

    ooFT0TotAmplA.setObject(new TH2F(ooFT0TotAmplA.label.c_str(), "FT0A total amplitude;FT0A amplitude (ADC)", nBinsTotAmpl, totAmplMin, totAmplMax, conditions.size(), 0, conditions.size()));
    objs[&ooFT0TotAmplA] = [&]() { return ft0TotAmplA; };

    ooFT0TotAmplC.setObject(new TH2F(ooFT0TotAmplC.label.c_str(), "FT0C total amplitude;FT0C amplitude (ADC)", nBinsTotAmpl, totAmplMin, totAmplMax, conditions.size(), 0, conditions.size()));
    objs[&ooFT0TotAmplC] = [&]() { return ft0TotAmplC; };

    ooFV0TotAmpl.setObject(new TH2F(ooFV0TotAmpl.label.c_str(), "FV0 total amplitude;FV0 amplitude (ADC)", nBinsTotAmpl, totAmplMin, totAmplMax, conditions.size(), 0, conditions.size()));
    objs[&ooFV0TotAmpl] = [&]() { return fv0TotAmpl; };

    ooFDDTotAmpl.setObject(new TH2F(ooFDDTotAmpl.label.c_str(), "FDD total amplitude;FDD amplitude (ADC)", nBinsTotAmpl, totAmplMin, totAmplMax, conditions.size(), 0, conditions.size()));
    objs[&ooFDDTotAmpl] = [&]() { return fddTotAmpl; };

    ooFDDTotAmplA.setObject(new TH2F(ooFDDTotAmplA.label.c_str(), "FDD-A total amplitude;FDD-A amplitude (ADC)", nBinsTotAmpl, totAmplMin, totAmplMax, conditions.size(), 0, conditions.size()));
    objs[&ooFDDTotAmplA] = [&]() { return fddTotAmplA; };

    ooFDDTotAmplC.setObject(new TH2F(ooFDDTotAmplC.label.c_str(), "FDD-C total amplitude;FDD-C amplitude (ADC)", nBinsTotAmpl, totAmplMin, totAmplMax, conditions.size(), 0, conditions.size()));
    objs[&ooFDDTotAmplC] = [&]() { return fddTotAmplC; };

    // Derived quantities
    // TODO: some should maybe be calculated in AO2D tables?

    ooPVns.setObject(new TH2F(ooPVns.label.c_str(), "Primary vertex;Primary vertex z position (ns)", nBinsT, tMin, tMax, conditions.size(), 0, conditions.size()));
    objs[&ooPVns] = [&]() { return cm2ns(objs[&ooPV]()); };

    ooFT0VtxNS.setObject(new TH2F(ooFT0VtxNS.label.c_str(), "FT0 vertex;$(\\langle t_{\\text{FT0C}} \\rangle - \\langle t_{\\text{FT0A}} \\rangle)/2 \\text{ (ns)}$", nBinsT, tMin, tMax, conditions.size(), 0, conditions.size()));
    objs[&ooFT0VtxNS] = [&]() { return cm2ns(objs[&ooFT0Vtx]()); };

    ooFV0TimeCorr.setObject(new TH2F(ooFV0TimeCorr.label.c_str(), "PV corrected FV0 time;$\\langle t_{\\text{FV0}} \\rangle + \\text{PV} \\text{ (ns)}$", nBinsT, tMin, tMax, conditions.size(), 0, conditions.size()));
    objs[&ooFV0TimeCorr] = [&]() { return objs[&ooFV0Time]() + objs[&ooPVns](); };

    ooFDDTime.setObject(new TH2F(ooFDDTime.label.c_str(), "FDD time;$(\\langle t_{\\text{FDDA}} \\rangle + \\langle t_{\\text{FDDC}} \\rangle)/2 \\text{ (ns)}$", nBinsTFDD, tMinFDD, tMaxFDD, conditions.size(), 0, conditions.size()));
    objs[&ooFDDTime] = [&]() { return collTime(objs[&ooFDDTimeA](), objs[&ooFDDTimeC]()); };

    ooFDDVtx.setObject(new TH2F(ooFDDVtx.label.c_str(), "FDD vertex;$(\\langle t_{\\text{FDDC}} \\rangle - \\langle t_{\\text{FDDA}} \\rangle)/2 \\text{ (cm)}$", nBinsVtxFDD, vtxMinFDD, vtxMaxFDD, conditions.size(), 0, conditions.size()));
    objs[&ooFDDVtx] = [&]() { return ns2cm(objs[&ooFDDVtxNS]()); };

    ooFDDVtxNS.setObject(new TH2F(ooFDDVtxNS.label.c_str(), "FDD vertex;$(\\langle t_{\\text{FDDC}} \\rangle - \\langle t_{\\text{FDDA}} \\rangle)/2 \\text{ (ns)}$", nBinsTFDD, tMinFDD, tMaxFDD, conditions.size(), 0, conditions.size()));
    objs[&ooFDDVtxNS] = [&]() { return vtxTime(objs[&ooFDDTimeA](), objs[&ooFDDTimeC]()); };

    ooFDDTimeACorr.setObject(new TH2F(ooFDDTimeACorr.label.c_str(), "PV corrected FDDA time;$\\langle t_{\\text{FDDA}} \\rangle + \\text{PV} \\text{ (ns)}$", nBinsTFDD, tMinFDD, tMaxFDD, conditions.size(), 0, conditions.size()));
    objs[&ooFDDTimeACorr] = [&]() { return objs[&ooFDDTimeA]() + objs[&ooPVns](); };

    ooFDDTimeCCorr.setObject(new TH2F(ooFDDTimeCCorr.label.c_str(), "PV corrected FDDC time;$\\langle t_{\\text{FDDC}} \\rangle - \\text{PV} \\text{ (ns)}$", nBinsTFDD, tMinFDD, tMaxFDD, conditions.size(), 0, conditions.size()));
    objs[&ooFDDTimeCCorr] = [&]() { return objs[&ooFDDTimeC]() - objs[&ooPVns](); };

    ooFT0TimeFV0TimeDiff.setObject(new TH2F(ooFT0TimeFV0TimeDiff.label.c_str(), "FT0 time - FV0 time;$(\\langle t_{\\text{FT0A}} \\rangle + \\langle t_{\\text{FT0C}} \\rangle)/2 - \\langle t_{\\text{FV0}} \\rangle \\text{ (ns)}$", nBinsTRes, tResMin, tResMax, conditions.size(), 0, conditions.size()));
    objs[&ooFT0TimeFV0TimeDiff] = [&]() { return objs[&ooFT0Time]() - objs[&ooFV0Time](); };

    ooFT0TimeAFV0TimeDiff.setObject(new TH2F(ooFT0TimeAFV0TimeDiff.label.c_str(), "FT0A time - FV0 time;$\\langle t_{\\text{FT0A}} \\rangle - \\langle t_{\\text{FV0}} \\rangle \\text{ (ns)}$", nBinsTRes, tResMin, tResMax, conditions.size(), 0, conditions.size()));
    objs[&ooFT0TimeAFV0TimeDiff] = [&]() { return objs[&ooFT0TimeA]() - objs[&ooFV0Time](); };

    ooFT0TimeFDDTimeDiff.setObject(new TH2F(ooFT0TimeFDDTimeDiff.label.c_str(), "FT0 time - FDD time;$(\\langle t_{\\text{FT0A}} \\rangle + \\langle t_{\\text{FT0C}} \\rangle)/2 - (\\langle t_{\\text{FDDA}} \\rangle + \\langle t_{\\text{FDDC}} \\rangle)/2 \\text{ (ns)}$", nBinsTResFDD, tResMinFDD, tResMaxFDD, conditions.size(), 0, conditions.size()));
    objs[&ooFT0TimeFDDTimeDiff] = [&]() { return objs[&ooFT0Time]() - objs[&ooFDDTime](); };

    ooFV0TimeFDDTimeDiff.setObject(new TH2F(ooFV0TimeFDDTimeDiff.label.c_str(), "FV0 time - FDD time;$\\langle t_{\\text{FV0}} \\rangle - (\\langle t_{\\text{FDDA}} \\rangle + \\langle t_{\\text{FDDC}} \\rangle)/2 \\text{ (ns)}$", nBinsTResFDD, tResMinFDD, tResMaxFDD, conditions.size(), 0, conditions.size()));
    objs[&ooFV0TimeFDDTimeDiff] = [&]() { return objs[&ooFV0Time]() - objs[&ooFDDTime](); };

    // TODO: plot FV0-FT0C vertex separately?
    ooPVFV0FT0CVtxDiffNS.setObject(new TH2F(ooPVFV0FT0CVtxDiffNS.label.c_str(), "PV - FV0-FT0C vertex;$\\text{PV} - (\\langle t_{\\text{FT0C}} \\rangle - \\langle t_{\\text{FV0}} \\rangle)/2 \\text{ (ns)}$", nBinsTRes, tResMin, tResMax, conditions.size(), 0, conditions.size()));
    objs[&ooPVFV0FT0CVtxDiffNS] = [&]() { return objs[&ooPVns]() - vtxTime(objs[&ooFV0Time](), objs[&ooFT0TimeC]()); };

    ooPVFDDVtxDiffNS.setObject(new TH2F(ooPVFDDVtxDiffNS.label.c_str(), "PV - FDD vertex;$\\text{PV} - (\\langle t_{\\text{FDDC}} \\rangle - \\langle t_{\\text{FDDA}} \\rangle)/2 \\text{ (ns)}$", nBinsTResFDD, tResMinFDD, tResMaxFDD, conditions.size(), 0, conditions.size()));
    objs[&ooPVFDDVtxDiffNS] = [&]() { return objs[&ooPVns]() - objs[&ooFDDVtxNS](); };

    ooFT0VtxFDDVtxDiffNS.setObject(new TH2F(ooFT0VtxFDDVtxDiffNS.label.c_str(), "FT0 vertex - FDD vertex;$(\\langle t_{\\text{FT0C}} \\rangle - \\langle t_{\\text{FT0A}} \\rangle)/2 - (\\langle t_{\\text{FDDC}} \\rangle - \\langle t_{\\text{FDDA}} \\rangle)/2 \\text{ (ns)}$", nBinsTResFDD, tResMinFDD, tResMaxFDD, conditions.size(), 0, conditions.size()));
    objs[&ooFT0VtxFDDVtxDiffNS] = [&]() { return objs[&ooFT0VtxNS]() - objs[&ooFDDVtxNS](); };

    // 2D quantities

    ooFT0TimeVsFT0Vtx.setObject(new TH3F(ooFT0TimeVsFT0Vtx.label.c_str(), "FT0 time vs FT0 vertex;$(\\langle t_{\\text{FT0C}} \\rangle - \\langle t_{\\text{FT0A}} \\rangle)/2 \\text{ (cm)}$;$(\\langle t_{\\text{FT0A}} \\rangle + \\langle t_{\\text{FT0C}} \\rangle)/2 \\text{ (ns)}$", nBinsVtx, vtxMin, vtxMax, nBinsT, tMin, tMax, conditions.size(), 0, conditions.size()));
    objs2D[&ooFT0TimeVsFT0Vtx] = [&]() { return std::make_pair(objs[&ooFT0Vtx](), objs[&ooFT0Time]()); };

    ooPVvsFT0Vtx.setObject(new TH3F(ooPVvsFT0Vtx.label.c_str(), "PV vs FT0 vertex;$(\\langle t_{\\text{FT0C}} \\rangle - \\langle t_{\\text{FT0A}} \\rangle)/2 \\text{ (cm)}$;Primary vertex z position (cm)", nBinsVtx, vtxMin, vtxMax, nBinsVtx, vtxMin, vtxMax, conditions.size(), 0, conditions.size()));
    objs2D[&ooPVvsFT0Vtx] = [&]() { return std::make_pair(objs[&ooFT0Vtx](), objs[&ooPV]()); };

    ooFT0TimeResVsNContrib.setObject(new TH3F(ooFT0TimeResVsNContrib.label.c_str(), "FT0 time resolution vs number of contributors;Number of contributors to primary vertex;$\\text{PV} - (\\langle t_{\\text{FT0C}} \\rangle - \\langle t_{\\text{FT0A}} \\rangle)/2 \\text{ (ns)}$", nBinsNContrib / 10, nContribMin, nContribMax, nBinsTRes, tResMin, tResMax, conditions.size(), 0, conditions.size()));
    objs2D[&ooFT0TimeResVsNContrib] = [&]() { return std::make_pair(objs[&ooNcontrib](), objs[&ooFT0TimeRes]()); };

    ooFT0TimeResVsFT0TotAmpl.setObject(new TH3F(ooFT0TimeResVsFT0TotAmpl.label.c_str(), "FT0 time resolution vs FT0 total amplitude;FT0 total amplitude (ADC);$\\text{PV} - (\\langle t_{\\text{FT0C}} \\rangle - \\langle t_{\\text{FT0A}} \\rangle)/2 \\text{ (ns)}$", nBinsTotAmpl / 10, totAmplMin, totAmplMax, nBinsTRes, tResMin, tResMax, conditions.size(), 0, conditions.size()));
    objs2D[&ooFT0TimeResVsFT0TotAmpl] = [&]() { return std::make_pair(objs[&ooFT0TotAmpl](), objs[&ooFT0TimeRes]()); };

    ooPVFV0FT0CVtxDiffNSVsNContrib.setObject(new TH3F(ooPVFV0FT0CVtxDiffNSVsNContrib.label.c_str(), "PV - FV0-FT0C vertex vs number of contributors;Number of contributors to primary vertex;$\\text{PV} - (\\langle t_{\\text{FT0C}} \\rangle - \\langle t_{\\text{FV0}} \\rangle)/2 \\text{ (ns)}$", nBinsNContrib / 10, nContribMin, nContribMax, nBinsTRes, tResMin, tResMax, conditions.size(), 0, conditions.size()));
    objs2D[&ooPVFV0FT0CVtxDiffNSVsNContrib] = [&]() { return std::make_pair(objs[&ooNcontrib](), objs[&ooPVFV0FT0CVtxDiffNS]()); };

    ooPVFV0FT0CVtxDiffNSVsFV0TotAmpl.setObject(new TH3F(ooPVFV0FT0CVtxDiffNSVsFV0TotAmpl.label.c_str(), "PV - FV0-FT0C vertex vs FV0 total amplitude;FV0 total amplitude (ADC);$\\text{PV} - (\\langle t_{\\text{FT0C}} \\rangle - \\langle t_{\\text{FV0}} \\rangle)/2 \\text{ (ns)}$", nBinsTotAmpl / 10, totAmplMin, totAmplMax, nBinsTRes, tResMin, tResMax, conditions.size(), 0, conditions.size()));
    objs2D[&ooPVFV0FT0CVtxDiffNSVsFV0TotAmpl] = [&]() { return std::make_pair(objs[&ooFV0TotAmpl](), objs[&ooPVFV0FT0CVtxDiffNS]()); };

    ooFDDTimeVsFDDVtx.setObject(new TH3F(ooFDDTimeVsFDDVtx.label.c_str(), "FDD collision time vs FDD vertex;$(\\langle t_{\\text{FDDC}} \\rangle - \\langle t_{\\text{FDDA}} \\rangle)/2 \\text{ (cm)}$;$(\\langle t_{\\text{FDDA}} \\rangle + \\langle t_{\\text{FDDC}} \\rangle)/2 \\text{ (ns)}$", nBinsVtxFDD, vtxMinFDD, vtxMaxFDD, nBinsTFDD, tMinFDD, tMaxFDD, conditions.size(), 0, conditions.size()));
    objs2D[&ooFDDTimeVsFDDVtx] = [&]() { return std::make_pair(objs[&ooFDDVtx](), objs[&ooFDDTime]()); };

    ooPVvsFDDVtx.setObject(new TH3F(ooPVvsFDDVtx.label.c_str(), "PV vs FDD vertex;$(\\langle t_{\\text{FDDC}} \\rangle - \\langle t_{\\text{FDDA}} \\rangle)/2 \\text{ (cm)}$;Primary vertex z position (cm)", nBinsVtxFDD, vtxMinFDD, vtxMaxFDD, nBinsVtx, vtxMin, vtxMax, conditions.size(), 0, conditions.size()));
    objs2D[&ooPVvsFDDVtx] = [&]() { return std::make_pair(objs[&ooFDDVtx](), objs[&ooPV]()); };

    ooPVFDDVtxDiffNSVsNContrib.setObject(new TH3F(ooPVFDDVtxDiffNSVsNContrib.label.c_str(), "PV - FDD vertex vs number of contributors;Number of contributors to primary vertex;$\\text{PV} - (\\langle t_{\\text{FDDC}} \\rangle - \\langle t_{\\text{FDDA}} \\rangle)/2 \\text{ (ns)}$", nBinsNContrib / 10, nContribMin, nContribMax, nBinsTResFDD, tResMinFDD, tResMaxFDD, conditions.size(), 0, conditions.size()));
    objs2D[&ooPVFDDVtxDiffNSVsNContrib] = [&]() { return std::make_pair(objs[&ooNcontrib](), objs[&ooPVFDDVtxDiffNS]()); };

    // ooPVFDDVtxDiffNSVsFDDTotAmpl.setObject(new TH3F(ooPVFDDVtxDiffNSVsFDDTotAmpl.label.c_str(), "PV - FDD vertex vs FDD total amplitude;FDD total amplitude (ADC);$\\text{PV} - (\\langle t_{\\text{FDDC}} \\rangle - \\langle t_{\\text{FDDA}} \\rangle)/2 \\text{ (ns)}$", nBinsTotAmpl / 10, totAmplMin, totAmplMax, nBinsTResFDD, tResMinFDD, tResMaxFDD, conditions.size(), 0, conditions.size()));
    // objs2D[&ooPVFDDVtxDiffNSVsFDDTotAmpl] = [&]() { return std::make_pair(objs[&ooFDDTotAmpl](), objs[&ooPVFDDVtxDiffNS]()); };

    ooFT0TimeAFV0TimeDiffVsNContrib.setObject(new TH3F(ooFT0TimeAFV0TimeDiffVsNContrib.label.c_str(), "FT0A time - FV0 time vs number of contributors;Number of contributors to primary vertex;$\\langle t_{\\text{FT0A}} \\rangle - \\langle t_{\\text{FV0}} \\rangle \\text{ (ns)}$", nBinsNContrib / 10, nContribMin, nContribMax, nBinsTRes, tResMin, tResMax, conditions.size(), 0, conditions.size()));
    objs2D[&ooFT0TimeAFV0TimeDiffVsNContrib] = [&]() { return std::make_pair(objs[&ooNcontrib](), objs[&ooFT0TimeAFV0TimeDiff]()); };

    ooFT0TimeAFV0TimeDiffVsFV0TotAmpl.setObject(new TH3F(ooFT0TimeAFV0TimeDiffVsFV0TotAmpl.label.c_str(), "FT0A time - FV0 time vs FV0 total amplitude;FV0 total amplitude (ADC);$\\langle t_{\\text{FT0A}} \\rangle - \\langle t_{\\text{FV0}} \\rangle \\text{ (ns)}$", nBinsTotAmpl / 10, totAmplMin, totAmplMax, nBinsTRes, tResMin, tResMax, conditions.size(), 0, conditions.size()));
    objs2D[&ooFT0TimeAFV0TimeDiffVsFV0TotAmpl] = [&]() { return std::make_pair(objs[&ooFV0TotAmpl](), objs[&ooFT0TimeAFV0TimeDiff]()); };

    ooFT0TotAmplVsFV0TotAmpl.setObject(new TH3F(ooFT0TotAmplVsFV0TotAmpl.label.c_str(), "FT0 total amplitude vs FV0 total amplitude;FV0 total amplitude (ADC);FT0 total amplitude (ADC)", nBinsTotAmpl / 10, totAmplMin, totAmplMax, nBinsTotAmpl / 10, totAmplMin, totAmplMax, conditions.size(), 0, conditions.size()));
    objs2D[&ooFT0TotAmplVsFV0TotAmpl] = [&]() { return std::make_pair(objs[&ooFV0TotAmpl](), objs[&ooFT0TotAmpl]()); };

    ooFT0TotAmplAVsFV0TotAmpl.setObject(new TH3F(ooFT0TotAmplAVsFV0TotAmpl.label.c_str(), "FT0A total amplitude vs FV0 total amplitude;FV0 total amplitude (ADC);FT0A total amplitude (ADC)", nBinsTotAmpl / 10, totAmplMin, totAmplMax, nBinsTotAmpl / 10, totAmplMin, totAmplMax, conditions.size(), 0, conditions.size()));
    objs2D[&ooFT0TotAmplAVsFV0TotAmpl] = [&]() { return std::make_pair(objs[&ooFV0TotAmpl](), objs[&ooFT0TotAmplA]()); };

    ooFT0TotAmplCVsFV0TotAmpl.setObject(new TH3F(ooFT0TotAmplCVsFV0TotAmpl.label.c_str(), "FT0C total amplitude vs FV0 total amplitude;FV0 total amplitude (ADC);FT0C total amplitude (ADC)", nBinsTotAmpl / 10, totAmplMin, totAmplMax, nBinsTotAmpl / 10, totAmplMin, totAmplMax, conditions.size(), 0, conditions.size()));
    objs2D[&ooFT0TotAmplCVsFV0TotAmpl] = [&]() { return std::make_pair(objs[&ooFV0TotAmpl](), objs[&ooFT0TotAmplC]()); };

    ooFT0TotAmplVsFDDTotAmpl.setObject(new TH3F(ooFT0TotAmplVsFDDTotAmpl.label.c_str(), "FT0 total amplitude vs FDD total amplitude;FDD total amplitude (ADC);FT0 total amplitude (ADC)", nBinsTotAmpl / 10, totAmplMin, totAmplMax, nBinsTotAmpl / 10, totAmplMin, totAmplMax, conditions.size(), 0, conditions.size()));
    objs2D[&ooFT0TotAmplVsFDDTotAmpl] = [&]() { return std::make_pair(objs[&ooFDDTotAmpl](), objs[&ooFT0TotAmpl]()); };

    // ooFT0TotAmplVsFDDTotAmplA.setObject(new TH3F(ooFT0TotAmplVsFDDTotAmplA.label.c_str(), "FT0 total amplitude vs FDD-A total amplitude;FDD-A total amplitude (ADC);FT0 total amplitude (ADC)", nBinsTotAmpl / 10, totAmplMin, totAmplMax, nBinsTotAmpl / 10, totAmplMin, totAmplMax, conditions.size(), 0, conditions.size()));
    // objs2D[&ooFT0TotAmplVsFDDTotAmplA] = [&]() { return std::make_pair(objs[&ooFDDTotAmplA](), objs[&ooFT0TotAmpl]()); };

    // ooFT0TotAmplVsFDDTotAmplC.setObject(new TH3F(ooFT0TotAmplVsFDDTotAmplC.label.c_str(), "FT0 total amplitude vs FDD-C total amplitude;FDD-C total amplitude (ADC);FT0 total amplitude (ADC)", nBinsTotAmpl / 10, totAmplMin, totAmplMax, nBinsTotAmpl / 10, totAmplMin, totAmplMax, conditions.size(), 0, conditions.size()));
    // objs2D[&ooFT0TotAmplVsFDDTotAmplC] = [&]() { return std::make_pair(objs[&ooFDDTotAmplC](), objs[&ooFT0TotAmpl]()); };

    // ooFT0TotAmplAVsFDDTotAmpl.setObject(new TH3F(ooFT0TotAmplAVsFDDTotAmpl.label.c_str(), "FT0A total amplitude vs FDD total amplitude;FDD total amplitude (ADC);FT0A total amplitude (ADC)", nBinsTotAmpl / 10, totAmplMin, totAmplMax, nBinsTotAmpl / 10, totAmplMin, totAmplMax, conditions.size(), 0, conditions.size()));
    // objs2D[&ooFT0TotAmplAVsFDDTotAmpl] = [&]() { return std::make_pair(objs[&ooFDDTotAmpl](), objs[&ooFT0TotAmplA]()); };

    // ooFT0TotAmplAVsFDDTotAmplA.setObject(new TH3F(ooFT0TotAmplAVsFDDTotAmplA.label.c_str(), "FT0A total amplitude vs FDD-A total amplitude;FDD-A total amplitude (ADC);FT0A total amplitude (ADC)", nBinsTotAmpl / 10, totAmplMin, totAmplMax, nBinsTotAmpl / 10, totAmplMin, totAmplMax, conditions.size(), 0, conditions.size()));
    // objs2D[&ooFT0TotAmplAVsFDDTotAmplA] = [&]() { return std::make_pair(objs[&ooFDDTotAmplA](), objs[&ooFT0TotAmplA]()); };

    // ooFT0TotAmplAVsFDDTotAmplC.setObject(new TH3F(ooFT0TotAmplAVsFDDTotAmplC.label.c_str(), "FT0A total amplitude vs FDD-C total amplitude;FDD-C total amplitude (ADC);FT0A total amplitude (ADC)", nBinsTotAmpl / 10, totAmplMin, totAmplMax, nBinsTotAmpl / 10, totAmplMin, totAmplMax, conditions.size(), 0, conditions.size()));
    // objs2D[&ooFT0TotAmplAVsFDDTotAmplC] = [&]() { return std::make_pair(objs[&ooFDDTotAmplC](), objs[&ooFT0TotAmplA]()); };

    // ooFT0TotAmplCVsFDDTotAmpl.setObject(new TH3F(ooFT0TotAmplCVsFDDTotAmpl.label.c_str(), "FT0C total amplitude vs FDD total amplitude;FDD total amplitude (ADC);FT0C total amplitude (ADC)", nBinsTotAmpl / 10, totAmplMin, totAmplMax, nBinsTotAmpl / 10, totAmplMin, totAmplMax, conditions.size(), 0, conditions.size()));
    // objs2D[&ooFT0TotAmplCVsFDDTotAmpl] = [&]() { return std::make_pair(objs[&ooFDDTotAmpl](), objs[&ooFT0TotAmplC]()); };

    // ooFT0TotAmplCVsFDDTotAmplA.setObject(new TH3F(ooFT0TotAmplCVsFDDTotAmplA.label.c_str(), "FT0C total amplitude vs FDD-A total amplitude;FDD-A total amplitude (ADC);FT0C total amplitude (ADC)", nBinsTotAmpl / 10, totAmplMin, totAmplMax, nBinsTotAmpl / 10, totAmplMin, totAmplMax, conditions.size(), 0, conditions.size()));
    // objs2D[&ooFT0TotAmplCVsFDDTotAmplA] = [&]() { return std::make_pair(objs[&ooFDDTotAmplA](), objs[&ooFT0TotAmplC]()); };

    // ooFT0TotAmplCVsFDDTotAmplC.setObject(new TH3F(ooFT0TotAmplCVsFDDTotAmplC.label.c_str(), "FT0C total amplitude vs FDD-C total amplitude;FDD-C total amplitude (ADC);FT0C total amplitude (ADC)", nBinsTotAmpl / 10, totAmplMin, totAmplMax, nBinsTotAmpl / 10, totAmplMin, totAmplMax, conditions.size(), 0, conditions.size()));
    // objs2D[&ooFT0TotAmplCVsFDDTotAmplC] = [&]() { return std::make_pair(objs[&ooFDDTotAmplC](), objs[&ooFT0TotAmplC]()); };

    ooFV0TotAmplVsFDDTotAmpl.setObject(new TH3F(ooFV0TotAmplVsFDDTotAmpl.label.c_str(), "FV0 total amplitude vs FDD total amplitude;FDD total amplitude (ADC);FV0 total amplitude (ADC)", nBinsTotAmpl / 10, totAmplMin, totAmplMax, nBinsTotAmpl / 10, totAmplMin, totAmplMax, conditions.size(), 0, conditions.size()));
    objs2D[&ooFV0TotAmplVsFDDTotAmpl] = [&]() { return std::make_pair(objs[&ooFDDTotAmpl](), objs[&ooFV0TotAmpl]()); };

    // ooFV0TotAmplVsFDDTotAmplA.setObject(new TH3F(ooFV0TotAmplVsFDDTotAmplA.label.c_str(), "FV0 total amplitude vs FDD-A total amplitude;FDD-A total amplitude (ADC);FV0 total amplitude (ADC)", nBinsTotAmpl / 10, totAmplMin, totAmplMax, nBinsTotAmpl / 10, totAmplMin, totAmplMax, conditions.size(), 0, conditions.size()));
    // objs2D[&ooFV0TotAmplVsFDDTotAmplA] = [&]() { return std::make_pair(objs[&ooFDDTotAmplA](), objs[&ooFV0TotAmpl]()); };

    // ooFV0TotAmplVsFDDTotAmplC.setObject(new TH3F(ooFV0TotAmplVsFDDTotAmplC.label.c_str(), "FV0 total amplitude vs FDD-C total amplitude;FDD-C total amplitude (ADC);FV0 total amplitude (ADC)", nBinsTotAmpl / 10, totAmplMin, totAmplMax, nBinsTotAmpl / 10, totAmplMin, totAmplMax, conditions.size(), 0, conditions.size()));
    // objs2D[&ooFV0TotAmplVsFDDTotAmplC] = [&]() { return std::make_pair(objs[&ooFDDTotAmplC](), objs[&ooFV0TotAmpl]()); };

    // 2D quantities per X bin

    ooFT0AmplPerCh.setObject(new TH3F(ooFT0AmplPerCh.label.c_str(), "FT0 channel amplitudes;Channel ID;Amplitude (ADC)", nChFT0, 0, static_cast<float>(nChFT0), nADC, 0.f, static_cast<float>(nADC), conditions.size(), 0, conditions.size()));
    objs2DPerXBin[&ooFT0AmplPerCh] = [&](int chID) { return ft0ChAmpl[chID]; };

    ooFV0AmplPerCh.setObject(new TH3F(ooFV0AmplPerCh.label.c_str(), "FV0 channel amplitudes;Channel ID;Amplitude (ADC)", nChFV0, 0, static_cast<float>(nChFV0), nADC, 0.f, static_cast<float>(nADC), conditions.size(), 0, conditions.size()));
    objs2DPerXBin[&ooFV0AmplPerCh] = [&](int chID) { return fv0ChAmpl[chID]; };

    ooFDDAmplPerCh.setObject(new TH3F(ooFDDAmplPerCh.label.c_str(), "FDD channel amplitudes;Channel ID;Amplitude (ADC)", nChFDD, 0, static_cast<float>(nChFDD), nADC, 0.f, static_cast<float>(nADC), conditions.size(), 0, conditions.size()));
    objs2DPerXBin[&ooFDDAmplPerCh] = [&](int chID) { return fddChAmpl[chID]; };

    // Set event selection condition names on as bin labels

    for (auto &h : objs) {
      for (size_t c = 0; c < conditions.size(); c++) {
        h.first->object->GetYaxis()->SetBinLabel(c + 1, conditions[c].title.c_str());
      }
    }
    for (auto &h : objsft0) {
      for (size_t c = 0; c < ft0Conditions.size(); c++) {
        h.first->object->GetYaxis()->SetBinLabel(c + 1, ft0Conditions[c].title.c_str());
      }
    }
    for (auto &h : objsfv0) {
      for (size_t c = 0; c < fv0Conditions.size(); c++) {
        h.first->object->GetYaxis()->SetBinLabel(c + 1, fv0Conditions[c].title.c_str());
      }
    }
    for (auto &h : objsfdd) {
      for (size_t c = 0; c < fddConditions.size(); c++) {
        h.first->object->GetYaxis()->SetBinLabel(c + 1, fddConditions[c].title.c_str());
      }
    }
    for (auto &h : objs2D) {
      for (size_t c = 0; c < conditions.size(); c++) {
        h.first->object->GetZaxis()->SetBinLabel(c + 1, conditions[c].title.c_str());
      }
    }
    for (auto &h : objs2DPerXBin) {
      for (size_t c = 0; c < conditions.size(); c++) {
        h.first->object->GetZaxis()->SetBinLabel(c + 1, conditions[c].title.c_str());
      }
    }
  }

  // TODO: add configurable to enable processing of Extras
  // NOTE: this processing of FIT data is only for collisions matched to a BC
  //       -> bias as we don't process all FIT data
  void process(soa::Join<aod::Collisions, aod::EvSels, aod::FT0sCorrected>::iterator const& collision,
               aod::FT0s const&, aod::FV0As const&, aod::FDDs const&) {
    resetVariables();
    
    // TODO: validity checks to force own default values?

    isSel8 = collision.sel8();
    hasFT0 = collision.has_foundFT0();
    hasFV0 = collision.has_foundFV0();
    hasFDD = collision.has_foundFDD();

    pv = collision.posZ();
    nContrib = collision.numContrib();
    ft0timeACorr = collision.t0ACorrected();
    ft0timeCCorr = collision.t0ACorrected();
    ft0timeRes = collision.t0resolution();

    if (hasFT0) {
      auto ft0 = collision.foundFT0();
      ft0Triggers = ft0.triggerMask();
      // isFT0VTX = ft0Triggers[o2::fit::Triggers::bitVertex];
      // isFT0CE = ft0Triggers[o2::fit::Triggers::bitCen];
      // isFT0SCE = ft0Triggers[o2::fit::Triggers::bitSCen];

      ft0timeA = ft0.timeA();
      ft0timeC = ft0.timeC();
      ft0time = ft0.collTime();
      ft0vtx = ft0.posZ();

      for (size_t i = 0; i < ft0.amplitudeA().size(); i++) {
        ft0ChAmpl[ft0.channelA()[i]] = ft0.amplitudeA()[i];
        ft0TotAmpl += ft0.amplitudeA()[i];
        ft0TotAmplA += ft0.amplitudeA()[i];
      }
      for (size_t i = 0; i < ft0.amplitudeC().size(); i++) {
        ft0ChAmpl[ft0.channelC()[i] + nChFT0A] = ft0.amplitudeC()[i]; // Channel IDs in the C-side array start from zero in AO2D (JIRA AFIT-129)
        ft0TotAmpl += ft0.amplitudeC()[i];
        ft0TotAmplC += ft0.amplitudeC()[i];
      }
      
    }

    if (hasFV0) {
      auto fv0 = collision.foundFV0();
      fv0Triggers = fv0.triggerMask();
      // isFV0ORA = fv0Triggers[o2::fit::Triggers::bitA];
      // isFV0CH = fv0Triggers[o2::fit::Triggers::bitTrgNchan];
      // isFV0IN = fv0Triggers[o2::fit::Triggers::bitAIn];

      fv0time = fv0.time();

      for (size_t i = 0; i < fv0.amplitude().size(); i++) {
        fv0ChAmpl[fv0.channel()[i]] = fv0.amplitude()[i];
        fv0TotAmpl += fv0.amplitude()[i];
      }
    }

    if (hasFDD) {
      auto fdd = collision.foundFDD();
      fddTriggers = fdd.triggerMask();
      // isFDDORA = fddTriggers[o2::fit::Triggers::bitA];
      // isFDDORC = fddTriggers[o2::fit::Triggers::bitC];
      // isFDDVTX = fddTriggers[o2::fit::Triggers::bitVertex];

      fddtimeA = fdd.timeA();
      fddtimeC = fdd.timeC();

      // TODO: don't hard code?
      for (size_t i = 0; i < 8; i++) {
        fddChAmpl[i + 8] = fdd.chargeA()[i];
        fddTotAmpl += fdd.chargeA()[i];
        fddTotAmplA += fdd.chargeA()[i];
      }
      for (size_t i = 0; i < 8; i++) {
        fddChAmpl[i] = fdd.chargeC()[i];
        fddTotAmpl += fdd.chargeC()[i];
        fddTotAmplC += fdd.chargeC()[i];
      }
    }

    // isFT0VTXandFV0ORA = isFT0VTX && isFV0ORA;
    // isFT0VTXandFDDVTX = isFT0VTX && isFDDVTX;

    /* Evaluate conditions and fill stats histo */

    std::vector<bool> condResults(conditions.size(), false);
    std::vector<bool> ft0CondResults(ft0Conditions.size(), false);
    std::vector<bool> fv0CondResults(fv0Conditions.size(), false);
    std::vector<bool> fddCondResults(fddConditions.size(), false);

    for (size_t c = 0; c < conditions.size(); c++) {
      condResults[c] = conditions[c].eval();
      if (condResults[c]) {
        ooStats->Fill(c + 0.5);
      }
    }

    for (size_t c = 0; c < ft0Conditions.size(); c++) {
      ft0CondResults[c] = ft0Conditions[c].eval();
    }
    for (size_t c = 0; c < fv0Conditions.size(); c++) {
      fv0CondResults[c] = fv0Conditions[c].eval();
    }
    for (size_t c = 0; c < fddConditions.size(); c++) {
      fddCondResults[c] = fddConditions[c].eval();
    }

    /* Fill histograms */

    for (auto &h : objs) {
      float v = h.second();
      for (size_t c = 0; c < conditions.size(); c++) {
        if (condResults[c]) {
          (*h.first)->Fill(v, c + 0.5);
        }
      }
    }

    for (auto &h : objsft0) {
      float v = h.second();
      for (size_t c = 0; c < ft0Conditions.size(); c++) {
        if (ft0CondResults[c]) {
          (*h.first)->Fill(v, c + 0.5f);
        }
      }
    }

    for (auto &h : objsfv0) {
      float v = h.second();
      for (size_t c = 0; c < fv0Conditions.size(); c++) {
        if (fv0CondResults[c]) {
          (*h.first)->Fill(v, c + 0.5f);
        }
      }
    }

    for (auto &h : objsfdd) {
      float v = h.second();
      for (size_t c = 0; c < fddConditions.size(); c++) {
        if (fddCondResults[c]) {
          (*h.first)->Fill(v, c + 0.5f);
        }
      }
    }

    for (auto &h : objs2D) {
      auto v = h.second();
      for (size_t c = 0; c < conditions.size(); c++) {
        if (condResults[c]) {
          (*h.first)->Fill(v.first, v.second, c + 0.5f);
        }
      }
    }

    for (auto &h : objs2DPerXBin) {
      for (int xBin = 0; xBin < h.first->object->GetXaxis()->GetNbins(); xBin++) {
        float vy = h.second(xBin);
        for (size_t c = 0; c < conditions.size(); c++) {
          if (condResults[c] && vy >= 0.f) { // TODO: allow negative, set some default value in lambda instead
            (*h.first)->Fill(xBin + 0.5f, vy, c + 0.5f);
          }
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<fitQa>(cfgc, TaskName{"fit-qa"})};
}