#pragma once

#include <vector>
#include <string>

#include "component.h"

extern bool startsWith(const std::string &str, const std::string &prefix);
extern std::string trim(const std::string& s);

struct InputReader
{
  InputReader(const std::string fileName);

  enum class SimulationType
  {
    Breakthrough = 0,
    MixturePrediction = 1,
    Fitting = 2,
    Test = 3
  };

  std::vector<Component> components;
  size_t numberOfCarrierGases{ 0 };
  size_t carrierGasComponent{ 0 };
  size_t maxIsothermTerms{ 0 };

  SimulationType simulationType{ SimulationType::Breakthrough };
  size_t mixturePredictionMethod{ 0 };
  size_t IASTMethod{ 0 };
  std::string displayName{"Column"};
  double temperature{ 433.0 };
  double columnVoidFraction{ 0.4 };
  double particleDensity{ 1000.0 };
  double totalPressure{ 1.0e6 };
  double pressureGradient{ 0.0 };
  double columnEntranceVelocity{ 0.1 };
  double columnLength { 0.3 };

  size_t numberOfTimeSteps{ 0 };
  bool autoNumberOfTimeSteps{ true };
  double timeStep{ 0.0005 };
  bool pulseBreakthrough{ false };
  double pulseTime{ 0.0 };
  size_t printEvery{ 10000 };
  size_t writeEvery{ 10000 };
  size_t numberOfGridPoints{ 100 };

  double pressureStart{ -1.0 };
  double pressureEnd{ -1.0 };
  size_t numberOfPressurePoints{ 100 };
  size_t pressureScale{ 0 };

  size_t columnPressure{ 0 };
  size_t columnLoading{ 1 };
  size_t columnError{ 2 };
};
