#pragma once

#include <string>
#include <vector>

#include "component.h"

extern bool startsWith(const std::string &str, const std::string &prefix);
extern std::string trim(const std::string &s);

/**
 * \brief Parses input files and stores simulation parameters.
 *
 * The InputReader struct is responsible for reading and parsing input files containing simulation parameters.
 * It stores all the necessary data required to set up and run simulations, including components, simulation types,
 * and various parameters related to the simulation environment.
 */
struct InputReader
{
  /**
   * \brief Constructs an InputReader and parses the given input file.
   *
   * \param fileName The name of the input file to parse.
   */
  InputReader(const std::string fileName);

  /**
   * \brief Enumerates the types of simulations supported.
   */
  enum class SimulationType
  {
    Breakthrough = 0,       ///< Breakthrough simulation.
    MixturePrediction = 1,  ///< Mixture prediction simulation.
    Fitting = 2,            ///< Fitting simulation.
    Test = 3                ///< Test simulation.
  };

  std::vector<Component> components;  ///< The list of components involved in the simulation.
  size_t numberOfCarrierGases{0};     ///< The number of carrier gas components.
  size_t carrierGasComponent{0};      ///< The index of the carrier gas component.
  size_t maxIsothermTerms{0};         ///< The maximum number of isotherm terms among all components.

  SimulationType simulationType{SimulationType::Breakthrough};  ///< The type of simulation to perform.
  size_t mixturePredictionMethod{0};                            ///< The method used for mixture prediction.
  size_t IASTMethod{0};                                         ///< The method used for IAST calculations.
  std::string displayName{"Column"};                            ///< The display name for the simulation.
  double temperature{433.0};                                    ///< The simulation temperature in Kelvin.
  double columnVoidFraction{0.4};                               ///< The void fraction of the column.
  double particleDensity{1000.0};                               ///< The density of the particles in kg/m^3.
  double totalPressure{1.0e6};                                  ///< The total pressure in the system in Pa.
  double pressureGradient{0.0};                                 ///< The pressure gradient in the column.
  double columnEntranceVelocity{0.1};                           ///< The entrance velocity of the column in m/s.
  double columnLength{0.3};                                     ///< The length of the column in meters.

  size_t numberOfTimeSteps{0};       ///< The number of time steps in the simulation.
  bool autoNumberOfTimeSteps{true};  ///< Whether to automatically determine the number of time steps.
  double timeStep{0.0005};           ///< The time step size in seconds.
  bool pulseBreakthrough{false};     ///< Whether to use pulse breakthrough mode.
  double pulseTime{0.0};             ///< The duration of the pulse in seconds.
  size_t printEvery{10000};          ///< The interval at which to print output.
  size_t writeEvery{10000};          ///< The interval at which to write output.
  size_t numberOfGridPoints{100};    ///< The number of grid points in the column.

  double pressureStart{-1.0};          ///< The starting pressure for isotherm calculations.
  double pressureEnd{-1.0};            ///< The ending pressure for isotherm calculations.
  size_t numberOfPressurePoints{100};  ///< The number of pressure points to calculate.
  size_t pressureScale{0};             ///< The scale for pressure calculations (0 for log, 1 for linear).

  size_t columnPressure{0};  ///< The index of the column for pressure data.
  size_t columnLoading{1};   ///< The index of the column for loading data.
  size_t columnError{2};     ///< The index of the column for error data.
};
