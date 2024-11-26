#include <cstddef>
#include <tuple>
#include <vector>

#include "component.h"
#include "inputreader.h"
#include "mixture_prediction.h"

#ifdef PYBUILD
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
namespace py = pybind11;
#endif  // PYBUILD

/**
 * \brief Simulates a breakthrough process in an adsorption column.
 *
 * The Breakthrough struct encapsulates the parameters and methods required to simulate
 * the breakthrough of gases in an adsorption column. It handles the initialization of
 * simulation parameters, computation of time steps, and generation of output scripts
 * for plotting and visualization.
 */
struct Breakthrough
{
 public:
  /**
   * \brief Constructs a Breakthrough simulation using an InputReader.
   *
   * Initializes a Breakthrough object with parameters specified in the InputReader.
   *
   * \param inputreader Reference to an InputReader containing simulation parameters.
   */
  Breakthrough(const InputReader &inputreader);

  /**
   * \brief Constructs a Breakthrough simulation with specified parameters.
   *
   * Initializes a Breakthrough object with the provided simulation parameters.
   *
   * \param _displayName Name of the simulation for display purposes.
   * \param _components Vector of components involved in the simulation.
   * \param _carrierGasComponent Index of the carrier gas component.
   * \param _numberOfGridPoints Number of grid points in the column.
   * \param _printEvery Frequency of printing time steps to the screen.
   * \param _writeEvery Frequency of writing data to files.
   * \param _temperature Simulation temperature in Kelvin.
   * \param _p_total Total pressure in the column in Pascals.
   * \param _columnVoidFraction Void fraction of the column.
   * \param _pressureGradient Pressure gradient in the column.
   * \param _particleDensity Particle density in kg/m³.
   * \param _columnEntranceVelocity Interstitial velocity at the beginning of the column in m/s.
   * \param _columnLength Length of the column in meters.
   * \param _timeStep Time step for the simulation.
   * \param _numberOfTimeSteps Total number of time steps.
   * \param _autoSteps Flag to use automatic number of steps.
   * \param _pulse Flag to indicate pulsed inlet condition.
   * \param _pulseTime Pulse time.
   * \param _mixture MixturePrediction object for mixture predictions.
   */
  Breakthrough(std::string _displayName, std::vector<Component> _components, size_t _carrierGasComponent,
               size_t _numberOfGridPoints, size_t _printEvery, size_t _writeEvery, double _temperature, double _p_total,
               double _columnVoidFraction, double _pressureGradient, double _particleDensity,
               double _columnEntranceVelocity, double _columnLength, double _timeStep, size_t _numberOfTimeSteps,
               bool _autoSteps, bool _pulse, double _pulseTime, const MixturePrediction _mixture);

  /**
   * \brief Prints the representation of the Breakthrough object to the console.
   */
  void print() const;

  /**
   * \brief Returns a string representation of the Breakthrough object.
   *
   * \return A string representing the Breakthrough object.
   */
  std::string repr() const;

  /**
   * \brief Initializes the Breakthrough simulation.
   *
   * Sets up initial conditions and precomputes factors required for the simulation.
   */
  void initialize();

  /**
   * \brief Runs the Breakthrough simulation.
   *
   * Executes the simulation over the specified number of time steps.
   */
  void run();

  /**
   * \brief Creates a Gnuplot script for plotting breakthrough curves.
   *
   * Generates a Gnuplot script to visualize the simulation results.
   */
  void createPlotScript();

  /**
   * \brief Creates scripts for generating movies of the simulation.
   *
   * Generates scripts to create movies visualizing the simulation over time.
   */
  void createMovieScripts();

#ifdef PYBUILD
  /**
   * \brief Computes the Breakthrough simulation and returns the results.
   *
   * Executes the simulation and returns a NumPy array containing the simulation data.
   *
   * \return A NumPy array of simulation results.
   */
  py::array_t<double> compute();

  /**
   * \brief Sets the component parameters for the simulation.
   *
   * Updates the mole fractions and isotherm parameters for each component.
   *
   * \param molfracs Vector of mole fractions for each component.
   * \param params Vector of isotherm parameters for the components.
   */
  void setComponentsParameters(std::vector<double> molfracs, std::vector<double> params);

  /**
   * \brief Retrieves the component parameters used in the simulation.
   *
   * Returns the isotherm parameters for each component.
   *
   * \return A vector containing the isotherm parameters for the components.
   */
  std::vector<double> getComponentsParameters();
#endif  // PYBUILD

 private:
  const std::string displayName;      ///< Name of the simulation for display purposes.
  std::vector<Component> components;  ///< Vector of components involved in the simulation.
  size_t carrierGasComponent{0};      ///< Index of the carrier gas component.
  size_t Ncomp;                       ///< Number of components.
  size_t Ngrid;                       ///< Number of grid points.

  size_t printEvery;  ///< Frequency of printing time steps to the screen.
  size_t writeEvery;  ///< Frequency of writing data to files.

  double T;        ///< Absolute temperature in Kelvin.
  double p_total;  ///< Total pressure column [Pa].
  double dptdx;    ///< Pressure gradient [N/m³].
  double epsilon;  ///< Void-fraction of the column [-].
  double rho_p;    ///< Particle density [kg/m³].
  double v_in;     ///< Interstitial velocity at the beginning of the column [m/s].

  double L;                                         ///< Length of the column.
  double dx;                                        ///< Spacing in spatial direction.
  double dt;                                        ///< Time step for integration.
  size_t Nsteps;                                    ///< Total number of steps.
  bool autoSteps;                                   ///< Flag to use automatic number of steps.
  bool pulse;                                       ///< Pulsed inlet condition for breakthrough.
  double tpulse;                                    ///< Pulse time.
  MixturePrediction mixture;                        ///< MixturePrediction object for mixture predictions.
  size_t maxIsothermTerms;                          ///< Maximum number of isotherm terms.
  std::pair<size_t, size_t> iastPerformance{0, 0};  ///< Performance metrics for IAST calculations.

  // vector of size 'Ncomp'
  std::vector<double> prefactor;  ///< Precomputed factors for mass transfer.
  std::vector<double> Yi;         ///< Ideal gas mole fractions for each component.
  std::vector<double> Xi;         ///< Adsorbed mole fractions for each component.
  std::vector<double> Ni;         ///< Number of molecules for each component.

  // vector of size '(Ngrid + 1)'
  std::vector<double> V;     ///< Interstitial gas velocity along the column.
  std::vector<double> Vnew;  ///< Updated interstitial gas velocities.
  std::vector<double> Pt;    ///< Total pressure along the column.

  // vector of size '(Ngrid + 1) * Ncomp', for each grid point, data per component (contiguous)
  std::vector<double> P;          ///< Partial pressure at every grid point for each component.
  std::vector<double> Pnew;       ///< Updated partial pressures.
  std::vector<double> Q;          ///< Volume-averaged adsorption amount at every grid point for each component.
  std::vector<double> Qnew;       ///< Updated adsorption amounts.
  std::vector<double> Qeq;        ///< Equilibrium adsorption amount at every grid point for each component.
  std::vector<double> Qeqnew;     ///< Updated equilibrium adsorption amounts.
  std::vector<double> Dpdt;       ///< Derivative of P with respect to time.
  std::vector<double> Dpdtnew;    ///< Updated derivative of P with respect to time.
  std::vector<double> Dqdt;       ///< Derivative of Q with respect to time.
  std::vector<double> Dqdtnew;    ///< Updated derivative of Q with respect to time.
  std::vector<double> cachedP0;   ///< Cached hypothetical pressure.
  std::vector<double> cachedPsi;  ///< Cached reduced grand potential over the column.

  enum class IntegrationScheme
  {
    SSP_RK = 0,    ///< Strong Stability Preserving Runge-Kutta method.
    Iterative = 1  ///< Iterative integration scheme.
  };

  /**
   * \brief Computes the first derivatives of concentrations and pressures.
   *
   * Calculates the derivatives Dq/dt and Dp/dt along the column.
   *
   * \param dqdt Output vector for the derivatives of Q with respect to time.
   * \param dpdt Output vector for the derivatives of P with respect to time.
   * \param q_eq Equilibrium adsorption amounts.
   * \param q Current adsorption amounts.
   * \param v Interstitial gas velocities.
   * \param p Partial pressures.
   */
  void computeFirstDerivatives(std::vector<double> &dqdt, std::vector<double> &dpdt, const std::vector<double> &q_eq,
                               const std::vector<double> &q, const std::vector<double> &v,
                               const std::vector<double> &p);

  /**
   * \brief Computes a single simulation step.
   *
   * Advances the simulation by one time step.
   *
   * \param step The current time step index.
   */
  void computeStep(size_t step);

  /**
   * \brief Computes the equilibrium loadings for the current time step.
   *
   * Calculates the equilibrium adsorption amounts based on current pressures.
   */
  void computeEquilibriumLoadings();

  /**
   * \brief Computes the interstitial gas velocities along the column.
   *
   * Updates the velocities based on current pressures and adsorption amounts.
   */
  void computeVelocity();

  /**
   * \brief Creates a script to generate a movie for the interstitial gas velocity.
   */
  void createMovieScriptColumnV();

  /**
   * \brief Creates a script to generate a movie for the total pressure along the column.
   */
  void createMovieScriptColumnPt();

  /**
   * \brief Creates a script to generate a movie for the adsorption amounts Q.
   */
  void createMovieScriptColumnQ();

  /**
   * \brief Creates a script to generate a movie for the equilibrium adsorption amounts Qeq.
   */
  void createMovieScriptColumnQeq();

  /**
   * \brief Creates a script to generate a movie for the partial pressures P.
   */
  void createMovieScriptColumnP();

  /**
   * \brief Creates a script to generate a movie for the derivatives of pressure Dp/dt.
   */
  void createMovieScriptColumnDpdt();

  /**
   * \brief Creates a script to generate a movie for the derivatives of adsorption amounts Dq/dt.
   */
  void createMovieScriptColumnDqdt();

  /**
   * \brief Creates a script to generate a movie for the normalized partial pressures.
   */
  void createMovieScriptColumnPnormalized();
};
