#pragma once

#include <tuple>
#include <vector>

#ifdef PYBUILD
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
namespace py = pybind11;
#endif  // PYBUILD

#include "component.h"
#include "inputreader.h"

/**
 * \brief Class for predicting mixture adsorption isotherms.
 *
 * The MixturePrediction class provides methods to predict mixture adsorption isotherms using various methods like IAST,
 * SIAST, and explicit isotherm models. It handles the computation of adsorbed phase mole fractions and loadings based
 * on gas phase compositions and pressure.
 */
class MixturePrediction
{
 public:
  /**
   * \brief Enum class for prediction methods.
   *
   * Specifies the method used for predicting mixture adsorption isotherms.
   */
  enum class PredictionMethod
  {
    IAST = 0,   ///< Ideal Adsorbed Solution Theory
    SIAST = 1,  ///< Segregated Ideal Adsorbed Solution Theory
    EI = 2,     ///< Explicit Isotherm
    SEI = 3     ///< Segregated Explicit Isotherm
  };

  /**
   * \brief Enum class for IAST methods.
   *
   * Specifies the method used for solving IAST equations.
   */
  enum class IASTMethod
  {
    FastIAST = 0,            ///< Fast IAST algorithm
    NestedLoopBisection = 1  ///< Nested Loop Bisection method
  };

  /**
   * \brief Constructs a MixturePrediction object from an InputReader.
   *
   * Initializes the MixturePrediction instance using the parameters provided by the InputReader.
   *
   * \param inputreader The InputReader containing the simulation parameters.
   */
  MixturePrediction(const InputReader &inputreader);

  /**
   * \brief Constructs a MixturePrediction object with specified parameters.
   *
   * Initializes the MixturePrediction instance using the provided parameters.
   *
   * \param _displayName The display name for the simulation.
   * \param _components A vector of Component objects representing the mixture components.
   * \param _numberOfCarrierGases The number of carrier gases in the mixture.
   * \param _carrierGasComponent The index of the carrier gas component.
   * \param _temperature The temperature of the system.
   * \param _pressureStart The starting pressure for the simulation.
   * \param _pressureEnd The ending pressure for the simulation.
   * \param _numberOfPressurePoints The number of pressure points in the simulation.
   * \param _pressureScale The pressure scale (0 for Log, 1 for Normal).
   * \param _predictionMethod The prediction method to use.
   * \param _iastMethod The IAST method to use.
   */
  MixturePrediction(std::string _displayName, std::vector<Component> _components, size_t _numberOfCarrierGases,
                    size_t _carrierGasComponent, double _temperature, double _pressureStart, double _pressureEnd,
                    size_t _numberOfPressurePoints, size_t _pressureScale, size_t _predictionMethod,
                    size_t _iastMethod);

  /**
   * \brief Prints the mixture prediction data.
   *
   * Outputs the mixture prediction data to the standard output.
   */
  void print() const;

  /**
   * \brief Returns a string representation of the MixturePrediction object.
   *
   * \return A string representing the MixturePrediction object.
   */
  std::string repr() const;

  /**
   * \brief Runs the mixture prediction simulation.
   *
   * Performs the mixture prediction calculations and writes the results to output files.
   */
  void run();

  /**
   * \brief Creates a Gnuplot script for pure component isotherms.
   *
   * Generates a Gnuplot script to plot pure component isotherms.
   */
  void createPureComponentsPlotScript();

  /**
   * \brief Creates a Gnuplot script for mixture prediction.
   *
   * Generates a Gnuplot script to plot mixture prediction results.
   */
  void createMixturePlotScript();

  /**
   * \brief Creates a Gnuplot script for mixture adsorbed mol fractions.
   *
   * Generates a Gnuplot script to plot adsorbed mol fractions in the mixture.
   */
  void createMixtureAdsorbedMolFractionPlotScript();

  /**
   * \brief Creates plot scripts for the mixture prediction.
   *
   * Generates the necessary Gnuplot scripts to plot the simulation results.
   */
  void createPlotScript();

#ifdef PYBUILD
  /**
   * \brief Computes the mixture prediction.
   *
   * Performs the mixture prediction calculations and returns the results as a NumPy array.
   *
   * \return A NumPy array containing the mixture prediction results.
   */
  py::array_t<double> compute();

  /**
   * \brief Sets the pressure range for the simulation.
   *
   * Updates the starting and ending pressures for the simulation.
   *
   * \param _pressureStart The new starting pressure.
   * \param _pressureEnd The new ending pressure.
   */
  void setPressure(double _pressureStart, double _pressureEnd);

  /**
   * \brief Sets the components' parameters.
   *
   * Updates the molar fractions and isotherm parameters of the components.
   *
   * \param molfracs A vector of molar fractions for each component.
   * \param params A vector of isotherm parameters for the components.
   */
  void setComponentsParameters(std::vector<double> molfracs, std::vector<double> params);

  /**
   * \brief Gets the components' parameters.
   *
   * Retrieves the isotherm parameters of the components.
   *
   * \return A vector containing the isotherm parameters of the components.
   */
  std::vector<double> getComponentsParameters();
#endif  // PYBUILD

  /**
   * \brief Gets the maximum number of isotherm terms.
   *
   * \return The maximum number of isotherm terms.
   */
  size_t getMaxIsothermTerms() const { return maxIsothermTerms; }

  /**
   * \brief Predicts the mixture adsorption isotherm.
   *
   * Computes the adsorbed phase mole fractions and loadings based on the gas phase compositions and pressure.
   *
   * \param Yi The gas phase mole fractions.
   * \param P The total pressure.
   * \param Xi The adsorbed phase mole fractions (output).
   * \param Ni The number of adsorbed molecules of each component (output).
   * \param cachedP0 An array to cache intermediate pressure calculations.
   * \param cachedPsi An array to cache intermediate psi calculations.
   * \return A pair containing the number of IAST steps and a status code.
   */
  std::pair<size_t, size_t> predictMixture(const std::vector<double> &Yi, const double &P, std::vector<double> &Xi,
                                           std::vector<double> &Ni, double *cachedP0, double *cachedPsi);

 private:
  std::string displayName;                  ///< The display name for the simulation.
  std::vector<Component> components;        ///< The vector of components in the mixture.
  std::vector<Component> sortedComponents;  ///< Components sorted according to specific criteria.
  const size_t Ncomp;                       ///< The total number of components.
  const size_t Nsorted;                     ///< The number of sorted components.
  size_t numberOfCarrierGases;              ///< The number of carrier gases in the mixture.
  size_t carrierGasComponent;               ///< The index of the carrier gas component.
  PredictionMethod predictionMethod;        ///< The method used for predicting mixture adsorption isotherms.
  IASTMethod iastMethod;                    ///< The method used for solving IAST equations.
  size_t maxIsothermTerms;                  ///< The maximum number of isotherm terms.
  std::vector<std::vector<Component>>
      segregatedSortedComponents;  ///< Segregated and sorted components for SIAST/SEI methods.

  std::vector<double> alpha1;      ///< Intermediate calculation vector for explicit isotherms.
  std::vector<double> alpha2;      ///< Intermediate calculation vector for explicit isotherms.
  std::vector<double> alpha_prod;  ///< Product of alpha values for explicit isotherms.
  std::vector<double> x;           ///< Intermediate calculation vector.

  std::vector<double> pstar;  ///< Hypothetical pressures in IAST calculations.
  std::vector<double> psi;    ///< Reduced grand potential values.
  std::vector<double> G;      ///< Intermediate calculation vector in IAST.
  std::vector<double> delta;  ///< Correction vector in IAST.
  std::vector<double> Phi;    ///< Jacobian matrix in IAST calculations.

  /**
   * \brief Enum class for pressure scales.
   *
   * Specifies the scale to use for pressure in the simulation.
   */
  enum class PressureScale
  {
    Log = 0,    ///< Logarithmic pressure scale
    Normal = 1  ///< Linear pressure scale
  };
  double temperature{300.0};                        ///< The temperature of the system.
  double pressureStart{1e3};                        ///< The starting pressure for the simulation.
  double pressureEnd{1e8};                          ///< The ending pressure for the simulation.
  size_t numberOfPressurePoints{100};               ///< The number of pressure points in the simulation.
  PressureScale pressureScale{PressureScale::Log};  ///< The pressure scale to use.

  /**
   * \brief Initializes the pressure points for the simulation.
   *
   * Generates a vector of pressure points based on the starting and ending pressures and the pressure scale.
   *
   * \return A vector containing the pressure points for the simulation.
   */
  std::vector<double> initPressures();

  /**
   * \brief Sorts the components based on specific criteria.
   *
   * Sorts the components to optimize calculations in prediction methods.
   */
  void sortComponents();

  /**
   * \brief Computes mixture prediction using Fast IAST method.
   *
   * \param Yi The gas phase mole fractions.
   * \param P The total pressure.
   * \param Xi The adsorbed phase mole fractions (output).
   * \param Ni The number of adsorbed molecules of each component (output).
   * \param cachedP0 An array to cache intermediate pressure calculations.
   * \param cachedPsi An array to cache intermediate psi calculations.
   * \return A pair containing the number of IAST steps and a status code.
   */
  std::pair<size_t, size_t> computeFastIAST(const std::vector<double> &Yi, const double &P, std::vector<double> &Xi,
                                            std::vector<double> &Ni, double *cachedP0, double *cachedPsi);

  /**
   * \brief Computes mixture prediction using Fast SIAST method.
   *
   * \param Yi The gas phase mole fractions.
   * \param P The total pressure.
   * \param Xi The adsorbed phase mole fractions (output).
   * \param Ni The number of adsorbed molecules of each component (output).
   * \param cachedP0 An array to cache intermediate pressure calculations.
   * \param cachedPsi An array to cache intermediate psi calculations.
   * \return A pair containing the number of IAST steps and a status code.
   */
  std::pair<size_t, size_t> computeFastSIAST(const std::vector<double> &Yi, const double &P, std::vector<double> &Xi,
                                             std::vector<double> &Ni, double *cachedP0, double *cachedPsi);

  /**
   * \brief Computes mixture prediction for a specific term using Fast SIAST method.
   *
   * \param term The index of the isotherm term.
   * \param Yi The gas phase mole fractions.
   * \param P The total pressure.
   * \param Xi The adsorbed phase mole fractions (output).
   * \param Ni The number of adsorbed molecules of each component (output).
   * \param cachedP0 An array to cache intermediate pressure calculations.
   * \param cachedPsi An array to cache intermediate psi calculations.
   * \return A pair containing the number of IAST steps and a status code.
   */
  std::pair<size_t, size_t> computeFastSIAST(size_t term, const std::vector<double> &Yi, const double &P,
                                             std::vector<double> &Xi, std::vector<double> &Ni, double *cachedP0,
                                             double *cachedPsi);

  /**
   * \brief Computes mixture prediction using IAST with nested loop bisection method.
   *
   * \param Yi The gas phase mole fractions.
   * \param P The total pressure.
   * \param Xi The adsorbed phase mole fractions (output).
   * \param Ni The number of adsorbed molecules of each component (output).
   * \param cachedP0 An array to cache intermediate pressure calculations.
   * \param cachedPsi An array to cache intermediate psi calculations.
   * \return A pair containing the number of IAST steps and a status code.
   */
  std::pair<size_t, size_t> computeIASTNestedLoopBisection(const std::vector<double> &Yi, const double &P,
                                                           std::vector<double> &Xi, std::vector<double> &Ni,
                                                           double *cachedP0, double *cachedPsi);

  /**
   * \brief Computes mixture prediction using SIAST with nested loop bisection method.
   *
   * \param Yi The gas phase mole fractions.
   * \param P The total pressure.
   * \param Xi The adsorbed phase mole fractions (output).
   * \param Ni The number of adsorbed molecules of each component (output).
   * \param cachedP0 An array to cache intermediate pressure calculations.
   * \param cachedPsi An array to cache intermediate psi calculations.
   * \return A pair containing the number of IAST steps and a status code.
   */
  std::pair<size_t, size_t> computeSIASTNestedLoopBisection(const std::vector<double> &Yi, const double &P,
                                                            std::vector<double> &Xi, std::vector<double> &Ni,
                                                            double *cachedP0, double *cachedPsi);

  /**
   * \brief Computes mixture prediction for a specific term using SIAST with nested loop bisection method.
   *
   * \param term The index of the isotherm term.
   * \param Yi The gas phase mole fractions.
   * \param P The total pressure.
   * \param Xi The adsorbed phase mole fractions (output).
   * \param Ni The number of adsorbed molecules of each component (output).
   * \param cachedP0 An array to cache intermediate pressure calculations.
   * \param cachedPsi An array to cache intermediate psi calculations.
   * \return A pair containing the number of IAST steps and a status code.
   */
  std::pair<size_t, size_t> computeSIASTNestedLoopBisection(size_t term, const std::vector<double> &Yi, const double &P,
                                                            std::vector<double> &Xi, std::vector<double> &Ni,
                                                            double *cachedP0, double *cachedPsi);

  /**
   * \brief Computes mixture prediction using explicit isotherm model.
   *
   * \param Yi The gas phase mole fractions.
   * \param P The total pressure.
   * \param Xi The adsorbed phase mole fractions (output).
   * \param Ni The number of adsorbed molecules of each component (output).
   * \return A pair containing the number of steps and a status code.
   */
  std::pair<size_t, size_t> computeExplicitIsotherm(const std::vector<double> &Yi, const double &P,
                                                    std::vector<double> &Xi, std::vector<double> &Ni);

  /**
   * \brief Computes mixture prediction using segregated explicit isotherm model.
   *
   * \param Yi The gas phase mole fractions.
   * \param P The total pressure.
   * \param Xi The adsorbed phase mole fractions (output).
   * \param Ni The number of adsorbed molecules of each component (output).
   * \return A pair containing the number of steps and a status code.
   */
  std::pair<size_t, size_t> computeSegratedExplicitIsotherm(const std::vector<double> &Yi, const double &P,
                                                            std::vector<double> &Xi, std::vector<double> &Ni);

  /**
   * \brief Computes mixture prediction for a specific term using segregated explicit isotherm model.
   *
   * \param site The index of the isotherm site.
   * \param Yi The gas phase mole fractions.
   * \param P The total pressure.
   * \param Xi The adsorbed phase mole fractions (output).
   * \param Ni The number of adsorbed molecules of each component (output).
   * \return A pair containing the number of steps and a status code.
   */
  std::pair<size_t, size_t> computeSegratedExplicitIsotherm(size_t site, const std::vector<double> &Yi, const double &P,
                                                            std::vector<double> &Xi, std::vector<double> &Ni);

  /**
   * \brief Prints error status for debugging purposes.
   *
   * Outputs the current state of variables when an error occurs in IAST calculations.
   *
   * \param psi The current psi value.
   * \param sum The current sum of mole fractions.
   * \param P The total pressure.
   * \param Yi The gas phase mole fractions.
   * \param cachedP0 An array of cached pressure values.
   */
  void printErrorStatus(double psi, double sum, double P, const std::vector<double> Yi, double cachedP0[]);
};
