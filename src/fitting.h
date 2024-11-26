#pragma once

#include <array>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

#include "component.h"
#include "inputreader.h"
#include "multi_site_isotherm.h"

#ifdef PYBUILD
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
namespace py = pybind11;
#endif  // PYBUILD

/**
 * \brief Class for fitting isotherm models to adsorption data.
 *
 * The Fitting class implements genetic algorithms and the Nelder-Mead simplex method
 * to optimize isotherm parameters based on input data.
 */
struct Fitting
{
  /**
   * \brief Structure representing an individual in the genetic algorithm.
   */
  struct DNA
  {
    /**
     * \brief Constructs a DNA object with specified genotype, phenotype, and fitness.
     * \param g Genotype represented as a string.
     * \param p Phenotype represented as a MultiSiteIsotherm.
     * \param f Fitness value of the DNA.
     */
    DNA(std::string g, MultiSiteIsotherm p, double f)
        : genotype(g), phenotype(p), fitness(f), hash(std::hash<std::string>{}(g))
    {
    }
    DNA() noexcept = default;

    std::string genotype;         ///< Genotype represented as a bitstring.
    MultiSiteIsotherm phenotype;  ///< Phenotype of the individual.
    double fitness;               ///< Fitness value.
    size_t hash;                  ///< Hash value for uniqueness.
  };

  /**
   * \brief Enum representing pressure scale options.
   */
  enum class PressureScale
  {
    Log = 0,    ///< Logarithmic pressure scale.
    Normal = 1  ///< Linear pressure scale.
  };

  /**
   * \brief Constructs a Fitting object from input parameters.
   * \param inputreader InputReader containing simulation parameters.
   */
  Fitting(const InputReader &inputreader);

  /**
   * \brief Reads data for a specific component.
   * \param ID Index of the component.
   */
  void readData(size_t ID);

  /**
   * \brief Prints the solution for a specific component.
   * \param ID Index of the component.
   */
  void printSolution(size_t ID);

  /**
   * \brief Runs the fitting process for all components.
   */
  void run();

  /**
   * \brief Creates plot scripts for a specific component.
   * \param citizen DNA of the optimized individual.
   * \param ID Index of the component.
   */
  void createPlotScripts(const DNA &citizen, size_t ID);

  /**
   * \brief Creates a master plot script.
   */
  void createPlotScript();

#ifdef PYBUILD
  /**
   * \brief Constructs a Fitting object for Python integration.
   * \param _displayName Display name of the fitting session.
   * \param _components Vector of components to fit.
   * \param _fullData Full dataset for all components.
   * \param _pressureScale Pressure scale type.
   */
  Fitting(std::string _displayName, std::vector<Component> _components, std::vector<std::vector<double>> _fullData,
          size_t _pressureScale);

  /**
   * \brief Extracts data slices for a specific component.
   * \param ID Index of the component.
   */
  void sliceData(size_t ID);

  std::vector<std::vector<double>> fullData;  ///< Full dataset for all components.

  /**
   * \brief Computes the fitting parameters.
   * \return Vector of optimized parameters.
   */
  std::vector<double> compute();

  /**
   * \brief Evaluates the fitted isotherm models.
   * \return NumPy array of evaluated values.
   */
  py::array_t<double> evaluate();

#endif  // PYBUILD

  /**
   * \brief Generates a new individual for the genetic algorithm.
   * \param ID Index of the component.
   * \return A new DNA object.
   */
  DNA newCitizen(size_t ID);

  /**
   * \brief Updates the fitness of a DNA object.
   * \param citizen DNA object to update.
   */
  void updateCitizen(DNA &citizen);

  /**
   * \brief Calculates the fitness of a phenotype.
   * \param phenotype Phenotype to evaluate.
   * \return Fitness value.
   */
  double fitness(const MultiSiteIsotherm &phenotype);

  /**
   * \brief Calculates the correlation coefficient R.
   * \param phenotype Phenotype to evaluate.
   * \return Correlation coefficient.
   */
  double RCorrelation(const MultiSiteIsotherm &phenotype);

  /**
   * \brief Measures biodiversity in the population.
   * \param citizens Vector of DNA objects.
   * \return Biodiversity metric.
   */
  size_t biodiversity(const std::vector<DNA> &citizens);

  /**
   * \brief Simulates a nuclear disaster to introduce new genetic material.
   * \param ID Index of the component.
   */
  void nuclearDisaster(size_t ID);

  /**
   * \brief Performs elitism selection in the genetic algorithm.
   */
  void elitism();

  /**
   * \brief Mutates a DNA object.
   * \param Mutant DNA object to mutate.
   */
  void mutate(DNA &Mutant);

  /**
   * \brief Performs crossover between individuals.
   * \param ID Index of the component.
   * \param s1 Start index for children.
   * \param s2 End index for children.
   * \param i1 Start index for parent1.
   * \param i2 End index for parent1.
   * \param j1 Start index for parent2.
   * \param j2 End index for parent2.
   */
  void crossover(size_t ID, size_t s1, size_t s2, size_t i1, size_t i2, size_t j1, size_t j2);

  /**
   * \brief Randomly selects parent indices.
   * \param kk1 Minimum index for parent1.
   * \param kk2 Maximum index for parent1.
   * \param jj1 Minimum index for parent2.
   * \param jj2 Maximum index for parent2.
   * \param ii1 Output index for parent1.
   * \param ii2 Output index for parent2.
   */
  void chooseRandomly(size_t kk1, size_t kk2, size_t jj1, size_t jj2, size_t &ii1, size_t &ii2);

  /**
   * \brief Mates individuals to produce the next generation.
   * \param ID Index of the component.
   */
  void mate(size_t ID);

  /**
   * \brief Sorts the population by fitness.
   */
  void sortByFitness();

  /**
   * \brief Writes information about a DNA object.
   * \param citizen Index of the citizen in the population.
   * \param id Component index.
   * \param step Current optimization step.
   * \param variety Biodiversity metric.
   * \param fullfilledCondition Condition fulfillment counter.
   */
  void writeCitizen(size_t citizen, size_t id, size_t step, size_t variety, size_t fullfilledCondition);

  /**
   * \brief Fits the isotherm model to data for a specific component.
   * \param ID Index of the component.
   * \return Best DNA object found.
   */
  DNA fit(size_t ID);

  /**
   * \brief Optimizes a DNA object using the Nelder-Mead simplex method.
   * \param citizen DNA object to optimize.
   * \param scale Scaling factor for the simplex.
   * \return Optimized DNA object.
   */
  const DNA simplex(DNA citizen, double scale);

  size_t Ncomp;                                     ///< Number of components.
  std::vector<Component> components;                ///< Components involved in fitting.
  std::string displayName;                          ///< Display name for the fitting session.
  std::vector<std::string> componentName;           ///< Names of the components.
  std::vector<std::string> filename;                ///< Filenames for input data.
  std::vector<MultiSiteIsotherm> isotherms;         ///< Isotherm models for each component.
  size_t columnPressure{0};                         ///< Column index for pressure data.
  size_t columnLoading{1};                          ///< Column index for loading data.
  size_t columnError{2};                            ///< Column index for error data.
  double maximumLoading{0.0};                       ///< Maximum loading observed.
  PressureScale pressureScale{PressureScale::Log};  ///< Pressure scale type.

  std::vector<std::pair<double, double>> rawData;  ///< Raw data points (pressure, loading).

  bool fittingFlag{false};             ///< Flag indicating if fitting is active.
  bool physicalConstrainsFlag{false};  ///< Flag for physical constraints.
  bool seedFlag{false};                ///< Flag for seed initialization.
  bool pressureRangeFlag{false};       ///< Flag for pressure range usage.
  bool refittingFlag{false};           ///< Flag for refitting.

  std::pair<double, double> pressureRange;     ///< Range of pressures.
  std::pair<double, double> logPressureRange;  ///< Logarithmic pressure range.

  size_t GA_Size;             ///< Genetic algorithm population size.
  double GA_MutationRate;     ///< Mutation rate in genetic algorithm.
  double GA_EliteRate;        ///< Proportion of elite individuals.
  double GA_MotleyCrowdRate;  ///< Proportion of diverse individuals.
  double GA_DisasterRate;     ///< Rate of introducing new genetic material.
  size_t GA_Elitists;         ///< Number of elite individuals.
  size_t GA_Motleists;        ///< Number of diverse individuals.

  std::vector<DNA> popAlpha;   ///< First population buffer.
  std::vector<DNA> popBeta;    ///< Second population buffer.
  std::vector<DNA> &parents;   ///< Reference to current parent population.
  std::vector<DNA> &children;  ///< Reference to current child population.
};
