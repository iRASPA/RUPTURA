#pragma once

#include <array>
#include <cstddef>
#include <iostream>
#include <string>
#include <tuple>
#include <vector>

#include "hash_combine.h"
#include "isotherm.h"

/**
 * \brief Represents a collection of adsorption isotherm sites.
 *
 * The MultiSiteIsotherm struct encapsulates multiple individual isotherm sites,
 * allowing the modeling of complex adsorption behavior as a sum of multiple isotherm contributions.
 * It manages a collection of Isotherm objects and provides methods to evaluate the combined adsorption values,
 * inverse computations, and parameter handling. This struct facilitates operations such as parameter randomization,
 * fitness evaluation, and generating representations for plotting.
 */
struct MultiSiteIsotherm
{
  size_t numberOfSites{0};        ///< The number of isotherm sites included in the model.
  std::vector<Isotherm> sites{};  ///< A vector containing the individual isotherm site objects.

  size_t numberOfParameters{0};  ///< The total number of parameters across all isotherm sites.
  std::vector<std::pair<size_t, size_t>>
      parameterIndices{};                    ///< Mapping from parameter index to site and parameter within site.
  std::vector<size_t> siteParameterIndex{};  ///< Indices indicating the starting parameter index for each site.

  /**
   * \brief Accesses the parameter at the specified index.
   *
   * Provides a reference to the parameter value corresponding to the given global parameter index, allowing for
   * modification.
   *
   * \param i The global parameter index.
   * \return Reference to the parameter value.
   */
  double &parameters(size_t i)
  {
    std::pair<size_t, size_t> index = parameterIndices[i];
    return sites[index.first].parameters[index.second];
  }

  /**
   * \brief Accesses the parameter at the specified index (const version).
   *
   * Provides a const reference to the parameter value corresponding to the given global parameter index.
   *
   * \param i The global parameter index.
   * \return Const reference to the parameter value.
   */
  const double &parameters(size_t i) const
  {
    std::pair<size_t, size_t> index = parameterIndices[i];
    return sites[index.first].parameters[index.second];
  }

  /**
   * \brief Adds an isotherm site to the collection.
   *
   * Incorporates a new Isotherm object into the MultiSiteIsotherm, updating the parameter indices and counts
   * accordingly.
   *
   * \param isotherm The Isotherm object to be added.
   */
  void add(const Isotherm &isotherm);

  /**
   * \brief Prints the string representation of the MultiSiteIsotherm.
   *
   * Outputs the result of repr() to the standard output.
   */
  void print() const;

  /**
   * \brief Generates a string representation of the MultiSiteIsotherm.
   *
   * Creates a string that includes the number of sites and the representations of each site.
   *
   * \return A string representing the MultiSiteIsotherm.
   */
  std::string repr() const;

  /**
   * \brief Generates a randomized version of the MultiSiteIsotherm.
   *
   * Creates a copy of the current MultiSiteIsotherm with randomized parameters for each site,
   * within the specified maximum loading.
   *
   * \param maximumLoading The maximum loading value used for randomization.
   * \return A randomized MultiSiteIsotherm object.
   */
  MultiSiteIsotherm randomized(double maximumLoading)
  {
    MultiSiteIsotherm copy(*this);
    for (size_t i = 0; i < numberOfSites; ++i)
    {
      copy.sites[i].randomize(maximumLoading);
    }
    return copy;
  }

  /**
   * \brief Sets the parameters of all isotherm sites.
   *
   * Updates the parameters of the isotherm sites using the provided vector of parameter values.
   *
   * \param params A vector containing the new parameter values.
   */
  void setParameters(std::vector<double> params);

  /**
   * \brief Retrieves all parameters of the isotherm sites.
   *
   * Collects the parameters from all isotherm sites into a single vector.
   *
   * \return A vector containing all parameter values.
   */
  std::vector<double> getParameters();

  /**
   * \brief Computes the total adsorption value at a given pressure.
   *
   * Sums the adsorption values from all isotherm sites at the specified pressure.
   *
   * \param pressure The pressure at which to evaluate the adsorption.
   * \return The total adsorption value.
   */
  inline double value(double pressure) const
  {
    double sum = 0.0;
    for (size_t i = 0; i < numberOfSites; ++i)
    {
      sum += sites[i].value(pressure);
    }
    return sum;
  }

  /**
   * \brief Computes the adsorption value for a specific site at a given pressure.
   *
   * Evaluates the adsorption value from a specified isotherm site at the given pressure.
   *
   * \param site The index of the isotherm site.
   * \param pressure The pressure at which to evaluate the adsorption.
   * \return The adsorption value for the specified site, or 0.0 if the site index is invalid.
   */
  inline double value(size_t site, double pressure) const
  {
    if (site < numberOfSites)
    {
      return sites[site].value(pressure);
    }
    return 0.0;
  }

  /**
   * \brief Computes the reduced grand potential at a given pressure.
   *
   * Sums the reduced grand potential contributions from all isotherm sites at the specified pressure.
   *
   * \param pressure The pressure at which to compute the reduced grand potential.
   * \return The total reduced grand potential.
   */
  inline double psiForPressure(double pressure) const
  {
    double sum = 0.0;
    for (size_t i = 0; i < numberOfSites; ++i)
    {
      sum += sites[i].psiForPressure(pressure);
    }
    return sum;
  }

  /**
   * \brief Computes the reduced grand potential for a specific site at a given pressure.
   *
   * Evaluates the reduced grand potential from a specified isotherm site at the given pressure.
   *
   * \param site The index of the isotherm site.
   * \param pressure The pressure at which to compute the reduced grand potential.
   * \return The reduced grand potential for the specified site, or 0.0 if the site index is invalid.
   */
  inline double psiForPressure(size_t site, double pressure) const
  {
    if (site < numberOfSites)
    {
      return sites[site].psiForPressure(pressure);
    }
    return 0.0;
  }

  /**
   * \brief Computes the inverse pressure corresponding to a given reduced grand potential.
   *
   * Calculates the inverse pressure (1/P) that corresponds to the specified total reduced grand potential
   * using a bisection algorithm.
   *
   * \param reduced_grand_potential The target reduced grand potential.
   * \param cachedP0 A reference to a cached pressure value for starting point optimization.
   * \return The inverse of the pressure corresponding to the reduced grand potential.
   */
  double inversePressureForPsi(double reduced_grand_potential, double &cachedP0) const;

  /**
   * \brief Computes the inverse pressure for a specific site corresponding to a given reduced grand potential.
   *
   * Calculates the inverse pressure (1/P) for a specified isotherm site that corresponds to the given reduced grand
   * potential.
   *
   * \param site The index of the isotherm site.
   * \param reduced_grand_potential The target reduced grand potential.
   * \param cachedP0 A reference to a cached pressure value for starting point optimization.
   * \return The inverse of the pressure for the specified site, or 0.0 if the site index is invalid.
   */
  double inversePressureForPsi(size_t site, double reduced_grand_potential, double &cachedP0) const
  {
    if (site < numberOfSites)
    {
      return sites[site].inversePressureForPsi(reduced_grand_potential, cachedP0);
    }
    return 0.0;
  }

  /**
   * \brief Evaluates the fitness of the MultiSiteIsotherm.
   *
   * Calculates a fitness value used in optimization, applying a penalty if any of the isotherm sites are unphysical.
   *
   * \return The fitness value, where higher values indicate worse fitness.
   */
  double fitness() const;

  /**
   * \brief Generates a gnuplot function string for plotting.
   *
   * Creates a string that represents the MultiSiteIsotherm as a gnuplot function, suitable for plotting purposes.
   *
   * \param s A character representing the parameter variable in the function string.
   * \return A string representing the MultiSiteIsotherm in gnuplot syntax.
   */
  std::string gnuplotFunctionString(char s) const;
};

namespace std
{
template <>
struct hash<MultiSiteIsotherm>
{
  size_t operator()(const MultiSiteIsotherm &k) const
  {
    std::size_t h = 0;
    for (const Isotherm &isotherm : k.sites)
    {
      for (size_t i = 0; i < isotherm.numberOfParameters; ++i)
      {
        hash_combine(h, isotherm.parameters[i]);
      }
    }
    return h;
  }
};
}  // namespace std
