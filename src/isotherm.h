#pragma once

#include <array>
#include <cstddef>
#include <vector>
#define _USE_MATH_DEFINES
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
#include <math.h>
#else
#include <cmath>
#endif
#include <iostream>
#include <string>

#include "random_numbers.h"
#include "special_functions.h"

/**
 * \brief Maximum number of terms supported in the isotherm calculations.
 */
constexpr size_t maxTerms = 5;

// Langmuir:
// parameter 0: K
// parameter 1: N
//
// Langmuir-Freundlich
// parameter 0: K
// parameter 1: N
// parameter 2: power

/**
 * \brief Represents an isotherm model for adsorption processes.
 *
 * The Isotherm struct encapsulates various isotherm models used to describe the adsorption equilibrium
 * between a fluid and a solid at a constant temperature. It supports different types of isotherm models,
 * such as Langmuir, Freundlich, BET, and others.
 */
struct Isotherm
{
  /**
   * \brief Enumeration of the different types of isotherm models.
   */
  enum class Type
  {
    Langmuir = 0,             ///< Langmuir isotherm model
    Anti_Langmuir = 1,        ///< Anti-Langmuir isotherm model
    BET = 2,                  ///< Brunauer–Emmett–Teller (BET) isotherm model
    Henry = 3,                ///< Henry's law isotherm model
    Freundlich = 4,           ///< Freundlich isotherm model
    Sips = 5,                 ///< Sips isotherm model
    Langmuir_Freundlich = 6,  ///< Langmuir-Freundlich isotherm model
    Redlich_Peterson = 7,     ///< Redlich-Peterson isotherm model
    Toth = 8,                 ///< Toth isotherm model
    Unilan = 9,               ///< Unilan isotherm model
    OBrien_Myers = 10,        ///< O'Brien and Myers isotherm model
    Quadratic = 11,           ///< Quadratic isotherm model
    Temkin = 12,              ///< Temkin isotherm model
    BingelWalton = 13         ///< Bingel and Walton isotherm model
  };

  /**
   * \brief Constructs an Isotherm with specified type and parameters.
   *
   * Initializes an Isotherm object with the given isotherm type, parameter values, and the number of parameters.
   *
   * \param type The type of the isotherm model.
   * \param values A vector of parameter values for the isotherm model.
   * \param numberOfValues The number of parameters.
   */
  Isotherm(Isotherm::Type type, const std::vector<double> &values, size_t numberOfValues);

  /**
   * \brief Constructs an Isotherm with specified type index and parameters.
   *
   * Initializes an Isotherm object using the type index, parameter values, and the number of parameters.
   *
   * \param t The index of the isotherm type.
   * \param values A vector of parameter values for the isotherm model.
   * \param numberOfValues The number of parameters.
   */
  Isotherm(size_t t, const std::vector<double> &values, size_t numberOfValues);

  Isotherm::Type type;             ///< The type of the isotherm model.
  std::vector<double> parameters;  ///< Parameter values for the isotherm model.
  size_t numberOfParameters;       ///< The number of parameters for the isotherm model.

  /**
   * \brief Prints a representation of the isotherm to standard output.
   */
  void print() const;

  /**
   * \brief Returns a string representation of the isotherm.
   *
   * \return A string representing the isotherm model and its parameters.
   */
  std::string repr() const;

  /**
   * \brief Computes the adsorption amount at a given pressure.
   *
   * Calculates the adsorption loading based on the isotherm model and parameters for the specified pressure.
   *
   * \param pressure The pressure at which to evaluate the isotherm.
   * \return The adsorption amount at the given pressure.
   */
  inline double value(double pressure) const
  {
    switch (type)
    {
      case Isotherm::Type::Langmuir:
      {
        double temp = parameters[1] * pressure;
        return parameters[0] * temp / (1.0 + temp);
      }
      case Isotherm::Type::Anti_Langmuir:
      {
        return parameters[0] * pressure / (1.0 - parameters[1] * pressure);
      }
      case Isotherm::Type::BET:
      {
        return parameters[0] * parameters[1] * pressure /
               ((1.0 - parameters[2] * pressure) * (1.0 - parameters[2] + parameters[1] * pressure));
      }
      case Isotherm::Type::Henry:
      {
        return parameters[0] * pressure;
      }
      case Isotherm::Type::Freundlich:
      {
        return parameters[0] * std::pow(pressure, 1.0 / parameters[1]);
      }
      case Isotherm::Type::Sips:
      {
        double temp = std::pow(parameters[1] * pressure, 1.0 / parameters[2]);
        return parameters[0] * temp / (1.0 + temp);
      }
      case Isotherm::Type::Langmuir_Freundlich:
      {
        double temp = parameters[1] * std::pow(pressure, parameters[2]);
        return parameters[0] * temp / (1.0 + temp);
      }
      case Isotherm::Type::Redlich_Peterson:
      {
        return parameters[0] * pressure / (1.0 + parameters[1] * std::pow(pressure, parameters[2]));
      }
      case Isotherm::Type::Toth:
      {
        double temp = parameters[1] * pressure;
        return parameters[0] * temp / std::pow(1.0 + std::pow(temp, parameters[2]), 1.0 / parameters[2]);
      }
      case Isotherm::Type::Unilan:
      {
        double temp1 = 1.0 + parameters[1] * std::exp(parameters[2]) * pressure;
        double temp2 = 1.0 + parameters[1] * std::exp(-parameters[2]) * pressure;
        return parameters[0] * (0.5 / parameters[2]) * std::log(temp1 / temp2);
      }
      case Isotherm::Type::OBrien_Myers:
      {
        double temp1 = parameters[1] * pressure;
        double temp2 = 1.0 + temp1;
        return parameters[0] *
               (temp1 / temp2 + parameters[2] * parameters[2] * temp1 * (1.0 - temp1) / (temp2 * temp2 * temp2));
      }
      case Isotherm::Type::Quadratic:
      {
        double temp1 = parameters[1] * pressure;
        double temp2 = parameters[2] * pressure * pressure;
        return parameters[0] * (temp1 + 2.0 * temp2) / (1.0 + temp1 + temp2);
      }
      case Isotherm::Type::Temkin:
      {
        double temp = parameters[1] * pressure;
        double temp1 = temp / (1.0 + temp);
        return parameters[0] * (temp1 + parameters[2] * temp1 * temp1 * (temp1 - 1.0));
      }
      case Isotherm::Type::BingelWalton:
      {
        return parameters[0] * (1.0 - std::exp(-(parameters[1] + parameters[2]) * pressure)) /
               (1.0 + (parameters[2] / parameters[1]) * std::exp(-(parameters[1] + parameters[2]) * pressure));
      }
      default:
        throw std::runtime_error("Error: unknown isotherm type");
    }
  }

  /**
   * \brief Computes the reduced grand potential (spreading pressure) at a given pressure.
   *
   * Calculates the reduced grand potential psi based on the isotherm model and parameters for the specified pressure.
   *
   * \param pressure The pressure at which to evaluate the reduced grand potential.
   * \return The reduced grand potential psi at the given pressure.
   */
  inline double psiForPressure(double pressure) const
  {
    switch (type)
    {
      case Isotherm::Type::Langmuir:
      {
        return parameters[0] * std::log(1.0 + parameters[1] * pressure);
      }
      case Isotherm::Type::Anti_Langmuir:
      {
        return -(parameters[0] / parameters[1]) * std::log(1.0 - parameters[1] * pressure);
      }
      case Isotherm::Type::BET:
      {
        return (parameters[0] * parameters[1]) *
               std::log((1.0 - parameters[2] + parameters[1] * pressure) /
                        ((1.0 - parameters[2]) * (1.0 - parameters[2] * pressure))) /
               (parameters[1] + parameters[2] - parameters[2] * parameters[2]);
      }
      case Isotherm::Type::Henry:
      {
        return parameters[0] * pressure;
      }
      case Isotherm::Type::Freundlich:
      {
        return parameters[0] * parameters[1] * std::pow(pressure, 1.0 / parameters[1]);
      }
      case Isotherm::Type::Sips:
      {
        return parameters[2] * parameters[0] * std::log(1.0 + std::pow(parameters[1] * pressure, 1.0 / parameters[2]));
      }
      case Isotherm::Type::Langmuir_Freundlich:
      {
        return (parameters[0] / parameters[2]) * std::log(1.0 + parameters[1] * std::pow(pressure, parameters[2]));
      }
      case Isotherm::Type::Redlich_Peterson:
      {
        if (parameters[1] * std::pow(pressure, parameters[2]) < 1.0)
        {
          return parameters[0] * pressure *
                 hypergeometric2F1(1.0, 1.0 / parameters[2], 1.0 + 1.0 / parameters[2],
                                   -parameters[1] * std::pow(pressure, parameters[2]));
        }
        else
        {
          double prefactor = parameters[0] / parameters[2];
          double temp = M_PI / (std::pow(parameters[1], 1.0 / parameters[2]) * std::sin(M_PI * 1.0 / parameters[2]));

          double term1 = -1.0 / (parameters[1] * std::pow(pressure, parameters[2]));
          double numerator = 1.0;
          double sum = 0.0;
          // quickly converging sum
          for (size_t k = 1; k <= 15; k++)
          {
            numerator *= term1;
            sum += numerator / (static_cast<double>(k) * parameters[2] - 1.0);
          }
          return prefactor * (temp + pressure * parameters[2] * sum);
        }
      }
      case Isotherm::Type::Toth:
      {
        double temp = parameters[1] * pressure;
        double theta = temp / std::pow(1.0 + std::pow(temp, parameters[2]), 1.0 / parameters[2]);
        double theta_pow = std::pow(theta, parameters[2]);
        double psi = parameters[0] * (theta - (theta / parameters[2]) * std::log(1.0 - theta_pow));

        // use the first 100 terms of the sum
        double temp1 = parameters[0] * theta;
        double temp2 = 0.0;
        for (size_t k = 1; k <= 100; ++k)
        {
          temp1 *= theta_pow;
          temp2 += parameters[2];
          psi -= temp1 / (temp2 * (temp2 + 1.0));
        }

        return psi;
      }
      case Isotherm::Type::Unilan:
      {
        return (0.5 * parameters[0] / parameters[2]) * (li2(-parameters[1] * std::exp(-parameters[2]) * pressure) -
                                                        li2(-parameters[1] * std::exp(parameters[2]) * pressure));
      }
      case Isotherm::Type::OBrien_Myers:
      {
        double temp1 = parameters[1] * pressure;
        double temp2 = 1.0 + temp1;
        return parameters[0] * (std::log(temp2) + 0.5 * parameters[2] * parameters[2] * temp1 / (temp2 * temp2));
      }
      case Isotherm::Type::Quadratic:
      {
        double temp1 = parameters[1] * pressure;
        double temp2 = parameters[2] * pressure * pressure;
        return parameters[0] * std::log(1.0 + temp1 + temp2);
      }
      case Isotherm::Type::Temkin:
      {
        double temp = parameters[1] * pressure;
        double temp1 = temp / (1.0 + temp);
        return parameters[0] * (std::log(1.0 + temp) - 0.5 * parameters[2] * temp1 * temp1);
      }
      case Isotherm::Type::BingelWalton:
      {
        double start = 1e-14;
        size_t max_steps = 100;
        double acc = 1e-6;

        // Romberg integration: https://en.wikipedia.org/wiki/Romberg%27s_method
        std::vector<double> R1(max_steps), R2(max_steps);                       // buffers
        double *Rp = &R1[0], *Rc = &R2[0];                                      // Rp is previous row, Rc is current row
        double h = pressure - start;                                            // step size
        Rp[0] = (value(start) / start + value(pressure) / pressure) * h * 0.5;  // first trapezoidal step

        for (size_t i = 1; i < max_steps; ++i)
        {
          h /= 2.0;
          double c = 0;
          size_t ep = size_t{1} << (i - 1);  // 2^(n-1)
          for (size_t j = 1; j <= ep; ++j)
          {
            c += value(start + static_cast<double>(2 * j - 1) * h) / (start + static_cast<double>(2 * j - 1) * h);
          }
          Rc[0] = h * c + 0.5 * Rp[0];  // R(i,0)

          for (size_t j = 1; j <= i; ++j)
          {
            double n_k = std::pow(4, j);
            Rc[j] = (n_k * Rc[j - 1] - Rp[j - 1]) / (n_k - 1);  // compute R(i,j)
          }

          if (i > 1 && std::fabs(Rp[i - 1] - Rc[i]) < acc)
          {
            return Rc[i];
          }

          double *rt = Rp;
          Rp = Rc;
          Rc = rt;
        }
        return Rp[max_steps - 1];  // return our best guess
      }
      default:
        throw std::runtime_error("Error: unknown isotherm type");
    }
  }

  /**
   * \brief Computes the inverse pressure corresponding to a given reduced grand potential psi.
   *
   * Calculates the pressure that corresponds to the specified reduced grand potential psi using the isotherm model.
   * This function may cache intermediate results to improve performance in repeated calculations.
   *
   * \param reduced_grand_potential The reduced grand potential psi.
   * \param cachedP0 A reference to a cached pressure value used to initialize the calculation.
   * \return The inverse of the pressure corresponding to the given psi.
   */
  inline double inversePressureForPsi(double reduced_grand_potential, double &cachedP0) const
  {
    switch (type)
    {
      case Isotherm::Type::Langmuir:
      {
        double denominator = std::exp(reduced_grand_potential / parameters[0]) - 1.0;
        return parameters[1] / denominator;
      }
      case Isotherm::Type::Anti_Langmuir:
      {
        double denominator = 1.0 - std::exp(-parameters[1] * reduced_grand_potential / parameters[0]);
        return parameters[1] / denominator;
      }
      case Isotherm::Type::Henry:
      {
        return parameters[0] / reduced_grand_potential;
      }
      case Isotherm::Type::Freundlich:
      {
        return std::pow((parameters[0] * parameters[1]) / reduced_grand_potential, parameters[1]);
      }
      case Isotherm::Type::Sips:
      {
        return parameters[1] /
               std::pow((std::exp(reduced_grand_potential / (parameters[2] * parameters[0])) - 1.0), parameters[2]);
      }
      case Isotherm::Type::Langmuir_Freundlich:
      {
        double denominator = std::exp(reduced_grand_potential * parameters[2] / parameters[0]) - 1.0;
        return std::pow(parameters[1] / denominator, 1.0 / parameters[2]);
      }
      default:
      {
        const double tiny = 1.0e-15;

        // from here on, work with pressure, and return 1.0 / pressure at the end of the routine
        double p_start;
        if (cachedP0 <= 0.0)
        {
          p_start = 5.0;
        }
        else
        {
          // use the last value of Pi0
          p_start = cachedP0;
        }

        // use bisection algorithm
        double s = psiForPressure(p_start);

        size_t nr_steps = 0;
        double left_bracket = p_start;
        double right_bracket = p_start;

        if (s < reduced_grand_potential)
        {
          // find the bracket on the right
          do
          {
            right_bracket *= 2.0;
            s = psiForPressure(right_bracket);

            ++nr_steps;
            if (nr_steps > 100000)
            {
              std::cout << "reduced_grand_potential: " << reduced_grand_potential << std::endl;
              std::cout << "psi: " << s << std::endl;
              std::cout << "p_start: " << p_start << std::endl;
              std::cout << "Left bracket: " << left_bracket << std::endl;
              std::cout << "Right bracket: " << right_bracket << std::endl;
              throw std::runtime_error(
                  "Error (Inverse bisection): initial bracketing (for sum < 1) does NOT converge\n");
            }
          } while (s < reduced_grand_potential);
        }
        else
        {
          // find the bracket on the left
          do
          {
            left_bracket *= 0.5;
            s = psiForPressure(left_bracket);

            ++nr_steps;
            if (nr_steps > 100000)
            {
              std::cout << "reduced_grand_potential: " << reduced_grand_potential << std::endl;
              std::cout << "psi: " << s << std::endl;
              std::cout << "p_start: " << p_start << std::endl;
              std::cout << "Left bracket: " << left_bracket << std::endl;
              std::cout << "Right bracket: " << right_bracket << std::endl;
              throw std::runtime_error(
                  "Error (Inverse bisection): initial bracketing (for sum > 1) does NOT converge\n");
            }
          } while (s > reduced_grand_potential);
        }

        do
        {
          double middle = 0.5 * (left_bracket + right_bracket);
          s = psiForPressure(middle);

          if (s > reduced_grand_potential)
            right_bracket = middle;
          else
            left_bracket = middle;

          ++nr_steps;
          if (nr_steps > 100000)
          {
            std::cout << "Left bracket: " << left_bracket << std::endl;
            std::cout << "Right bracket: " << right_bracket << std::endl;
            throw std::runtime_error("Error (Inverse bisection): initial bracketing (for sum < 1) does NOT converge\n");
          }
        } while (std::abs(left_bracket - right_bracket) / std::abs(left_bracket + right_bracket) > tiny);

        double middle = 0.5 * (left_bracket + right_bracket);

        //  Store the last value of Pi0
        cachedP0 = middle;

        return 1.0 / middle;
      }
    }
  }

  /**
   * \brief Randomizes the isotherm parameters within specified bounds.
   *
   * Sets the isotherm parameters to random values within physically meaningful ranges,
   * based on the specified maximum loading.
   *
   * \param maximumLoading The maximum adsorption loading used to scale the randomized parameters.
   */
  void randomize(double maximumLoading);

  /**
   * \brief Checks if the isotherm parameters are physically meaningful.
   *
   * Determines whether the isotherm parameters are within acceptable physical ranges.
   *
   * \return True if the parameters are unphysical; false otherwise.
   */
  bool isUnphysical() const;

  /**
   * \brief Generates a Gnuplot-compatible function string for the isotherm.
   *
   * Creates a string representing the isotherm function that can be used in Gnuplot scripts,
   * using the specified character and index for parameter substitution.
   *
   * \param s The character representing the parameter array in Gnuplot (e.g., 'c' or 'p').
   * \param i The starting index for the parameters.
   * \return A string representing the isotherm function for Gnuplot.
   */
  std::string gnuplotFunctionString(char s, size_t i) const;
};
