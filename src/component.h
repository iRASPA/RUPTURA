#pragma once

#include <iostream>

#include "multi_site_isotherm.h"
/**
 * \brief Represents a chemical component in the simulation.
 *
 * The Component struct encapsulates the properties and behaviors of a chemical component within the system.
 * It includes identifiers, names, isotherm data, and parameters related to mass transfer and diffusion.
 */
struct Component
{
  /**
   * \brief Constructs a Component with an id and name.
   *
   * Initializes a Component using the provided identifier and name.
   *
   * \param i Identifier for the component.
   * \param n Name of the component.
   */
  Component(size_t i, std::string n) : id(i), name(n) {}

  /**
   * \brief Constructs a Component with specified parameters.
   *
   * Initializes a Component with the provided id, name, isotherms, initial mol-fraction, mass transfer coefficient,
   * diffusion coefficient, and an optional flag indicating if it is a carrier gas.
   *
   * \param _id Identifier for the component.
   * \param _name Name of the component.
   * \param _isotherms Vector of Isotherm objects associated with the component.
   * \param _Yi0 Initial gas phase mol-fraction [-].
   * \param _Kl Mass transfer coefficient [1/s].
   * \param _D Axial dispersion coefficient [m^2/s].
   * \param _isCarrierGas Optional flag indicating if this is the carrier gas (default is false).
   */
  Component(size_t _id, std::string _name, std::vector<Isotherm> _isotherms, double _Yi0, double _Kl, double _D,
            bool _isCarrierGas = false);

  size_t id;                   ///< Identifier of the component.
  std::string name{};          ///< Name of the component.
  std::string filename{};      ///< Filename associated with the component data.
  MultiSiteIsotherm isotherm;  ///< Isotherm information for the component.
  double Yi0;                  ///< Gas phase mol-fraction [-].
  double Kl;                   ///< Mass transfer coefficient [1/s].
  double D;                    ///< Axial dispersion coefficient [m^2/s].
  bool isCarrierGas{false};    ///< Flag indicating if this is the carrier gas.

  /**
   * \brief Prints the component information to the console.
   *
   * Outputs a string representation of the component to the standard output.
   */
  void print() const;

  /**
   * \brief Returns a string representation of the Component.
   *
   * Generates a string that includes the component's properties and parameters.
   *
   * \return A string representing the Component.
   */
  std::string repr() const;
};
