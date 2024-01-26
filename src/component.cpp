#include "component.h"

std::string Component::repr() const
{
  std::string s;
  s += "Component id: " + std::to_string(id) + " [" + name + "]:\n";
  if (isCarrierGas)
  {
    s += "    carrier-gas\n";
    s += isotherm.repr();
  }
    s += "    mol-fraction in the gas:   " + std::to_string(Yi0) + " [-]\n";
    if (!isCarrierGas)
    {
      s += "    mas-transfer coefficient: " + std::to_string(Kl) + " [1/s]\n";
      s += "    diffusion coefficient:     " + std::to_string(D) + " [m^2/s]\n";
      s += isotherm.repr();
    }
    return s;
}
