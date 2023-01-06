#include "component.h"

void Component::print(size_t i) const
{
  std::cout << "Component " << i << " id: " << id << " [" << name << "]:\n";
  if(isCarrierGas)
  {
    std::cout << "    carrier-gas\n";

    isotherm.print();
  }
  std::cout << "    mol-fraction in the gas:   " << Yi0 << " [-]\n";
  if(!isCarrierGas)
  {
    std::cout << "    mass-transfer coefficient: " << Kl << " [1/s]\n";
    std::cout << "    diffusion coefficient:     " << D << " [m^2/s]\n";

    isotherm.print();
  }
}
