#include "component.h"

#include <string>

#include "isotherm.h"

Component::Component(size_t _id, std::string _name, std::vector<Isotherm> _isotherms, double _Yi0, double _Kl,
                     double _D, bool _isCarrierGas)
    : id(_id), name(_name), Yi0(_Yi0), Kl(_Kl), D(_D), isCarrierGas(_isCarrierGas)
{
  isotherm.numberOfSites = _isotherms.size();
  for (Isotherm it : _isotherms)
  {
    isotherm.add(it);
  }
}

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
      s += "    mass-transfer coefficient: " + std::to_string(Kl) + " [1/s]\n";
      s += "    diffusion coefficient:     " + std::to_string(D) + " [m^2/s]\n";
      s += isotherm.repr();
    }
    return s;
}

std::vector<Component>& normalize_molfracs(std::vector<Component>& components)
{
  double total = 0;
  for (auto comp : components)
  {
    total += comp.Yi0;
  }
  for (auto comp : components)
  {
    comp.Yi0 /= total;
  }
  return components;
}
