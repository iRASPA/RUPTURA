#include "isotherm.h"

#include <cstdlib>

Isotherm::Isotherm(Isotherm::Type t, const std::vector<double> &values, size_t numberOfValues):
    type(t),
    parameters(values),
    numberOfParameters(numberOfValues)
{
}

Isotherm::Isotherm(std::string t, const std::vector<double> &values, size_t numberOfValues)
    : parameters(values), numberOfParameters(numberOfValues)
{
  auto itr = Isotherm::isotherm_map.find(t);
  if (itr != isotherm_map.end())
  {
    type = itr->second;
  }
  else
  {
    std::cout << "Error: unknown isotherm type (" << t << ").";
    exit(0);
  }
}

std::string Isotherm::repr() const
{
  std::string s;
  switch(type)
  {
    case Isotherm::Type::Langmuir:
    {
      s += "    Langmuir isotherm\n";
      s += "        q_sat: " + std::to_string(parameters[0]) + "\n";
      s += "        b:     " + std::to_string(parameters[1]) + "\n";
      break;
    }
    case Isotherm::Type::Anti_Langmuir:
    {
      s += "    Anti-Langmuir isotherm\n";
      s += "        a:     " + std::to_string(parameters[0]) + "\n";
      s += "        b:     " + std::to_string(parameters[1]) + "\n";
      break;
    }
    case Isotherm::Type::BET:
    {
      s += "    BET isotherm\n";
      s += "        q_sat: " + std::to_string(parameters[0]) + "\n";
      s += "        b:     " + std::to_string(parameters[1]) + "\n";
      s += "        c:     " + std::to_string(parameters[2]) + "\n";
      break;
    }
    case Isotherm::Type::Henry:
    {
      s += "    Henry isotherm\n";
      s += "        a:     " + std::to_string(parameters[0]) + "\n";
      break;
    }
    case Isotherm::Type::Freundlich:
    {
      s += "    Freundlich isotherm\n";
      s += "        a:     " + std::to_string(parameters[0]) + "\n";
      s += "        nu:    " + std::to_string(parameters[1]) + "\n";
      break;
    }
    case Isotherm::Type::Sips:
    {
      s += "    Sips isotherm\n";
      s += "        q_sat: " + std::to_string(parameters[0]) + "\n";
      s += "        b:     " + std::to_string(parameters[1]) + "\n";
      s += "        nu:    " + std::to_string(parameters[2]) + "\n";
      break;
    }
    case Isotherm::Type::Langmuir_Freundlich:
    {
      s += "    Langmuir-Freundlich isotherm\n";
      s += "        q_sat: " + std::to_string(parameters[0]) + "\n";
      s += "        b:     " + std::to_string(parameters[1]) + "\n";
      s += "        nu:    " + std::to_string(parameters[2]) + "\n";
      break;
    }
    case Isotherm::Type::Redlich_Peterson:
    {
      s += "    Redlich-Peterson isotherm\n";
      s += "        a:     " + std::to_string(parameters[0]) + "\n";
      s += "        b:     " + std::to_string(parameters[1]) + "\n";
      s += "        nu:    " + std::to_string(parameters[2]) + "\n";
      break;
    }
    case Isotherm::Type::Toth:
    {
      s += "    Toth isotherm\n";
      s += "        q_sat: " + std::to_string(parameters[0]) + "\n";
      s += "        b:     " + std::to_string(parameters[1]) + "\n";
      s += "        nu:    " + std::to_string(parameters[2]) + "\n";
      break;
    }
    case Isotherm::Type::Unilan:
    {
      s += "    Unilan isotherm\n";
      s += "        q_sat: " + std::to_string(parameters[0]) + "\n";
      s += "        b:     " + std::to_string(parameters[1]) + "\n";
      s += "        eta:   " + std::to_string(parameters[2]) + "\n";
      break;
    }
    case Isotherm::Type::OBrien_Myers:
    {
      s += "    O'Brian & Myers isotherm\n";
      s += "        q_sat: " + std::to_string(parameters[0]) + "\n";
      s += "        b:     " + std::to_string(parameters[1]) + "\n";
      s += "        sigma: " + std::to_string(parameters[2]) + "\n";
      break;
    }
    case Isotherm::Type::Quadratic:
    {
      s += "    Quadratic isotherm\n";
      s += "        q_sat: " + std::to_string(parameters[0]) + "\n";
      s += "        b:     " + std::to_string(parameters[1]) + "\n";
      s += "        c:     " + std::to_string(parameters[2]) + "\n";
      break;
    }
    case Isotherm::Type::Temkin:
    {
      s += "    Temkin isotherm\n";
      s += "        q_sat: " + std::to_string(parameters[0]) + "\n";
      s += "        b:     " + std::to_string(parameters[1]) + "\n";
      s += "        c:     " + std::to_string(parameters[2]) + "\n";
      break;
    }
    case Isotherm::Type::BingelWalton:
    {
      s += "    Bingel&Walton isotherm\n";
      s += "        q_sat: " + std::to_string(parameters[0]) + "\n";
      s += "        a:     " + std::to_string(parameters[1]) + "\n";
      s += "        b:     " + std::to_string(parameters[2]) + "\n";
      break;
    }
    default:
      break;
  }
  return s;
}

bool Isotherm::isUnphysical() const
{
  switch(type)
  {
    case Isotherm::Type::Langmuir:
    {
      if(parameters[0] < 0 || parameters[0] > 1.0e20 || parameters[1] < 0.0 || parameters[1] > 1.0e10 ) return true;
      return false;
    }
    case Isotherm::Type::Anti_Langmuir:
    {
      if(parameters[0] < 0 || parameters[0] > 1.0e20 || parameters[1] < 0.0 || parameters[1] > 1.0e10 ) return true;
      return false;
    }
    case Isotherm::Type::BET:
    {
      return false;
    }
    case Isotherm::Type::Henry:
    {
      if(parameters[0] < 0.0) return true;
      return false;
    }
    case Isotherm::Type::Freundlich:
    {
      if(parameters[0] < 0.0 || parameters[1] < 0.0 || parameters[2] < 0.0) return true;
      return false;
    }
    case Isotherm::Type::Sips:
    {
      if(parameters[0] < 0 || parameters[0] > 1.0e20 || parameters[1] < 0.0 || parameters[1] > 1.0e10 || parameters[2] < 0.0 || parameters[2] > 100.0 ) return true;
      return false;
    }
    case Isotherm::Type::Langmuir_Freundlich:
    {
      if(parameters[0] < 1.0e-20 || parameters[0] > 1.0e20 || parameters[1] < 0.0 || parameters[1] > 1.0e10 || parameters[2] < 0.0 || parameters[2] > 100.0 ) return true;
      return false;
    }
    case Isotherm::Type::Redlich_Peterson:
    {
      if(parameters[0] < 0.0 || parameters[1] < 0.0) return true;
      return false;
    }
    case Isotherm::Type::Toth:
    {
      if(parameters[0] < 0 || parameters[1] < 0.0 || parameters[2] < 0.0 || parameters[2] > 100.0 ) return true;
      return false;
    }
    case Isotherm::Type::Unilan:
    {
      if(parameters[0] < 0.0 || parameters[1] < 0.0 || parameters[2] < 0.0) return true;
      return false;
    }
    case Isotherm::Type::OBrien_Myers:
    {
      if(parameters[0] < 0.0 || parameters[1] < 0.0 || parameters[2] < 0.0) return true;
      return false;
    }
    case Isotherm::Type::Quadratic:
    {
      if(parameters[0] < 0.0 || parameters[1] < 0.0 || parameters[2] < 0.0) return true;
      return false;
    }
    case Isotherm::Type::Temkin:
    {
      if(parameters[0] <= 0.0 || parameters[1] < 0.0 || parameters[2] < 0.0) return true;
      return false;
    }
    case Isotherm::Type::BingelWalton:
    {
      if(parameters[0] <= 0.0 || (parameters[1] + parameters[2]) < 1e-3) return true;
      return false;
    }
    default:
      throw std::runtime_error("Error: unkown isotherm type");
  }
}

void Isotherm::randomize(double maximumLoading)
{
  switch(type)
  {
    case Isotherm::Type::Langmuir:
    {
      parameters[0] = 1.1 * maximumLoading * RandomNumber::Uniform();
      parameters[1] = std::pow(RandomNumber::Uniform(), 10.0 * 2.0 * (RandomNumber::Uniform() - 1.0));
      break;
    }
    case Isotherm::Type::Anti_Langmuir:
    {
      parameters[0] = 1.1 * maximumLoading * RandomNumber::Uniform();
      parameters[1] = std::pow(RandomNumber::Uniform(), 10.0 * 2.0 * (RandomNumber::Uniform() - 1.0));
      break;
    }
    case Isotherm::Type::BET:
    {
      parameters[0] = 10.0 * RandomNumber::Uniform();
      parameters[1] = 10.0 * RandomNumber::Uniform();
      parameters[2] = 10.0 * RandomNumber::Uniform();
      break;
    }
    case Isotherm::Type::Henry:
    {
      parameters[0] = std::pow(RandomNumber::Uniform(), 10.0 * 2.0 * (RandomNumber::Uniform() - 1.0));
      break;
    }
    case Isotherm::Type::Freundlich:
    {
      parameters[0] = std::pow(RandomNumber::Uniform(), 10.0 * 2.0 * (RandomNumber::Uniform() - 1.0));
      parameters[1] = 0.1 + 2.0 * RandomNumber::Uniform();
      break;
    }
    case Isotherm::Type::Sips:
    {
      parameters[0] = 1.1 * maximumLoading * RandomNumber::Uniform();
      parameters[1] = std::pow(RandomNumber::Uniform(), 10.0 * 2.0 * (RandomNumber::Uniform() - 1.0));
      parameters[2] = 0.1 + 2.0 * RandomNumber::Uniform();
      break;
    }
    case Isotherm::Type::Langmuir_Freundlich:
    {
      parameters[0] = 1.1 * maximumLoading * RandomNumber::Uniform();
      parameters[1] = std::pow(RandomNumber::Uniform(), 10.0 * 2.0 * (RandomNumber::Uniform() - 1.0));
      parameters[2] = 0.1 + 2.0 * RandomNumber::Uniform();
      break;
    }
    case Isotherm::Type::Redlich_Peterson:
    {
      parameters[0] = 1.1 * maximumLoading * RandomNumber::Uniform();
      parameters[1] = std::pow(RandomNumber::Uniform(), 10.0 * 2.0 * (RandomNumber::Uniform() - 1.0));
      parameters[2] = 0.1 + 2.0 * RandomNumber::Uniform();
      break;
    }
    case Isotherm::Type::Toth:
    {
      parameters[0] = 1.1 * maximumLoading * RandomNumber::Uniform();
      parameters[1] = std::pow(RandomNumber::Uniform(), 10.0 * 2.0 * (RandomNumber::Uniform() - 1.0));
      parameters[2] = 0.1 + 2.0 * RandomNumber::Uniform();
      break;
    }
    case Isotherm::Type::Unilan:
    {
      parameters[0] = 1.1 * maximumLoading * RandomNumber::Uniform();
      parameters[1] = std::pow(RandomNumber::Uniform(), 10.0 * 2.0 * (RandomNumber::Uniform() - 1.0));
      parameters[2] = 0.1 + 2.0 * RandomNumber::Uniform();
      break;
    }
    case Isotherm::Type::OBrien_Myers:
    {
      parameters[0] = 1.1 * maximumLoading * RandomNumber::Uniform();
      parameters[1] = std::pow(RandomNumber::Uniform(), 10.0 * 2.0 * (RandomNumber::Uniform() - 1.0));
      parameters[2] = 0.1 + 2.0 * RandomNumber::Uniform();
      break;
    }
    case Isotherm::Type::Quadratic:
    {
      parameters[0] = 2.1 * maximumLoading * RandomNumber::Uniform();
      parameters[1] = 10.0 * RandomNumber::Uniform();
      parameters[2] = 10.0 * RandomNumber::Uniform();
      break;
    }
    case Isotherm::Type::Temkin:
    {
      parameters[0] = 1.1 * maximumLoading * RandomNumber::Uniform();
      parameters[1] = std::pow(RandomNumber::Uniform(), 10.0 * 2.0 * (RandomNumber::Uniform() - 1.0));
      parameters[2] = 0.1 + 2.0 * RandomNumber::Uniform();
      break;
    }
    case Isotherm::Type::BingelWalton:
    {
      parameters[0] = 1.1 * maximumLoading * RandomNumber::Uniform();
      parameters[1] = 0.1 + 2.0 * RandomNumber::Uniform();
      parameters[2] = 0.1 + 2.0 * RandomNumber::Uniform();
      break;
    }
    default:
      throw std::runtime_error("Error: unkown isotherm type");
  }
}

std::string Isotherm::gnuplotFunctionString(char c, size_t i) const
{
  char stringBuffer[1024];

  switch(type)
  {
    case Isotherm::Type::Langmuir:
    {
      snprintf(stringBuffer, 1024, "%c[%ld]*%c[%ld]*x/(1.0+%c[%ld]*x)", c, i, c, i+1, c, i+1);
      return stringBuffer;
    }
    case Isotherm::Type::Anti_Langmuir:
    {
      snprintf(stringBuffer, 1024, "%c[%ld]*x/(1.0-%c[%ld]*x)", c, i, c, i+1);
      return stringBuffer;
    }
    case Isotherm::Type::BET:
    {
      snprintf(stringBuffer, 1024, "%c[%ld]*%c[%ld]*x/((1.0-%c[%ld]*x)*(1.0-%c[%ld]+%c[%ld]*x))", 
              c, i, c, i+1, c, i+2, c, i+2, c, i+1);
      return stringBuffer;
    }
    case Isotherm::Type::Henry:
    {
      snprintf(stringBuffer, 1024, "%c[%ld]*x", c, i);
      return stringBuffer;
    }
    case Isotherm::Type::Freundlich:
    {
      snprintf(stringBuffer, 1024, "%c[%ld]*x**[%ld]", c, i, i+1);
      return stringBuffer;
    }
    case Isotherm::Type::Sips:
    {
      snprintf(stringBuffer, 1024, "%c[%ld]*((%c[%ld]*x)**(1.0/%c[%ld]))/(1.0+(%c[%ld]*x)**(1.0/%c[%ld]))",
              c, i, c, i+1, c, i+2, c, i+1, c, i+2);
      return stringBuffer;
    }
    case Isotherm::Type::Langmuir_Freundlich:
    {
      snprintf(stringBuffer, 1024, "%c[%ld]*%c[%ld]*x**%c[%ld]/(1.0+%c[%ld]*x**%c[%ld])", 
              c, i, c, i+1, c, i+2, c, i +1, c, i+2);
      return stringBuffer;
    }
    case Isotherm::Type::Redlich_Peterson:
    {
      snprintf(stringBuffer, 1024, "%c[%ld]*x/(1.0+%c[%ld]*x**%c[%ld])", c, i, c, i+1, c, i+2);
      return stringBuffer;
    }
    case Isotherm::Type::Toth:
    {
      snprintf(stringBuffer, 1024, "%c[%ld]*%c[%ld]*x/((1.0+(%c[%ld]*x)**%c[%ld])**(1.0/%c[%ld]))",
              c, i, c, i+1, c, i+1, c, i+2, c, i+2);
      return stringBuffer;
    }
    case Isotherm::Type::Unilan:
    {
      snprintf(stringBuffer, 1024, "(%c[%ld]/(2.0*%c[%ld]))*log((1.0+%c[%ld]*exp(%c[%ld])*x)/(1.0+%c[%ld]*exp(-%c[%ld])*x))",
              c, i, c, i+2, c, i+1, c, i+2, c, i+1, c, i+2);
      return stringBuffer;
    }
    case Isotherm::Type::OBrien_Myers:
    {
      snprintf(stringBuffer, 1024, "%c[%ld]*(%c[%ld]*x/(1.0+%c[%ld]*x) + (%c[%ld]**2)*%c[%ld]*x*(1.0-%c[%ld]*x)/(2.0*(1.0+%c[%ld]*x)**3))",
          c, i, c, i+1, c, i+1, c, i+2, c, i+1, c, i+1, c, i+1);
      return stringBuffer;
    }
    case Isotherm::Type::Quadratic:
    {
      snprintf(stringBuffer, 1024, "%c[%ld]*(%c[%ld]*x+2.0*%c[%ld]*x**2)/(1.0+%c[%ld]*x+%c[%ld]*x**2)",
              c, i, c, i+1, c, i+2, c, i+1, c, i+2);
      return stringBuffer;
    }
    case Isotherm::Type::Temkin:
    {
      snprintf(stringBuffer, 1024, "%c[%ld]*(%c[%ld]*x/(1.0+%c[%ld]*x))+%c[%ld]*%c[%ld]*((%c[%ld]*x/(1.0+%c[%ld]*x))**2)*(%c[%ld]*x/(1.0+%c[%ld]*x)-1.0)", 
          c, i, c, i+1, c, i+1, c, i, c, i+2, c, i+1, c, i+1, c, i+1, c, i+1);
      return stringBuffer;
    }
    case Isotherm::Type::BingelWalton:
    {
      snprintf(stringBuffer, 1024, "%c[%ld]*(1.0-exp(-(%c[%ld]+%c[%ld])*x))/(1.0+(%c[%ld]/%c[%ld])*exp(-(%c[%ld]+%c[%ld])*x))", 
              c, i, c, i+1, c, i+2, c, i+2, c, i+1, c, i +1, c, i+2);
      return stringBuffer;
    }
    default:
      throw std::runtime_error("Error: unkown isotherm type");
  }
}
