#include "inputreader.h"

#include <cstddef>
#include <cctype>
#include <iostream>
#include <sstream>
#include <fstream>

bool caseInSensStringCompare(const std::string& str1, const std::string& str2)
{
  return str1.size() == str2.size() && std::equal(str1.begin(), str1.end(), str2.begin(), 
        [](int a, int b) {return std::tolower(a) == std::tolower(b); });
}

bool startsWith(const std::string &str, const std::string &prefix) {
    return str.size() >= prefix.size() && str.substr(0, prefix.size()) == prefix;
}

std::string trim(const std::string& s)
{
  auto start = s.begin();
  while (start != s.end() && std::isspace(*start)) 
  {
    start++;
  }

  auto end = s.end();
  do 
  {
    end--;
  } while (std::distance(start, end) > 0 && std::isspace(*end));

  return std::string(start, end + 1);
}

template<class T>
T parse(const std::string& arguments, [[maybe_unused]] const std::string& keyword, [[maybe_unused]] size_t lineNumber)
{
  T value;

  std::string str;
  std::istringstream ss(arguments);

  ss >> value;

  return value;
}

template<typename T>
std::vector<T> parseListOfSystemValues(const std::string& arguments, const std::string& keyword, size_t lineNumber)
{
  std::vector<T> list{};

  std::string str;
  std::istringstream ss(arguments);

  std::string errorString = "No values could be read for keyword '" + keyword + 
                            "' at line: " + std::to_string(lineNumber) + "\n";

  while (ss >> str)
  {
    if (trim(str).rfind("//", 0) == 0)
    {
      if (list.empty())
      {
        throw std::runtime_error(errorString);
      }
      return list;
    }
    T value;
    std::istringstream s(str);
    if (s >> value)
    {
      list.push_back(value);
    }
    else
    {
      if (list.empty())
      {
        throw std::runtime_error(errorString);
      }
      return list;
    }

  };

  if (list.empty())
  {
    throw std::runtime_error(errorString);
  }
  return list;
}

double parseDouble(const std::string& arguments, const std::string& keyword, size_t lineNumber)
{
  double value{};

  std::string str;
  std::istringstream ss(arguments);

  if (ss >> value)
  {
      return value;
  };

  std::string errorString = "Numbers could not be read for keyword '" + keyword + 
                            "' at line: " + std::to_string(lineNumber) + "\n";
  throw std::runtime_error(errorString);
}

int parseBoolean(const std::string& arguments, const std::string& keyword, size_t lineNumber)
{
  bool value{};

  std::istringstream ss(arguments);

  if (ss >> std::boolalpha >> value)
  {
      return value;
  };

  std::string str;
  std::istringstream ss2(arguments);
  if (ss2 >> str)
  {
    if (caseInSensStringCompare(str, "yes")) return true;
    if (caseInSensStringCompare(str, "no")) return false;
  };

  std::string errorString = "Booleands could not be read for keyword '" + keyword + 
                            "' at line: " + std::to_string(lineNumber) + "\n";
  throw std::runtime_error(errorString);
}



InputReader::InputReader(const std::string fileName):
  components()
{
  components.reserve(16);

  std::ifstream fileInput{ fileName };
  std::string errorOpeningFile = "Required input file '" + fileName + "' does not exist";
  if (!fileInput) throw std::runtime_error(errorOpeningFile);

  std::string line{};
  std::string keyword{};
  std::string arguments{};
  size_t lineNumber{ 0 };
  size_t numberOfComponents{ 0 };

  while (std::getline(fileInput, line))
  {
    lineNumber += 1;
    if (!line.empty())
    {
      std::istringstream iss(line);

      iss >> keyword;
      keyword = trim(keyword);
      std::getline(iss, arguments);
      arguments = trim(arguments);

      if (caseInSensStringCompare(keyword, "SimulationType"))
      {
        std::string str;
        std::istringstream ss(arguments);
        if (ss >> str)
        {
          if (caseInSensStringCompare(str, "Breakthrough")) 
          {
            simulationType = SimulationType::Breakthrough;
            continue;
          }
          if (caseInSensStringCompare(str, "MixturePrediction")) 
          {
            simulationType = SimulationType::MixturePrediction;
            continue;
          }
          if (caseInSensStringCompare(str, "Fitting")) 
          {
            simulationType = SimulationType::Fitting;
            continue;
          }
          if (caseInSensStringCompare(str, "Test")) 
          {
            simulationType = SimulationType::Test;
            continue;
          }
        };
      }
      if (caseInSensStringCompare(keyword, "MixturePredictionMethod"))
      {
        std::string str;
        std::istringstream ss(arguments);
        if (ss >> str)
        {
          if (caseInSensStringCompare(str, "IAST")) 
          {
            mixturePredictionMethod = 0;
            continue;
          }
          if (caseInSensStringCompare(str, "SIAST")) 
          {
            mixturePredictionMethod = 1;
            continue;
          }
          if (caseInSensStringCompare(str, "EI")) 
          {
            mixturePredictionMethod = 2;
            continue;
          }
          if (caseInSensStringCompare(str, "SEI")) 
          {
            mixturePredictionMethod = 3;
            continue;
          }
        };
      }
      if (caseInSensStringCompare(keyword, "IASTMethod"))
      {
        std::string str;
        std::istringstream ss(arguments);
        if (ss >> str)
        {
          if (caseInSensStringCompare(str, "FastIAS")) 
          {
            IASTMethod = 0;
            continue;
          }
          if (caseInSensStringCompare(str, "Bisection")) 
          {
            IASTMethod = 1;
            continue;
          }
        };
      }
      if (caseInSensStringCompare(keyword, "DisplayName"))
      {
        std::string str;
        std::istringstream ss(arguments);
        if (ss >> str)
        {
          this->displayName = str;
          continue;
        }
      }

      if (caseInSensStringCompare(keyword, "Temperature"))
      {
        double value = parseDouble(arguments, keyword, lineNumber);
        this->temperature = value;
        continue;
      }
      if (caseInSensStringCompare(keyword, "ColumnVoidFraction"))
      {
        double value = parseDouble(arguments, keyword, lineNumber);
        this->columnVoidFraction = value;
        continue;
      }
      if (caseInSensStringCompare(keyword, "ParticleDensity"))
      {
        double value = parseDouble(arguments, keyword, lineNumber);
        this->particleDensity = value;
        continue;
      }
      if (caseInSensStringCompare(keyword, "TotalPressure"))
      {
        double value = parseDouble(arguments, keyword, lineNumber);
        this->totalPressure = value;
        continue;
      }
      if (caseInSensStringCompare(keyword, "PressureStart"))
      {
        double value = parseDouble(arguments, keyword, lineNumber);
        this->pressureStart = value;
        continue;
      }
      if (caseInSensStringCompare(keyword, "PressureEnd"))
      {
        double value = parseDouble(arguments, keyword, lineNumber);
        this->pressureEnd = value;
        continue;
      }
      if (caseInSensStringCompare(keyword, "NumberOfPressurePoints"))
      {
        size_t value = parse<size_t>(arguments, keyword, lineNumber);
        this->numberOfPressurePoints = value;
        continue;
      }
      if (caseInSensStringCompare(keyword, "PressureScale"))
      {
        std::string str;
        std::istringstream ss(arguments);
        if (ss >> str)
        {
          if (caseInSensStringCompare(str, "Log")) 
          {
            pressureScale = 0;
            continue;
          }
          if (caseInSensStringCompare(str, "Linear")) 
          {
            pressureScale = 1;
            continue;
          }
          if (caseInSensStringCompare(str, "Normal")) 
          {
            pressureScale = 1;
            continue;
          }
        };
      }
      if (caseInSensStringCompare(keyword, "PressureGradient"))
      {
        double value = parseDouble(arguments, keyword, lineNumber);
        this->pressureGradient = value;
        continue;
      }
      if (caseInSensStringCompare(keyword, "ColumnEntranceVelocity"))
      {
        double value = parseDouble(arguments, keyword, lineNumber);
        this->columnEntranceVelocity = value;
        continue;
      }

      if (caseInSensStringCompare(keyword, "NumberOfTimeSteps"))
      {
        std::string str;
        std::istringstream ss(arguments);
        if (ss >> str)
        {
          if (caseInSensStringCompare(str, "auto"))
          {
            autoNumberOfTimeSteps = true;
          }
          else 
          {
            size_t value = parse<size_t>(arguments, keyword, lineNumber);
            this->numberOfTimeSteps = value;
            autoNumberOfTimeSteps = false;
          }
          continue;
        }
      }

      if (caseInSensStringCompare(keyword, "PulseBreakthrough"))
      {
        bool value = parseBoolean(arguments, keyword, lineNumber);
        this->pulseBreakthrough = value;
        continue;
      }

      if (caseInSensStringCompare(keyword, "PulseTime"))
      {
        double value = parseDouble(arguments, keyword, lineNumber);
        this->pulseTime = value;
        continue;
      }     
      if (caseInSensStringCompare(keyword, "TimeStep"))
      {
        double value = parseDouble(arguments, keyword, lineNumber);
        this->timeStep = value;
        continue;
      }
      if (caseInSensStringCompare(keyword, "PrintEvery"))
      {
        size_t value = parse<size_t>(arguments, keyword, lineNumber);
        this->printEvery = value;
        continue;
      }
      if (caseInSensStringCompare(keyword, "WriteEvery"))
      {
        size_t value = parse<size_t>(arguments, keyword, lineNumber);
        this->writeEvery = value;
        continue;
      }
      if (caseInSensStringCompare(keyword, "ColumnLength"))
      {
        double value = parseDouble(arguments, keyword, lineNumber);
        this->columnLength = value;
        continue;
      }
      if (caseInSensStringCompare(keyword, "NumberOfGridPoints"))
      {
        size_t value = parse<size_t>(arguments, keyword, lineNumber);
        this->numberOfGridPoints = value;
        continue;
      }

      if (caseInSensStringCompare(keyword, "ColumnPressure"))
      {
        size_t value = parse<size_t>(arguments, keyword, lineNumber);
        this->columnPressure = value;
        continue;
      }

      if (caseInSensStringCompare(keyword, "ColumnLoading"))
      {
        size_t value = parse<size_t>(arguments, keyword, lineNumber);
        this->columnLoading = value;
        continue;
      }

      if (caseInSensStringCompare(keyword, "ColumnError"))
      {
        size_t value = parse<size_t>(arguments, keyword, lineNumber);
        this->columnError = value;
        continue;
      }

      if (caseInSensStringCompare(keyword, std::string("Component")))
      {
        std::istringstream ss(arguments);
        std::cout << "arguments: " << arguments << std::endl;
        std::string c, moleculeNameKeyword,remainder,componentName;
        ss >> c >> moleculeNameKeyword >> componentName;
        std::getline(ss, remainder);

        components.push_back(Component(numberOfComponents, componentName));
        numberOfComponents += 1;
        continue;
      }
      if (caseInSensStringCompare(keyword, "FileName"))
      {
        std::string str;
        std::istringstream ss(arguments);
        if (ss >> str)
        {
          components[numberOfComponents - 1].filename = str;
          continue;
        }
      }
      if (caseInSensStringCompare(keyword, "CarrierGas"))
      {
        bool value = parseBoolean(arguments, keyword, lineNumber);
        components[numberOfComponents - 1].isCarrierGas = value;
        continue;
      }
      if (caseInSensStringCompare(keyword, "GasPhaseMolFraction"))
      {
        double value = parseDouble(arguments, keyword, lineNumber);
        components[numberOfComponents - 1].Yi0 = value;
        continue;
      }
      if (caseInSensStringCompare(keyword, "MassTransferCoefficient"))
      {
        double value = parseDouble(arguments, keyword, lineNumber);
        components[numberOfComponents - 1].Kl = value;
        continue;
      }
      if (caseInSensStringCompare(keyword, "AxialDispersionCoefficient"))
      {
        double value = parseDouble(arguments, keyword, lineNumber);
        components[numberOfComponents - 1].D = value;
        continue;
      }

      if (caseInSensStringCompare(keyword, "NumberOfIsothermSites"))
      {
        size_t value = parse<size_t>(arguments, keyword, lineNumber);
        components[numberOfComponents - 1].isotherm.numberOfSites = value;
        continue;
      }
      if (caseInSensStringCompare(keyword, "Langmuir"))
      {
        std::vector<double> values = parseListOfSystemValues<double>(arguments, keyword, lineNumber);
        if(values.size() < 2)
        {
          throw std::runtime_error("Error: Langmuir requires two parameters");
        }
        values.resize(2);
        Isotherm isotherm = Isotherm(Isotherm::Type::Langmuir, values, 2);
        components[numberOfComponents - 1].isotherm.add(isotherm);
        continue;
      }
      if (caseInSensStringCompare(keyword, "Anti-Langmuir"))
      {
        std::vector<double> values = parseListOfSystemValues<double>(arguments, keyword, lineNumber);
        if(values.size() < 2)
        {
          throw std::runtime_error("Error: Anti-Langmuir requires two parameters");
        }
        values.resize(2);
        Isotherm isotherm = Isotherm(Isotherm::Type::Anti_Langmuir, values, 2);
        components[numberOfComponents - 1].isotherm.add(isotherm);
        continue;
      }
      if (caseInSensStringCompare(keyword, "BET"))
      {
        std::vector<double> values = parseListOfSystemValues<double>(arguments, keyword, lineNumber);
        if(values.size() < 3)
        {
          throw std::runtime_error("Error: BET requires three parameters");
        }
        values.resize(3);
        Isotherm isotherm = Isotherm(Isotherm::Type::BET, values, 3);
        components[numberOfComponents - 1].isotherm.add(isotherm);
        continue;
      }
      if (caseInSensStringCompare(keyword, "Henry"))
      {
        std::vector<double> values = parseListOfSystemValues<double>(arguments, keyword, lineNumber);
        if(values.size() < 1)
        {
          throw std::runtime_error("Error: Henry requires one parameter");
        }
        values.resize(1);
        Isotherm isotherm = Isotherm(Isotherm::Type::Henry, values, 1);
        components[numberOfComponents - 1].isotherm.add(isotherm);
        continue;
      }
      if (caseInSensStringCompare(keyword, "Freundlich"))
      {
        std::vector<double> values = parseListOfSystemValues<double>(arguments, keyword, lineNumber);
        if(values.size() < 2)
        {
          throw std::runtime_error("Error: Freundlich requires two parameters");
        }
        values.resize(2);
        Isotherm isotherm = Isotherm(Isotherm::Type::Freundlich, values, 2);
        components[numberOfComponents - 1].isotherm.add(isotherm);
        continue;
      }
      if (caseInSensStringCompare(keyword, "Sips"))
      {
        std::vector<double> values = parseListOfSystemValues<double>(arguments, keyword, lineNumber);
        if(values.size() < 3)
        {
          throw std::runtime_error("Error: Sips requires three parameters");
        }
        values.resize(3);
        Isotherm isotherm = Isotherm(Isotherm::Type::Sips, values, 3);
        components[numberOfComponents - 1].isotherm.add(isotherm);
        continue;
      }
      if (caseInSensStringCompare(keyword, "Langmuir-Freundlich"))
      {
        std::vector<double> values = parseListOfSystemValues<double>(arguments, keyword, lineNumber);
        if(values.size() < 3)
        {
          throw std::runtime_error("Error: Langmuir-Freundlich requires three parameters");
        }
        values.resize(3);
        Isotherm isotherm = Isotherm(Isotherm::Type::Langmuir_Freundlich, values, 3);
        components[numberOfComponents - 1].isotherm.add(isotherm);
        continue;
      }
      if (caseInSensStringCompare(keyword, "Redlich-Peterson"))
      {
        std::vector<double> values = parseListOfSystemValues<double>(arguments, keyword, lineNumber);
        if(values.size() < 3)
        {
          throw std::runtime_error("Error: Redlich-Peterson requires three parameters");
        }
        values.resize(3);
        Isotherm isotherm = Isotherm(Isotherm::Type::Redlich_Peterson, values, 3);
        components[numberOfComponents - 1].isotherm.add(isotherm);
        continue;
      }
      if (caseInSensStringCompare(keyword, "Toth"))
      {
        std::vector<double> values = parseListOfSystemValues<double>(arguments, keyword, lineNumber);
        if(values.size() < 3)
        {
          throw std::runtime_error("Error: Toth requires three parameters");
        }
        values.resize(3);
        Isotherm isotherm = Isotherm(Isotherm::Type::Toth, values, 3);
        components[numberOfComponents - 1].isotherm.add(isotherm);
        continue;
      }
      if (caseInSensStringCompare(keyword, "Unilan"))
      {
        std::vector<double> values = parseListOfSystemValues<double>(arguments, keyword, lineNumber);
        if(values.size() < 3)
        {
          throw std::runtime_error("Error: Unilan requires three parameters");
        }
        values.resize(3);
        Isotherm isotherm = Isotherm(Isotherm::Type::Unilan, values, 3);
        components[numberOfComponents - 1].isotherm.add(isotherm);
        continue;
      }
      if (caseInSensStringCompare(keyword, "O'Brian&Myers"))
      {
        std::vector<double> values = parseListOfSystemValues<double>(arguments, keyword, lineNumber);
        if(values.size() < 3)
        {
          throw std::runtime_error("Error: O'Brien&Myers requires three parameters");
        }
        values.resize(3);
        Isotherm isotherm = Isotherm(Isotherm::Type::OBrien_Myers, values, 3);
        components[numberOfComponents - 1].isotherm.add(isotherm);
        continue;
      }
      if (caseInSensStringCompare(keyword, "Quadratic"))
      {
        std::vector<double> values = parseListOfSystemValues<double>(arguments, keyword, lineNumber);
        if(values.size() < 3)
        {
          throw std::runtime_error("Error: Quadratic requires three parameters");
        }
        values.resize(3);
        Isotherm isotherm = Isotherm(Isotherm::Type::Quadratic, values, 3);
        components[numberOfComponents - 1].isotherm.add(isotherm);
        continue;
      }
      if (caseInSensStringCompare(keyword, "Temkin"))
      {
        std::vector<double> values = parseListOfSystemValues<double>(arguments, keyword, lineNumber);
        if(values.size() < 3)
        {
          throw std::runtime_error("Error: Temkin requires three parameters");
        }
        values.resize(3);
        Isotherm isotherm = Isotherm(Isotherm::Type::Temkin, values, 3);
        components[numberOfComponents - 1].isotherm.add(isotherm);
        continue;
      }
      if (caseInSensStringCompare(keyword, "Bingel&Walton"))
      {
        std::vector<double> values = parseListOfSystemValues<double>(arguments, keyword, lineNumber);
        if(values.size() < 3)
        {
          throw std::runtime_error("Error: Bingel&Walton requires three parameters");
        }
        values.resize(3);
        Isotherm isotherm = Isotherm(Isotherm::Type::BingelWalton, values, 3);
        components[numberOfComponents - 1].isotherm.add(isotherm);
        continue;
      }

      if(!(startsWith(keyword, "//") || startsWith(keyword, "#")))
      {
        std::cout << "Error: unknown keyword (" << keyword << ") with arguments (" << arguments << ")" << std::endl;
        exit(0);
      }
    }
  }

  // normalize gas-phase mol-fractions to unity
  if(simulationType != SimulationType::Fitting)
  {
    double sum = 0.0;
    for(size_t j = 0; j < components.size(); ++j)
    {
      sum += components[j].Yi0;
    }
    if(std::abs(sum-1.0)>1e-15)
    {
      std::cout << "Normalizing: Gas-phase molfractions did not sum exactly to unity!\n\n";
      for(size_t j = 0; j < components.size(); ++j)
      {
        components[j].Yi0 /= sum;
      }
    }
  }

  numberOfCarrierGases = 0;
  carrierGasComponent = 0;
  for(size_t j = 0; j < components.size(); ++j)
  {
    if(components[j].isCarrierGas) 
    {
      carrierGasComponent = j;
      std::vector<double> values{1.0, 0.0};
      Isotherm isotherm = Isotherm(Isotherm::Type::Langmuir, values, 2);
      components[carrierGasComponent].isotherm.add(isotherm);
      components[carrierGasComponent].isotherm.numberOfSites = 1;

      ++numberOfCarrierGases;
    }
  }

  if((mixturePredictionMethod == 2) || (mixturePredictionMethod == 3))
  {
    for(size_t i = 0; i < components.size(); ++i)
    {
      for(size_t j = 0; j < components[i].isotherm.numberOfSites; ++j)
      {
        if( components[i].isotherm.sites[j].type != Isotherm::Type::Langmuir)
        {
          throw std::runtime_error("Error: Explicit mixture prediction must use single Langmuir isotherms");
        }
      }
    }
  }


  maxIsothermTerms = 0;
  if(!components.empty())
  {
    std::vector<Component>::iterator maxIsothermTermsIterator = std::max_element(components.begin(), components.end(),
              [] (Component& lhs, Component& rhs) {
                  return lhs.isotherm.numberOfSites < rhs.isotherm.numberOfSites;
              });
    maxIsothermTerms = maxIsothermTermsIterator->isotherm.numberOfSites;
  }

  if(simulationType == SimulationType::Breakthrough)
  {
    if(numberOfCarrierGases == 0)
    {
      throw std::runtime_error("Error: no carrier gas component present");
    }
    if(numberOfCarrierGases > 1)
    {
      throw std::runtime_error("Error: multiple carrier gas component present (there can be only one)");
    }
    if(temperature < 0.0)
    {
      throw std::runtime_error("Error: temperature not set (Use e.g.: 'Temperature 300'");
    }
    if(columnVoidFraction < 0.0)
    {
      throw std::runtime_error("Error: void-fraction of the colum not set (Use e.g.: 'ColumnVoidFraction 0.4'");
    }
    if(particleDensity < 0.0)
    {
      throw std::runtime_error("Error: particle density not set (Use e.g.: 'ParticleDensity 1408.2'");
    }
    if(totalPressure < 0.0)
    {
      throw std::runtime_error("Error: total pressure bot set (Use e.g.: 'TotalPressure 1e5'");
    }
    if(columnEntranceVelocity < 0.0)
    {
      throw std::runtime_error("Error: column entrance velocity not set (Use e.g.: 'columnEntranceVelocity 300'");
    }

    if((numberOfTimeSteps == 0) && (!autoNumberOfTimeSteps))
    {
      throw std::runtime_error("Error: number of time steps not set (Use e.g.: 'NumberOfTimeSteps 5000000'");
    }
    if(numberOfGridPoints == 0)
    {
      throw std::runtime_error("Error: number of grid points not set (Use e.g.: 'NumberOfGridPoints 50'");
    }
    if(columnLength < 0)
    {
      throw std::runtime_error("Error: column length not set (Use e.g.: 'ColumnLength 0.3'");
    }
  }
}
