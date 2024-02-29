#include "fitting.h"
#include "special_functions.h"
#include "random_numbers.h"

#include <fstream>
#include <sstream>
#include <iterator>
#include <map>
#include <algorithm>
#include <exception>
#include <cmath>
#include <cstdlib>
#include <bitset>
#include <cstring>
#include <climits>
#include <unordered_set>
#if __cplusplus >= 201703L && __has_include(<filesystem>)
  #include <filesystem>
#elif __cplusplus >= 201703L && __has_include(<experimental/filesystem>)
  #include <experimental/filesystem>
#else
  #include <sys/stat.h>
#endif

#ifdef PYBUILD
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
namespace py = pybind11;
#endif  // PYBUILDa

Fitting::Fitting(const InputReader &inputreader)
    : Ncomp(inputreader.components.size()),
      displayName(inputreader.displayName),
      components(inputreader.components),
      filename(Ncomp),
      columnPressure(inputreader.columnPressure - 1),
      columnLoading(inputreader.columnLoading - 1),
      columnError(inputreader.columnError - 1),
      pressureScale(PressureScale(inputreader.pressureScale)),
      GA_Size(static_cast<size_t>(std::pow(2.0, 12.0))),
      GA_MutationRate(1.0 / 3.0),
      GA_EliteRate(0.15),
      GA_MotleyCrowdRate(0.25),
      GA_DisasterRate(0.001),
      GA_Elitists(static_cast<size_t>(static_cast<double>(GA_Size) * GA_EliteRate)),
      GA_Motleists(static_cast<size_t>(static_cast<double>(GA_Size) * (1.0 - GA_MotleyCrowdRate))),
      popAlpha(static_cast<size_t>(std::pow(2.0, 12.0))),
      popBeta(static_cast<size_t>(std::pow(2.0, 12.0))),
      parents(popAlpha),
      children(popBeta)
{
  for(size_t i = 0 ; i < Ncomp; ++i)
  {
    filename[i] = inputreader.components[i].filename;
  }
}

Fitting::Fitting(std::string _displayName, std::vector<Component> _components, size_t _pressureScale)
    : Ncomp(_components.size()),
      displayName(_displayName),
      components(_components),
      pressureScale(PressureScale(_pressureScale)),
      GA_Size(static_cast<size_t>(std::pow(2.0, 12.0))),
      GA_MutationRate(1.0 / 3.0),
      GA_EliteRate(0.15),
      GA_MotleyCrowdRate(0.25),
      GA_DisasterRate(0.001),
      GA_Elitists(static_cast<size_t>(static_cast<double>(GA_Size) * GA_EliteRate)),
      GA_Motleists(static_cast<size_t>(static_cast<double>(GA_Size) * (1.0 - GA_MotleyCrowdRate))),
      popAlpha(static_cast<size_t>(std::pow(2.0, 12.0))),
      popBeta(static_cast<size_t>(std::pow(2.0, 12.0))),
      parents(popAlpha),
      children(popBeta)
{
}

void Fitting::readData(size_t ID)
{
  std::ifstream fileInput{ filename[ID] };
  std::string errorOpeningFile = "File '" + filename[ID] + "' exists, but error opening file";
  if (!fileInput) throw std::runtime_error(errorOpeningFile);

  std::cout << "Reading: " << filename[ID] << "\n";

  std::string line{};

  maximumLoading = 0.0;
  rawData.clear();
  while (std::getline(fileInput, line))
  {
    std::string trimmedLine = trim(line);
    if(!startsWith(trimmedLine, "#"))
    {
      if (!line.empty())
      {
        std::istringstream iss(line);

        std::vector<std::string> results((std::istream_iterator<std::string>(iss)),
                                 std::istream_iterator<std::string>());
        if(columnPressure < results.size() &&
           columnLoading < results.size())
        {
          double pressure;
          double loading;
          std::istringstream s(results[columnPressure]);
          s >> pressure;
          std::istringstream t(results[columnLoading]);
          t >> loading;
          if(loading > maximumLoading)
          {
            maximumLoading = loading;
          }
          rawData.push_back(std::make_pair(pressure, loading));
        }
      }
    }
  }

  if(rawData.empty())
  {
    throw std::runtime_error("Error: no pressure points found");
  }

  // sort the pressures
  std::sort(rawData.begin(), rawData.end());

  pressureRange = std::make_pair(rawData.front().first, rawData.back().first);
  logPressureRange = std::make_pair(std::log(pressureRange.first), std::log(pressureRange.second));

  std::cout << "Found " << rawData.size() << " data points\n";
  for(const std::pair<double, double> &data : rawData)
  {
    std::cout << data.first << " " << data.second << std::endl;
  }
  std::cout << "\n";
  std::cout << "Lowest pressure: " << pressureRange.first << std::endl;
  std::cout << "Highest pressure: " << pressureRange.second << std::endl;
  std::cout << "Log lowest pressure: " << logPressureRange.first << std::endl;
  std::cout << "Log highest pressure: " << logPressureRange.second << std::endl;

  for(size_t i = 0; i < Ncomp; ++i)
  {
    std::cout << "Number of isotherm parameters: " << components[i].isotherm.numberOfParameters << std::endl;
    std::cout << components[i].isotherm.repr();
  }
}

void Fitting::run()
{
  std::cout << "STARTING FITTING\n";
  for(size_t i = 0; i < Ncomp; ++i)
  {
    readData(i);
    const DNA bestCitizen = fit(i);
    const DNA optimizedBestCitizen = simplex(bestCitizen, 1.0);
    std::cout << optimizedBestCitizen.phenotype.repr();
    createPlotScripts(optimizedBestCitizen, i);
  }
  createPlotScript();
}

// create a new citizen in the Ensemble
Fitting::DNA Fitting::newCitizen(size_t ID)
{
  DNA citizen;

  citizen.phenotype = components[ID].isotherm.randomized(maximumLoading);

  citizen.genotype.clear();
  citizen.genotype.reserve((sizeof(double) * CHAR_BIT) *
                  citizen.phenotype.numberOfParameters);
  for(size_t i = 0; i < citizen.phenotype.numberOfParameters; ++i)
  {
    // convert from double to bitset
    uint64_t p;
    std::memcpy(&p, &citizen.phenotype.parameters(i), sizeof(double));
    std::bitset<sizeof(double) * CHAR_BIT> bitset(p);

    // add the bit-string to the genotype representation
    citizen.genotype += bitset.to_string();
  }

  citizen.hash = std::hash<MultiSiteIsotherm>{}(citizen.phenotype);
  citizen.fitness = fitness(citizen.phenotype);

  return citizen;
}

void Fitting::updateCitizen(DNA &citizen)
{
  citizen.fitness = fitness(citizen.phenotype);
}

inline bool my_isnan(double val) 
{
  union { double f; uint64_t x; } u = { val };
  return (u.x << 1) > (0x7ff0000000000000u << 1);
}


double Fitting::fitness(const MultiSiteIsotherm &phenotype)
// For evaluating isotherm goodness-of-fit:
// Residual Root Mean Square Error (RMSE)
{
  double fitnessValue = phenotype.fitness();
  size_t m = rawData.size();              // number of observations
  size_t p = phenotype.numberOfParameters; // number of adjustable parameters
  for(std::pair<double, double> dataPoint: rawData)
  {
    double pressure = dataPoint.first;
    double loading = dataPoint.second;
    double difference = loading - phenotype.value(pressure);
    // double weight = 1.0/(1.0+loading);
    double weight = 1.0;
    fitnessValue += weight * difference * difference;
  }
  fitnessValue = sqrt(fitnessValue / static_cast<double>(m - p));

  if(my_isnan(fitnessValue)) fitnessValue = 99999999.999999;
  if(fitnessValue==0.0000000000) fitnessValue = 99999999.999999;

  return fitnessValue;
}

double Fitting::RCorrelation(const MultiSiteIsotherm &phenotype)
{
  double RCorrelationValue = phenotype.fitness();
  size_t m = rawData.size();
  double loading_avg_o = 0.0;
  double loading_avg_e = 0.0;
  double tmp1 = 0.0;
  double tmp2 = 0.0;
  double tmp3 = 0.0;

  for(std::pair<double, double> dataPoint: rawData)
  {
    double pressure = dataPoint.first;
    double loading = dataPoint.second;
    loading_avg_o += loading / static_cast<double>(m);
    loading_avg_e += phenotype.value(pressure) / static_cast<double>(m);
  }

  for(std::pair<double, double> dataPoint: rawData)
  {
    double pressure = dataPoint.first;
    double loading = dataPoint.second; 
    tmp1 += (loading-loading_avg_o)*(phenotype.value(pressure)-loading_avg_e);
    tmp2 += (loading-loading_avg_o)*(loading-loading_avg_o);
    tmp3 += (phenotype.value(pressure)-loading_avg_e)*(phenotype.value(pressure)-loading_avg_e);
  }
  RCorrelationValue = tmp1/sqrt(tmp2*tmp3);
  
  return RCorrelationValue;
}

size_t Fitting::biodiversity(const std::vector<DNA> &citizens)
{
  std::map<size_t, size_t> counts;
  for(const DNA &dna: citizens) 
  {
    if(counts.find(dna.hash) != counts.end())
    {
      ++counts[dna.hash];
    }
    else 
    {
      counts[dna.hash] = 1;
    }
  }
  size_t biodiversity = 0;
  for(const std::pair<size_t, size_t> value: counts)
  {
    if(value.second > 1)
    {
      biodiversity += value.second;
    }
  }

  return biodiversity;
}

void Fitting::nuclearDisaster(size_t ID)
{
  for(size_t i = 1; i < children.size(); ++i)
  {
    children[i] = newCitizen(ID);
  }
}

void Fitting::elitism()
{
  std::copy(parents.begin(), parents.begin() + static_cast<std::vector<DNA>::difference_type>(GA_Elitists), children.begin());
}

void Fitting::mutate(DNA &mutant)
{
  mutant.genotype.clear();
  mutant.genotype.reserve((sizeof(double) * CHAR_BIT) * mutant.phenotype.numberOfParameters);
  for(size_t i = 0; i < mutant.phenotype.numberOfParameters; ++i)
  {
    // convert from double to bitset
    uint64_t p;
    std::memcpy(&p, &mutant.phenotype.parameters(i), sizeof(double));
    std::bitset<sizeof(double) * CHAR_BIT> bitset(p);

    // mutation: randomly flip bit
    bitset.flip(std::size_t((sizeof(double) * CHAR_BIT) * RandomNumber::Uniform()));

    // convert from bitset to double
    p = bitset.to_ullong();
    std::memcpy(&mutant.phenotype.parameters(i), &p, sizeof(double));

    // add the bit-string to the genotype representation
    mutant.genotype += bitset.to_string();
  }

  // calculate the hah-value from the entire bit-string
  mutant.hash = std::hash<MultiSiteIsotherm>{}(mutant.phenotype);
}

// [s1:s2] range of children
// [i1:i2] range of parent1
// [j1:j2] range of parent2
// One-point crossover
//----------------------------------------------
//  parent1    parent2                children
//    *          *                      *
//  00|000000  11|111111        ->    00|111111
//----------------------------------------------
// Two-point
//----------------------------------------------
//  parent1     parent2               children
//    *    *      *    *                *    *
//  00|0000|00  11|1111|11       ->   00|1111|00
//----------------------------------------------
void Fitting::crossover(size_t ID, size_t s1,size_t s2, size_t i1, size_t i2, size_t j1, size_t j2)
{
  size_t k1,k2;
  double  tmp1;
  for(size_t i = s1; i < s2; ++i)
  {
    chooseRandomly(i1, i2, j1, j2, k1, k2);
    tmp1 = RandomNumber::Uniform();
    // choose between single cross-over using bit-strings or random parameter-swap
    if(tmp1 < 0.490)
      // One-point crossover:
      // --------------------
    {
      // remove the extreme values 0 and 32*Npar - 1 (they are not valid for crossover)
      size_t bitStringSize = (sizeof(double) * CHAR_BIT) * components[ID].isotherm.numberOfParameters;
      size_t spos = RandomNumber::Integer(1, bitStringSize - 2);
      children[i].genotype = parents[k1].genotype.substr(0, spos) +
                             parents[k2].genotype.substr(spos, bitStringSize - spos);

      // convert the bit-strings to doubles
      for(size_t j = 0; j < children[i].phenotype.numberOfParameters; ++j)
      {
        size_t pos = j * (sizeof(double) * CHAR_BIT);
        size_t size = sizeof(double) * CHAR_BIT;
        std::bitset<sizeof(double) * CHAR_BIT> bitset(children[i].genotype, pos, size);
        uint64_t p = bitset.to_ullong();
        std::memcpy(&children[i].phenotype.parameters(j), &p, sizeof(double));
      }
    }
    else if ( tmp1 < 0.499)
      // Two-point crossover:
      // --------------------
      {
        size_t bitStringSize = (sizeof(double) * CHAR_BIT) * components[ID].isotherm.numberOfParameters;
        size_t spos1 = RandomNumber::Integer(1, bitStringSize - 3);
        size_t spos2 = RandomNumber::Integer(spos1, bitStringSize - 2);
        children[i].genotype = parents[k1].genotype.substr(0, spos1) +
                               parents[k2].genotype.substr(spos1, spos2 - spos1) +
                               parents[k1].genotype.substr(spos2, bitStringSize - spos2);
        // convert the bit-strings to doubles
        for (size_t j = 0; j < children[i].phenotype.numberOfParameters; ++j)
        {
          size_t pos = j * (sizeof(double) * CHAR_BIT);
          size_t size = sizeof(double) * CHAR_BIT;
          std::bitset<sizeof(double) * CHAR_BIT> bitset(children[i].genotype, pos, size);
          uint64_t p = bitset.to_ullong();
          std::memcpy(&children[i].phenotype.parameters(j), &p, sizeof(double));
        }
      }
      else if (tmp1 < 0.500)
      {
        // Uniform crossover:
        // ------------------
        size_t bitStringSize = (sizeof(double) * CHAR_BIT) * components[ID].isotherm.numberOfParameters;
        size_t rolling_k = k1;
        for (size_t j = 0; j < bitStringSize; j++)
        {
          if (RandomNumber::Uniform() < 0.25)
          {
            if (rolling_k == k1)
            {
              rolling_k = k2;
            }
            else
            {
              rolling_k = k1;
            }
          }
          children[i].genotype.substr(j, 1) = parents[rolling_k].genotype.substr(j, 1);
        }
        // convert the bit-strings to doubles
        for (size_t j = 0; j < children[i].phenotype.numberOfParameters; ++j)
        {
          size_t pos = j * (sizeof(double) * CHAR_BIT);
          size_t size = sizeof(double) * CHAR_BIT;
          std::bitset<sizeof(double) * CHAR_BIT> bitset(children[i].genotype, pos, size);
          uint64_t p = bitset.to_ullong();
          std::memcpy(&children[i].phenotype.parameters(j), &p, sizeof(double));
        }
      }
      else
      {
        children[i].genotype.clear();
        for (size_t j = 0; j < children[i].phenotype.numberOfParameters; ++j)
        {
          // randomly choose whether the parameter comes from parent k1 or k2
          if (RandomNumber::Uniform() < 0.5)
          {
            children[i].phenotype.parameters(j) = parents[k1].phenotype.parameters(j);
          }
          else
          {
            children[i].phenotype.parameters(j) = parents[k2].phenotype.parameters(j);
          }

          // convert from double to bitString
          uint64_t p;
          std::memcpy(&p, &children[i].phenotype.parameters(j), sizeof(double));
          std::bitset<sizeof(double) * CHAR_BIT> bitset(p);
          children[i].genotype += bitset.to_string();
        }
      }

      children[i].hash = std::hash<MultiSiteIsotherm>{}(children[i].phenotype);
  }
}

void Fitting::chooseRandomly(size_t kk1, size_t kk2, size_t jj1, size_t jj2, size_t &ii1, size_t &ii2)
{
  ii1 = RandomNumber::Integer(kk1, kk2);
  ii2 = RandomNumber::Integer(jj1, jj2);
  while (ii1 == ii2)
  {
    ii2 = RandomNumber::Integer(jj1, jj2);
  };
}

void Fitting::mate(size_t ID)
{
  // retain the first 25% of the children
  elitism();

  // mates from GA_Elitists to (GA_Size - GA_Elitists)
  crossover(ID, GA_Elitists, GA_Elitists + static_cast<size_t>(static_cast<double>(GA_Motleists) * 0.5), 0, GA_Elitists,
            0, GA_Elitists);
  crossover(ID, GA_Elitists + static_cast<size_t>(static_cast<double>(GA_Motleists) * 0.5) + 1, GA_Size - GA_Elitists,
            0, GA_Elitists, GA_Elitists, GA_Size - 1);

  // mutation from GA_Elitists to (GA_Size - GA_Elitists) with "GA_MutationRate" probability
  for (size_t i = GA_Elitists; i < GA_Size - GA_Elitists; ++i)
  {
    if (RandomNumber::Uniform() < GA_MutationRate)
    {
      mutate(children[i]);
    }
    updateCitizen(children[i]);
  }

  // replace the last GA_Elitists (the worst) of the children by new children
  for (size_t i = GA_Size - GA_Elitists; i < GA_Size; ++i)
  {
    children[i] = newCitizen(ID);
  }

  // replace the last (GA_Size - 1) children by new children
  if (RandomNumber::Uniform() < GA_DisasterRate)
  {
    nuclearDisaster(ID);
  }
}

bool DNA_Fitness_Sorter(Fitting::DNA const &lhs, Fitting::DNA const &rhs) { return lhs.fitness < rhs.fitness; }

void Fitting::sortByFitness() { std::sort(parents.begin(), parents.end(), &DNA_Fitness_Sorter); }

void Fitting::writeCitizen(size_t citizen, size_t id, size_t step, size_t variety, size_t fullfilledCondition)
{
  char info[256];
  if (fullfilledCondition > 0)
  {
    snprintf(info, 256,
             "mol: %2ld  step: %5ld  Fitness: %10.6lf R^2: %10.6lf Similarity: %5ld/%-5ld Finishing: %3ld/%-3d\n", id,
             step, parents[citizen].fitness, pow(RCorrelation(parents[citizen].phenotype), 2), variety, GA_Size,
             fullfilledCondition, 100);
  }
  else
  {
    snprintf(info, 256, "mol: %2ld  step: %5ld  Fitness: %10.6lf R^2: %10.6lf Similarity: %5ld/%-5ld\n", id, step,
             parents[citizen].fitness, pow(RCorrelation(parents[citizen].phenotype), 2), variety, GA_Size);
  }
  std::cout << info;
  std::cout << "number of parameters: " << parents[citizen].phenotype.numberOfParameters << std::endl;
  for (size_t i = 0; i < parents[citizen].phenotype.numberOfParameters; ++i)
  {
    std::cout << "      genotype: " << parents[citizen].genotype.substr(64 * i, 64)
              << " parameter: " << parents[citizen].phenotype.parameters(i) << "\n";
  }
  std::cout << std::endl;
}

Fitting::DNA Fitting::fit(size_t ID)
{
  size_t optimisationStep{0};
  const size_t maxOptimisationStep{1000};
  size_t fullFilledConditionStep{0};
  const size_t maxFullfilledConditionStep{100};
  const double minimumFitness{5.0e-1};
  double tempFitnessValue{999.0};
  size_t tempVarietyValue{0};
  const double toleranceEqualFitness{1e-3};
  const size_t minstep{10};

  fullFilledConditionStep = 0;
  optimisationStep = 0;

  for (size_t i = 0; i < popAlpha.size(); ++i)
  {
    popAlpha[i] = newCitizen(ID);
    popBeta[i] = newCitizen(ID);
  }

  parents = popAlpha;
  children = popBeta;

  sortByFitness();

  // print Initial (and unsorted) population
  writeCitizen(0, ID, 0, 0.0, 0);

  if (refittingFlag)
  {
    std::cout << "Refitting activated\n";
    for (size_t citizen = 0; citizen < 2; ++citizen)
    {
      parents[citizen].genotype.clear();
      parents[citizen].genotype.reserve((sizeof(double) * CHAR_BIT) * parents[citizen].phenotype.numberOfParameters);
      for (size_t i = 0; i < parents[citizen].phenotype.numberOfParameters; ++i)
      {
        // convert from double to bitset
        uint64_t p;
        std::memcpy(&p, &parents[citizen].phenotype.parameters(i), sizeof(double));
        std::bitset<sizeof(double) * CHAR_BIT> bitset(p);

        // add the bit-string to the genotype representation
        parents[citizen].genotype += bitset.to_string();
      }
      parents[citizen].hash = std::hash<MultiSiteIsotherm>{}(parents[citizen].phenotype);
      updateCitizen(parents[citizen]);
    }
    std::copy(parents.begin(), parents.end(), children.begin());

    mate(ID);

    std::swap(parents, children);

    sortByFitness();
  }

  components[ID].isotherm = parents[0].phenotype;

  std::cout << "Starting Genetic Algorithm optimization\n";

  bool continueCondition = true;
  do
  {
    sortByFitness();

    tempVarietyValue = biodiversity(children);

    writeCitizen(0, ID, optimisationStep, tempVarietyValue, fullFilledConditionStep);

    if (optimisationStep >= minstep && parents[0].fitness <= minimumFitness &&
        std::abs(parents[0].fitness - tempFitnessValue) <= toleranceEqualFitness)
    {
      fullFilledConditionStep += 1;
    }
    else
    {
      fullFilledConditionStep = 0;
    }

    if (optimisationStep >= maxOptimisationStep || fullFilledConditionStep >= maxFullfilledConditionStep)
    {
      continueCondition = false;
    }

    // pairing
    mate(ID);

    std::swap(parents, children);

    // take the best fitness value:
    tempFitnessValue = parents[0].fitness;

    optimisationStep += 1;

  } while (continueCondition);

  writeCitizen(0, ID, optimisationStep, tempVarietyValue, fullFilledConditionStep);

  return parents[0];
}

// The Nelder-Mead method uses a simplex (a hyper-tetrahedron of n+1 vertices in n dimensions)
// Advantage: it does not use derivates, works well, can tolerate some noise
// Disadvantage: it is not garanteed to converge
// The steps of the method are:
// 1) Sort: according to the fitness
// 2) Reflect: get rid of the worst point, replace it by something better
// 3) Extend: if better then extend it even further
// 4) Contract: if not better then contract it
// 5) Shrink: if still not better we shrink towards the best performing point
// 6) Check convergence
const Fitting::DNA Fitting::simplex(DNA citizen, double scale)
{
  size_t n = citizen.phenotype.numberOfParameters;
  std::vector<std::vector<double>> v(n + 1, std::vector<double>(n));  // holds vertices of simplex
  std::vector<double> f(n + 1);                                       // value of function at each vertex
  std::vector<double> vr(n);                                          // reflection - coordinates
  std::vector<double> ve(n);                                          // expansion - coordinates
  std::vector<double> vc(n);                                          // contraction - coordinates
  std::vector<double> vm(n);                                          // centroid - coordinates
  std::vector<double> vtmp(n);                                        // temporary array passed to FUNC
  size_t vs;                                                          // vertex with the smallest value
  size_t vh;                                                          // vertex with next smallest value
  size_t vg;                                                          // vertex with largest value
  double fr;                                                          // value of function at reflection point
  double fe;                                                          // value of function at expansion point
  double fc;                                                          // value of function at contraction point
  size_t iprint{0};
  size_t MAX_IT = 1000000;
  double EPSILON = 1.0e-4;
  double ALPHA = 1.0;
  double BETA = 0.5;
  double GAMMA = 2.0;

  std::cout << "\nMinimising the cost function using the Nelder-Mead SIMPLEX method:\n\n";

  for (size_t i = 0; i < n; ++i)
  {
    v[0][i] = citizen.phenotype.parameters(i);
  }

  // values used to create initial simplex
  double pn = scale * (std::sqrt(static_cast<double>(n) + 1.0) - 1.0 + static_cast<double>(n)) /
              (static_cast<double>(n) * sqrt(2.0));
  double qn = scale * (std::sqrt(static_cast<double>(n) + 1.0) - 1.0) / (static_cast<double>(n) * sqrt(2.0));
  for (size_t i = 1; i <= n; ++i)
  {
    for (size_t j = 0; j < n; ++j)
    {
      if (i - 1 == j)
      {
        v[i][j] = pn + citizen.phenotype.parameters(j);
      }
      else
      {
        v[i][j] = qn + citizen.phenotype.parameters(j);
      }
    }
  }

  for (size_t i = 0; i <= n; ++i)
  {
    for (size_t j = 0; j < n; ++j)
    {
      citizen.phenotype.parameters(j) = v[i][j];
    }
    f[i] = fitness(citizen.phenotype);
  }

  // print out the initial simplex
  // print out the initial function values
  // find the index of the smallest value for printing
  if (iprint == 0)
  {
    vs = 0;
    for (size_t j = 0; j <= n; ++j)
    {
      if (f[j] < f[vs])
      {
        vs = j;
      }
    }

    std::cout << "Initial Values from genetic algorithm:\n";

    for (size_t j = 0; j < n; ++j)
    {
      std::cout << v[vs][j] << " ";
    }
    std::cout << "Fit: " << f[vs] << "\n\n";
  }

  for (size_t itr = 1; itr <= MAX_IT; ++itr)
  {
    // Step 1: Sort
    // ====================================================================
    std::vector<size_t> sortIndexes = sort_indexes(f);
    vs = sortIndexes[0];      // index of smallest
    vg = sortIndexes[n];      // index of largest
    vh = sortIndexes[n - 1];  // index of second largest

    // calculate the center point of every point except for the worst one
    for (size_t j = 0; j < n; ++j)
    {
      double cent = 0.0;
      for (size_t i = 0; i <= n; ++i)
      {
        if (i != vg)
        {
          cent = cent + v[i][j];
        }
      }
      vm[j] = cent / static_cast<double>(n);
    }

    // Step 2: Reflect vg to new vertex vr
    // ====================================================================
    for (size_t j = 0; j < n; ++j)
    {
      vr[j] = (1.0 + ALPHA) * vm[j] - ALPHA * v[vg][j];
      citizen.phenotype.parameters(j) = vr[j];
    }
    fr = fitness(citizen.phenotype);

    if ((fr <= f[vh]) && (fr > f[vs]))
    {
      for (size_t j = 0; j < n; ++j)
      {
        v[vg][j] = vr[j];
      }
      f[vg] = fr;
    }

    // Step 3: Extend a step further in this direction
    // ====================================================================
    if (fr <= f[vs])
    {
      for (size_t j = 0; j < n; ++j)
      {
        ve[j] = GAMMA * vr[j] + (1.0 - GAMMA) * vm[j];
        citizen.phenotype.parameters(j) = ve[j];
      }
      fe = fitness(citizen.phenotype);

      // by making fe < fr as opposed to fe < f(vs), Rosenbrocks function
      // takes 62 iterations as opposed to 64.

      if (fe < fr)
      {
        for (size_t j = 0; j < n; ++j)
        {
          v[vg][j] = ve[j];
        }
        f[vg] = fe;
      }
      else
      {
        for (size_t j = 0; j < n; ++j)
        {
          v[vg][j] = vr[j];
        }
        f[vg] = fr;
      }
    }

    if (fr > f[vh])
    {
      // Step 4: Contraction
      // ====================================================================
      for (size_t j = 0; j < n; ++j)
      {
        vc[j] = BETA * v[vg][j] + (1.0 - BETA) * vm[j];
        citizen.phenotype.parameters(j) = vc[j];
      }
      fc = fitness(citizen.phenotype);
      if (fc < f[vg])
      {
        for (size_t j = 0; j < n; ++j)
        {
          v[vg][j] = vc[j];
        }
        f[vg] = fc;
      }
      else
      {
        // Step 4: Shrink
        // ====================================================================
        // at this point the contraction is not successful,
        // we must halve the distance from vs to all the
        // vertices of the simplex and then continue.
        for (size_t row = 0; row <= n; ++row)
        {
          if (row != vs)
          {
            for (size_t j = 0; j < n; ++j)
            {
              v[row][j] = v[vs][j] + 0.5 * (v[row][j] - v[vs][j]);
            }
          }
        }
        for (size_t m = 0; m < n; ++m)
        {
          vtmp[m] = v[vg][m];
          citizen.phenotype.parameters(m) = vtmp[m];
        }
        f[vg] = fitness(citizen.phenotype);

        for (size_t m = 0; m < n; ++m)
        {
          vtmp[m] = v[vh][m];
          citizen.phenotype.parameters(m) = vtmp[m];
        }
        f[vh] = fitness(citizen.phenotype);
      }
    }

    // Step 6: Test for convergence
    // ====================================================================
    double fsum = 0.0;
    for (size_t j = 0; j <= n; ++j)
    {
      fsum = fsum + f[j];
    }
    double favg = fsum / (static_cast<double>(n) + 1.0);

    if (favg < EPSILON || itr == MAX_IT)
    {
      // print out the value at each iteration

      if (itr != MAX_IT)
      {
        std::cout << "Nelder-Mead has converged: " << favg << " < " << EPSILON << "\n\n";
      }
      else
      {
        std::cout << "Reached maximum number of steps: " << itr << " = " << MAX_IT << "\n\n";
      }

      // find the index of the smallest value
      sortIndexes = sort_indexes(f);
      vs = sortIndexes[0];  // index of smallest
      for (size_t m = 0; m < n; ++m)
      {
        citizen.phenotype.parameters(m) = v[vs][m];
      }
      double min = fitness(citizen.phenotype);

      std::cout << "Final Values: " << std::endl;
      for (size_t j = 0; j < n; ++j)
      {
        std::cout << v[vs][j] << " ";
      }
      std::cout << "Fit: " << min << " R2: " << pow(RCorrelation(citizen.phenotype), 2) << "\n\n";

      return citizen;
    }
  }

  return citizen;
}

void Fitting::createPlotScripts(const DNA &citizen, size_t ID)
{
  std::string plotFileName = "plot_fit_component_" + std::to_string(ID) + "_" + components[ID].name;
  std::ofstream stream(plotFileName);

  stream << "set encoding utf8\n";
  stream << "set xlabel 'Pressure, {/Helvetica-Italic p} / [Pa]' font \"Helvetica,18\"\n";
  stream << "set ylabel 'Absolute loading, {/Helvetica-Italic q} / [mol/kg]' offset 0.0,0 font \"Helvetica,18\"\n";
  stream << "set bmargin 4\n";
  stream << "set yrange[0:]\n";
  if (pressureScale == PressureScale::Log)
  {
    stream << "set log x\n";
  }

  stream << "set key  right bottom vertical samplen 2.5 height 0.5 spacing 1.5 font 'Helvetica, 10'\n";
  stream << "set key title '" << components[ID].name << "'\n";

  stream << "set output 'isotherms_fit_" << components[ID].name << ".pdf'\n";
  stream << "set term pdf color solid\n";

  // colorscheme from book 'gnuplot in action', listing 12.7
  stream << "set linetype 1 pt 5 ps 0.5 lw 2 lc rgb '0xee0000'\n";
  stream << "set linetype 2 pt 7 ps 0.5 lw 2 lc rgb '0x008b00'\n";
  stream << "set linetype 3 pt 9 ps 0.5 lw 2 lc rgb '0x0000cd'\n";
  stream << "set linetype 4 pt 11 ps 0.5 lw 2 lc rgb '0xff3fb3'\n";
  stream << "set linetype 5 pt 13 ps 0.5 lw 2 lc rgb '0x00cdcd'\n";
  stream << "set linetype 6 pt 15 ps 0.5 lw 2 lc rgb '0xcd9b1d'\n";
  stream << "set linetype 7 pt  4 ps 0.5 lw 2 lc rgb '0x8968ed'\n";
  stream << "set linetype 8 pt  6 ps 0.5 lw 2 lc rgb '0x8b8b83'\n";
  stream << "set linetype 9 pt  8 ps 0.5 lw 2 lc rgb '0x00bb00'\n";
  stream << "set linetype 10 pt 10 ps 0.5 lw 2 lc rgb '0x1e90ff'\n";
  stream << "set linetype 11 pt 12 ps 0.5 lw 2 lc rgb '0x8b2500'\n";
  stream << "set linetype 12 pt 14 ps 0.5 lw 2 lc rgb '0x000000'\n";

  stream << "array s[" << components[ID].isotherm.numberOfParameters << "]\n";
  for (size_t i = 0; i < components[ID].isotherm.numberOfParameters; ++i)
  {
    stream << "s[" << i + 1 << "]=" << components[ID].isotherm.parameters(i) << "\n";
  }
  stream << "array p[" << citizen.phenotype.numberOfParameters << "]\n";
  for (size_t i = 0; i < citizen.phenotype.numberOfParameters; ++i)
  {
    stream << "p[" << i + 1 << "]=" << citizen.phenotype.parameters(i) << "\n";
  }
  stream << "plot \\\n"
         << components[ID].isotherm.gnuplotFunctionString('s') << " title 'start f(x)' with li dt 2 lw 2,\\\n"
         << citizen.phenotype.gnuplotFunctionString('p') << " title 'fit f(x)' with li lw 2,\\\n";
  stream << "'" << filename[ID] << "' us " << columnPressure + 1 << ":" << columnLoading + 1
         << " title 'raw data' with po pt 5 ps 0.5\n";
}

void Fitting::createPlotScript()
{
  #if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
    std::ofstream stream_graphs("make_graphs.bat");
    stream_graphs << "set PATH=%PATH%;C:\\Program Files\\gnuplot\\bin;C:\\Program Files\\ffmpeg-master-latest-win64-gpl\\bin;C:\\Program Files\\ffmpeg\\bin\n";
    for(size_t i = 0; i < Ncomp; ++i)
    {
      std::string plotFileName = "plot_fit_component_" + std::to_string(i) + "_" + components[i].name;
      stream_graphs << "gnuplot.exe " << plotFileName << "\n";
    }
  #else
    std::ofstream stream_graphs("make_graphs");
    stream_graphs << "#!/bin/sh\n";
    stream_graphs << "export LC_ALL='en_US.UTF-8'";
    for(size_t i = 0; i < Ncomp; ++i)
    {
      std::string plotFileName = "plot_fit_component_" + std::to_string(i) + "_" + components[i].name;
      stream_graphs << "gnuplot " << plotFileName << "\n";
    }
  #endif

  #if (__cplusplus >= 201703L)
    std::filesystem::path path{"make_graphs"};
    std::filesystem::permissions(path, std::filesystem::perms::owner_exec, std::filesystem::perm_options::add);
  #else
    chmod("make_graphs", S_IRWXU);
  #endif

}

#ifdef PYBUILD

void Fitting::selectData(size_t ID, std::vector<std::vector<std::pair<double, double>>> data)
{
  rawData = data[ID];

  // get pressure range
  pressureRange = std::make_pair(rawData.front().first, rawData.back().first);
  logPressureRange = std::make_pair(std::log(pressureRange.first), std::log(pressureRange.second));

  maximumLoading = 0.0;
  for (const auto &pair : rawData)
  {
    if (pair.second > maximumLoading)
    {
      maximumLoading = pair.second;
    }
  }
}

std::vector<double> Fitting::compute(std::vector<std::vector<std::pair<double, double>>> data)
{
  std::vector<double> output;
  std::cout << "STARTING FITTING\n";

  for (size_t ID = 0; ID < Ncomp; ++ID)
  {
    selectData(ID, data);

    const DNA bestCitizen = fit(ID);
    DNA optimizedBestCitizen = simplex(bestCitizen, 1.0);

    // save optimized params to component, insert in output array
    for (size_t j = 0; j < optimizedBestCitizen.phenotype.numberOfParameters; j++)
    {
      components[ID].isotherm.setParameters(j, optimizedBestCitizen.phenotype.parameters(j));
      output.push_back(optimizedBestCitizen.phenotype.parameters(j));
    }
  }
  for (size_t ID = 0; ID < Ncomp; ++ID)
  {
    std::cout << components[ID].repr() << "\n";
  }
  return output;
}

py::array_t<double> Fitting::evaluate(std::vector<double> pressure)
{
  // initialize numpy array
  size_t Npress = pressure.size();
  std::array<size_t, 2> shape{{Npress, Ncomp}};
  py::array_t<double> output(shape);
  double *data = output.mutable_data();

  // add datapoints
  for (size_t i = 0; i < Npress; i++)
  {
    for (size_t j = 0; j < Ncomp; j++)
    {
      data[i * Ncomp + j] = components[j].isotherm.value(pressure[i]);
    }
  }
  return output;
}

#endif  // PYBUILD
