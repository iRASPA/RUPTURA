#pragma once

#include <tuple>
#include <array>
#include <vector>
#include <string>
#include <unordered_map>

#include "inputreader.h"
#include "component.h"
#include "multi_site_isotherm.h"

#ifdef PYBUILD
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
namespace py = pybind11;
#endif  // PYBUILD

struct Fitting
{
  struct DNA
  {
    DNA(std::string g, MultiSiteIsotherm p, double f):
        genotype(g),
        phenotype(p),
        fitness(f),
        hash(std::hash<std::string>{}(g))
    {
    }
    DNA() noexcept = default;

    std::string genotype;
    MultiSiteIsotherm phenotype;
    double fitness;
    size_t hash;
  };

  enum class PressureScale
  {
    Log = 0,
    Normal = 1
  };

  Fitting(const InputReader &inputreader);
  Fitting(std::string _displayName, std::vector<Component> _components, size_t _pressureScale);

  void readData(size_t ID);
  void run();
  void createPlotScripts(const DNA &citizen, size_t ID);
  void createPlotScript();

  DNA newCitizen(size_t ID);
  void updateCitizen(DNA &citizen);
  double fitness(const MultiSiteIsotherm &phenotype);
  double RCorrelation(const MultiSiteIsotherm &phenotype);
  size_t biodiversity(const std::vector<DNA> &citizens);
  void nuclearDisaster(size_t ID);
  void elitism();
  void mutate(DNA &Mutant);
  void crossover(size_t ID, size_t s1,size_t s2, size_t i1, size_t i2, size_t j1, size_t j2);
  void chooseRandomly(size_t kk1,size_t kk2,size_t jj1,size_t jj2, size_t &ii1, size_t &ii2);
  void mate(size_t ID);
  void sortByFitness();
  void writeCitizen(size_t citizen, size_t id, size_t step, size_t variety, size_t fullfilledCondition);
  DNA fit(size_t ID);
  const DNA simplex(DNA citizen, double scale);

  size_t Ncomp;
  std::string displayName;
  std::vector<Component> components;
  std::vector<std::string> filename;
  size_t columnPressure{ 0 };
  size_t columnLoading{ 1 };
  size_t columnError{ 2 };
  double maximumLoading{ 0.0 };
  PressureScale pressureScale{ PressureScale::Log };

  std::vector<std::pair<double, double>> rawData;

  bool fittingFlag{ false };
  bool physicalConstrainsFlag{ false };
  bool seedFlag{ false };
  bool pressureRangeFlag{ false };
  bool refittingFlag{ false };
  std::pair<double, double> pressureRange;
  std::pair<double, double> logPressureRange;

  size_t GA_Size;            // population size
  double GA_MutationRate;    // mutation rate
  double GA_EliteRate;       // elitists population rate
  double GA_MotleyCrowdRate; // pirates population rate
  double GA_DisasterRate;
  size_t GA_Elitists;        // number of elitists
  size_t GA_Motleists;       // number of pirates

  std::vector<DNA> popAlpha;
  std::vector<DNA> popBeta;
  std::vector<DNA> &parents;
  std::vector<DNA> &children;

#ifdef PYBUILD
  void selectData(size_t ID, std::vector<std::vector<std::pair<double, double>>> data);
  std::vector<double> compute(std::vector<std::vector<std::pair<double, double>>> data);
  py::array_t<double> evaluate(std::vector<double> pressure);
#endif
};
