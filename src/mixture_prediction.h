#pragma once

#include <vector>
#include <tuple>

#include "inputreader.h"
#include "component.h"

#ifdef PYBUILD
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
namespace py = pybind11;
#endif  // P

class MixturePrediction
{
  public:
    enum class PredictionMethod
    {
      IAST = 0,
      SIAST = 1,
      EI = 2,
      SEI = 3
    };

    enum class IASTMethod
    {
      FastIAST = 0,
      NestedLoopBisection = 1
    };

    MixturePrediction(const InputReader &inputreader);
    MixturePrediction(std::string _displayName, std::vector<Component> _components, size_t _numberOfCarrierGases,
                      size_t _carrierGasComponent, double _temperature, double _pressureStart, double _pressureEnd,
                      size_t _numberOfPressurePoints, size_t _pressureScale, size_t _predictionMethod,
                      size_t _iastMethod);

    std::string repr() const;
    void sortComponents();
    void run();
    std::vector<double> initPressures();
    void createPureComponentsPlotScript();
    void createMixturePlotScript();
    void createMixtureAdsorbedMolFractionPlotScript();
    void createPlotScript();

    // Yi  = gas phase mol-fraction
    // P   = total pressure
    // Xi  = adsorbed phase mol-fraction
    // Ni  = number of adsorbed molecules of component i
    std::pair<size_t, size_t> predictMixture(const std::vector<double> &Yi,
                                             const double &P,
                                             std::vector<double> &Xi,
                                             std::vector<double> &Ni,
                                             double *cachedP0,
                                             double *cachedPsi);

    // keep this non private for breakthrough
    size_t maxIsothermTerms;

#ifdef PYBUILD
    py::array_t<double> compute();
#endif  // PYBUILD

   private:
    std::string displayName;
    const std::vector<Component> components;
    std::vector<Component> sortedComponents;
    const size_t Ncomp;
    const size_t Nsorted;
    size_t numberOfCarrierGases;
    size_t carrierGasComponent;
    PredictionMethod predictionMethod;
    IASTMethod iastMethod;
    std::vector<std::vector<Component>> segregatedSortedComponents;

    std::vector<double> alpha1;
    std::vector<double> alpha2;
    std::vector<double> alpha_prod;
    std::vector<double> x;

    std::vector<double> pstar;
    std::vector<double> psi;
    std::vector<double> G;
    std::vector<double> delta;
    std::vector<double> Phi;

    enum class PressureScale
    {
      Log = 0,
      Normal = 1
    };
    double temperature{ 300.0 };
    double pressureStart{ 1e3 };
    double pressureEnd{ 1e8 };
    size_t numberOfPressurePoints{ 100 };
    PressureScale pressureScale{ PressureScale::Log };

    std::pair<size_t, size_t> computeFastIAST(const std::vector<double> &Yi,
                                      const double &P,
                                      std::vector<double> &Xi,
                                      std::vector<double> &Ni,
                                      double *cachedP0,
                                      double *cachedPsi);
    std::pair<size_t, size_t> computeFastSIAST(const std::vector<double> &Yi,
                                      const double &P,
                                      std::vector<double> &Xi,
                                      std::vector<double> &Ni,
                                      double *cachedP0,
                                      double *cachedPsi);
    std::pair<size_t, size_t> computeFastSIAST(size_t term,
                                      const std::vector<double> &Yi,
                                      const double &P,
                                      std::vector<double> &Xi,
                                      std::vector<double> &Ni,
                                      double *cachedP0,
                                      double *cachedPsi);

    std::pair<size_t, size_t> computeIASTNestedLoopBisection(const std::vector<double> &Yi,
                                      const double &P,
                                      std::vector<double> &Xi,
                                      std::vector<double> &Ni,
                                      double *cachedP0,
                                      double *cachedPsi);
    std::pair<size_t, size_t> computeSIASTNestedLoopBisection(const std::vector<double> &Yi,
                                      const double &P,
                                      std::vector<double> &Xi,
                                      std::vector<double> &Ni,
                                      double *cachedP0,
                                      double *cachedPsi);
    std::pair<size_t, size_t> computeSIASTNestedLoopBisection(size_t term,
                                      const std::vector<double> &Yi,
                                      const double &P,
                                      std::vector<double> &Xi,
                                      std::vector<double> &Ni,
                                      double *cachedP0,
                                      double *cachedPsi);
    std::pair<size_t, size_t> computeExplicitIsotherm(const std::vector<double> &Yi,
                                                      const double &P,
                                                      std::vector<double> &Xi,
                                                      std::vector<double> &Ni);
    std::pair<size_t, size_t> computeSegratedExplicitIsotherm(const std::vector<double> &Yi,
                                                      const double &P,
                                                      std::vector<double> &Xi,
                                                      std::vector<double> &Ni);
    std::pair<size_t, size_t> computeSegratedExplicitIsotherm(size_t site, const std::vector<double> &Yi,
                                                      const double &P,
                                                      std::vector<double> &Xi,
                                                      std::vector<double> &Ni);

    void printErrorStatus(double psi, double sum, double P, const std::vector<double> Yi, double cachedP0[]);
};

