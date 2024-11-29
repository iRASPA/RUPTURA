#include <cstddef>
#include <vector>
#include <tuple>

#include "component.h"
#include "inputreader.h"
#include "mixture_prediction.h"

#ifdef PYBUILD
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
namespace py = pybind11;
#endif  // PYBUILD

struct Breakthrough
{
  public:
    // Constructor overloading, Breakthrough object can be created using either fo the two following definitions,
    // compiler will decide constructor at runtime based on the number and types of input arguments
    Breakthrough(const InputReader &inputreader);
    Breakthrough(std::string _displayName, std::vector<Component> _components, size_t _carrierGasComponent,
                 size_t _numberOfGridPoints, size_t _printEvery, size_t _writeEvery, double _temperature,
                 double _p_total, double _columnVoidFraction, double _pressureGradient, double _particleDensity,
                 double _columnEntranceVelocity, double _columnLength, double _timeStep, size_t _numberOfTimeSteps,
                 bool _autoSteps, bool _pulse, double _pulseTime, const MixturePrediction _mixture);

    std::string repr() const;
    void initialize();
    void run();
    void computeStep(size_t step);

    void createPlotScript();
    void createMovieScripts();

#ifdef PYBUILD
    py::array_t<double> compute();
#endif  // PYBUILD

   private:
    const std::string displayName;
    const std::vector<Component> components;
    size_t carrierGasComponent{ 0 }; 
    size_t Ncomp;      // number of components
    size_t Ngrid;      // number of grid points

    size_t printEvery; // print time step to the screen every printEvery steps
    size_t writeEvery; // write data to files every writeEvery steps

    double T_gas;      // absolute temperature [K]
    double p_total;    // total pressure column [Pa]
    double dptdx;      // pressure gradient [N/m3]
    double epsilon;    // void-fraction of the column [-]
    double rho_p;      // particle density [kg/m3]
    double v_in;       // interstitial velocity at the begin of the column [m/s]
    
    double L;          // length of the column
    double dx;         // spacing in spatial direction
    double dt;         // timestep integration
    size_t Nsteps;     // total number of steps
    bool autoSteps;    // use automatic number of steps
    bool pulse;        // pulsed inlet condition for breakthrough
    double tpulse;     // pulse time

    
    MixturePrediction mixture;
    size_t maxIsothermTerms;
    std::pair<size_t, size_t> iastPerformance{ 0, 0 };
    
    // vector of size 'Ncomp'
    std::vector<double> prefactor;
    std::vector<double> Yi;        // ideal gas mol-fraction for each component
    std::vector<double> Xi;        // adsorbed mol-fraction for each component
    std::vector<double> Ni;        // number of molecules for each component
    
    // vector of size '(Ngrid + 1)'
    std::vector<double> V;         // interstitial gas velocity along the column
    std::vector<double> Vnew; 
    std::vector<double> Pt;        // total pressure along the column

    // Hassan Modification , vectors of size '(Ngrid + 1)'
    std::vector<double> T;       // Temperature along the column
    std::vector<double> Tnew;
    std::vector<double> DTdt;     // Temperature gradient along the column
    std::vector<double> DTdtnew;
    std::vector<double> P;        // Total pressure at every grid point
    std::vector<double> Pnew;
    std::vector<double> DPdt;     // derivative of P with respect to time
    std::vector<double> DPdtnew;


    // vector of size '(Ngrid + 1) * Ncomp', for each grid point, data per component (contiguous)
    std::vector<double> y;         // mole fraction at every grid point for each component
    std::vector<double> ynew;
    std::vector<double> Dydt;      // derivative of mole fraction at every grid point for each component
    std::vector<double> Dydtnew;
    std::vector<double> Q;         // volume-averaged adsorption amount at every grid point for each component
    std::vector<double> Qnew;
    std::vector<double> Qeq;       // equilibrium adsorption amount at every grid point for each component
    std::vector<double> Qeqnew;
    std::vector<double> Dqdt;      // derivative of Q with respect to time
    std::vector<double> Dqdtnew;
    std::vector<double> cachedP0;  // cached hypothetical pressure
    std::vector<double> cachedPsi; // cached reduced grand potential over the column
    
    // Properties and Parameters
    double K_z;                        // Thermal conductivity of gas [J/mol/K]
    double C_ps;                      // Heat capacity of adsorbent [J/kg/K]
    double C_pg;                       // Heat capacity of gas [J/mol/K]
    double C_pa;                       // Heat capacity of adsorbate [J/mol/K]
    double mu;                     // Viscoisty of gas [Pa.s]
    double r_p;                     // Radius of adsorbent particles [m]
    double Q_s0;                        // Scaling factor for adsorption loading [mol/kg]

    std::vector<double> MW;            // Molecular weights of components [kg/mol]
    
    // Isotherm Parameters
    std::vector <double> sat_q_b;
    std::vector <double> sat_q_d;
    std::vector <double> b_0;
    std::vector <double> d_0 ;
    std::vector<double> del_H;     // Delta Heat of components [mol/kg]
    double T_ref;
    
    enum class IntegrationScheme
    {
      SSP_RK = 0,
      Iterative = 1
    };

    void computeFirstDerivatives(std::vector<double> &dqdt,
                                 std::vector<double> &dpdt,
                                 std::vector<double> &dTdt,
                                 std::vector<double> &dydt,
                                 const std::vector<double> &q_eq,
                                 const std::vector<double> &q,
                                 const std::vector<double> &v,
                                 const std::vector<double> &p,
                                 const std::vector<double> &Temp,
                                 const std::vector<double> &y_vec);
                                //  const std::vector<double> &T_arg,
                                //  const std::vector<double> &y_arg);

    void computeEquilibriumLoadings();

    void computeVelocity();    
    
    void createMovieScriptColumnV();
    void createMovieScriptColumnPt();
    void createMovieScriptColumnQ();
    void createMovieScriptColumnQeq();
    void createMovieScriptColumnP();
    void createMovieScriptColumnDpdt();
    void createMovieScriptColumnDqdt();
    void createMovieScriptColumnPnormalized();
};
