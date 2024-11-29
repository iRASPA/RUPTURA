#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include <limits>
#include <algorithm>
#include <numeric>
#include <sstream>
#if __cplusplus >= 201703L && __has_include(<filesystem>)
  #include <filesystem>
#elif __cplusplus >= 201703L && __has_include(<experimental/filesystem>)
  #include <experimental/filesystem>
#else
  #include <sys/stat.h>
#endif

#include "breakthrough.h"

const double R=8.31446261815324;

inline double maxVectorDifference(const std::vector<double> &v, const std::vector<double> &w)
{
  if(v.empty() || w.empty()) return 0.0;
  if(v.size() != w.size()) throw std::runtime_error("Error: unequal vector size\n");

  double max = std::abs(v[0] - w[0]);
  for(size_t i = 1; i < v.size(); ++i)
  {
    double temp = std::abs(v[i] - w[i]);
    if(temp > max) max = temp;
  }
  return max;
}


// allow std::pairs to be added
template <typename T,typename U>                                                   
std::pair<T,U> operator+(const std::pair<T,U> & l,const std::pair<T,U> & r) {   
    return {l.first+r.first,l.second+r.second};
}
template <typename T, typename U>
std::pair<T,U> &operator+=(std::pair<T,U> & l, const std::pair<T,U> & r) {   
    l.first += r.first;
    l.second += r.second;
    return l;
}
// Constructor declaration with initializer list
Breakthrough::Breakthrough(const InputReader &inputReader):
    displayName(inputReader.displayName),
    components(inputReader.components),
    carrierGasComponent(inputReader.carrierGasComponent),
    Ncomp(components.size()),
    Ngrid(inputReader.numberOfGridPoints),
    printEvery(inputReader.printEvery),
    writeEvery(inputReader.writeEvery),
    T_gas(inputReader.temperature),
    p_total(inputReader.totalPressure),
    dptdx(inputReader.pressureGradient),
    epsilon(inputReader.columnVoidFraction),
    rho_p(inputReader.particleDensity),
    v_in(inputReader.columnEntranceVelocity),
    L(inputReader.columnLength),
    dx(L / static_cast<double>(Ngrid)),
    dt(inputReader.timeStep),
    Nsteps(inputReader.numberOfTimeSteps),
    autoSteps(inputReader.autoNumberOfTimeSteps),
    pulse(inputReader.pulseBreakthrough),
    tpulse(inputReader.pulseTime),
    mixture(inputReader),
    maxIsothermTerms(inputReader.maxIsothermTerms),
    prefactor(Ncomp),
    Yi(Ncomp),
    Xi(Ncomp),
    Ni(Ncomp),
    V(Ngrid+1),
    Vnew(Ngrid+1),
    Pt(Ngrid+1),
    T(Ngrid + 1),
    Tnew(Ngrid + 1),
    DTdt(Ngrid + 1),
    DTdtnew(Ngrid + 1),
    P(Ngrid + 1),
    Pnew(Ngrid + 1),
    DPdt(Ngrid + 1),
    DPdtnew(Ngrid + 1),
    y((Ngrid + 1) * Ncomp),
    ynew((Ngrid + 1) * Ncomp),
    Dydt((Ngrid + 1) * Ncomp),
    Dydtnew((Ngrid + 1) * Ncomp),
    Q((Ngrid + 1) * Ncomp),
    Qnew((Ngrid + 1) * Ncomp),
    Qeq((Ngrid + 1) * Ncomp),
    Qeqnew((Ngrid + 1) * Ncomp),
    Dqdt((Ngrid + 1) * Ncomp),
    Dqdtnew((Ngrid + 1) * Ncomp),
    cachedP0((Ngrid + 1) * Ncomp * maxIsothermTerms),
    cachedPsi((Ngrid + 1) * maxIsothermTerms) 
{
}
Breakthrough::Breakthrough(std::string _displayName, std::vector<Component> _components, size_t _carrierGasComponent,
                           size_t _numberOfGridPoints, size_t _printEvery, size_t _writeEvery, double _temperature,
                           double _p_total, double _columnVoidFraction, double _pressureGradient,
                           double _particleDensity, double _columnEntranceVelocity, double _columnLength,
                           double _timeStep, size_t _numberOfTimeSteps, bool _autoSteps, bool _pulse, double _pulseTime,
                           const MixturePrediction _mixture)
    : displayName(_displayName),
      components(normalize_molfracs(_components)),
      carrierGasComponent(_carrierGasComponent),
      Ncomp(_components.size()),
      Ngrid(_numberOfGridPoints),
      printEvery(_printEvery),
      writeEvery(_writeEvery),
      T_gas(_temperature),
      p_total(_p_total),
      dptdx(_pressureGradient),
      epsilon(_columnVoidFraction),
      rho_p(_particleDensity),
      v_in(_columnEntranceVelocity),
      L(_columnLength),
      dx(L / static_cast<double>(Ngrid)),
      dt(_timeStep),
      Nsteps(_numberOfTimeSteps),
      autoSteps(_autoSteps),
      pulse(_pulse),
      tpulse(_pulseTime),
      mixture(_mixture),
      maxIsothermTerms(mixture.maxIsothermTerms),
      prefactor(Ncomp),
      Yi(Ncomp),
      Xi(Ncomp),
      Ni(Ncomp),
      V(Ngrid + 1),
      Vnew(Ngrid + 1),
      Pt(Ngrid + 1),
      T(Ngrid + 1),
      Tnew(Ngrid + 1),
      DTdt(Ngrid + 1),
      DTdtnew(Ngrid + 1),
      P(Ngrid + 1),
      Pnew(Ngrid + 1),
      DPdt(Ngrid + 1),
      DPdtnew(Ngrid + 1),
      y((Ngrid + 1) * Ncomp),
      ynew((Ngrid + 1) * Ncomp),
      Dydt((Ngrid + 1) * Ncomp),
      Dydtnew((Ngrid + 1) * Ncomp),
      Q((Ngrid + 1) * Ncomp),
      Qnew((Ngrid + 1) * Ncomp),
      Qeq((Ngrid + 1) * Ncomp),
      Qeqnew((Ngrid + 1) * Ncomp),
      Dqdt((Ngrid + 1) * Ncomp),
      Dqdtnew((Ngrid + 1) * Ncomp),
      cachedP0((Ngrid + 1) * Ncomp * maxIsothermTerms),
      cachedPsi((Ngrid + 1) * maxIsothermTerms)
{
  // normally ran in main.cpp, now run by default
  initialize();
}

void Breakthrough::initialize()
{
  // Hassan Properties and Parameters
  K_z = 0.09;                        
  C_ps = 750.0;                      
  C_pg = 35.8;                       
  C_pa = 35.8;                       
  mu = 1.13e-05;                     
  r_p = 5.0e-03;                     
  Q_s0 = 3.0;
  MW = {0.004, 0.044, 0.028};

  sat_q_b = {0.00, 2.74, 3.12};
  sat_q_d = {0.00, 3.11, 0.00};
  b_0 = {0.00, 2.00e-3, 4.87e-6};
  d_0 = {0.00, 2.70e-5, 0.00};
  del_H = {0.0, -38.87e3, -20.79e3};
  T_ref = 273.00;

  // Initialize nondimensional time step and grid
  dt = dt * v_in/L;
  dx = dx / L;                     
  // precomputed factor for mass transfer
  // Hassan: Directly computed at runtime for every iteration because of variable temperature
  // for(size_t j = 0; j < Ncomp; ++j)
  // {
    // prefactor[j] = R * T * ((1.0 - epsilon) / epsilon) * rho_p * components[j].Kl;
  // }

  // set P and Q to zero
  // std::fill(P.begin(), P.end(), 0.0);
  std::fill(Q.begin(), Q.end(), 0.0);

  // Hassan: Initialize the temperature and mole fractions
  std::fill(T.begin(), T.end(), T_gas/T_gas);   // Equal to T_gas
  std::fill(y.begin(), y.end(), 0.0);

  // initial pressure along the column
  std::vector<double> pt_init(Ngrid + 1);

  // set the initial total pressure along the column assuming the pressure gradient is constant
  for(size_t i = 0; i < Ngrid + 1; ++i)
  {
    pt_init[i] = (p_total + dptdx * static_cast<double>(i) * dx*L) / p_total;
  }

  // initialize the interstitial gas velocity in the column
  for(size_t i = 0; i < Ngrid + 1; ++i)
  {
    // V[i] = v_in * p_total / pt_init[i];
    // V[i] = (v_in * 1 / pt_init[i]) / v_in;
    V[i] = 0.0;
  }
  V[0] = v_in/v_in;

  // set the molefraction of the carrier gas equal to 1.0, as column is initially filled with carrier gas only.
  // for the column except for the entrance (i=0)
  for(size_t i = 1; i < Ngrid + 1; ++i)
  {
    y[i * Ncomp + carrierGasComponent] = pt_init[i]/pt_init[i];
  }

  // at the column entrance, the mol-fractions of the components in the gas phase are fixed
  // the partial pressures of the components at the entrance are the mol-fractions times the 
  // total pressure
  for(size_t j = 0; j < Ncomp; ++j)
  {
    // P[0 * Ncomp + j] = p_total * components[j].Yi0;
    y[0 * Ncomp + j] = components[j].Yi0;
  }

  // at the entrance: mol-fractions Yi are the gas-phase mol-fractions
  // for the column: the initial mol-fraction of the carrier-gas is 1, and 0 for the other components
  //
  // the K of the carrier gas is chosen as zero 
  // so Qeq is zero for all components in the column after the entrance
  // only the values for Yi at the entrance are effected by adsorption
  for(size_t i = 0; i < Ngrid + 1; ++i)
  {
    // double sum = 0.0;
    // for(size_t j = 0; j < Ncomp; ++j)
    // {
    //   Yi[j] = std::max(P[i * Ncomp + j] / pt_init[i], 0.0);
    //   sum += Yi[j];
    // }
    for(size_t j = 0; j < Ncomp; ++j)
    {
      Yi[j] = std::max(y[i * Ncomp + j], 0.0);
    }

    // iastPerformance += mixture.predictMixture(Yi, pt_init[i], Xi, Ni, 
    //     &cachedP0[i * Ncomp * maxIsothermTerms], &cachedPsi[i * maxIsothermTerms]);

  
    iastPerformance += mixture.predictMixture(Yi, pt_init[i]*p_total, Xi, Ni, 
        &cachedP0[i * Ncomp * maxIsothermTerms], &cachedPsi[i * maxIsothermTerms]);

    for(size_t j = 0; j < Ncomp; ++j)
    {
      Qeq[i * Ncomp + j] = Ni[j];
    }
  }

  for(size_t i = 0; i < Ngrid + 1; ++i)
  {
    // Pt[i] = 0.0;
    // for(size_t j = 0; j < Ncomp; ++j)
    // {
    //   Pt[i] += std::max(0.0, P[i * Ncomp + j]);
    // }
    // Initial pressure profile is same as pt_init
    P[i] += std::max(pt_init[i], 0.0);
  }

  // check the MW vector size, it should not be less or more than Ncomp
  size_t length = MW.size();
  if (length != Ncomp)
  {
    throw std::runtime_error("Error: Mismatch in MW vector and Ncomp\n");
  }

  length = del_H.size();
  if (length != Ncomp)
  {
    throw std::runtime_error("Error: Mismatch in MW vector and Ncomp\n");
  }
}

void Breakthrough::run()
{
  // create the output files
  std::vector<std::ofstream> streams;
  for (size_t i = 0; i < Ncomp; i++)
  {
    std::string fileName = "component_" + std::to_string(i) + "_" + components[i].name + ".data";
    streams.emplace_back(std::ofstream{fileName});
  }

  std::ofstream movieStream("column.data");

  size_t column_nr = 1;
  movieStream << "# column " << column_nr++ << ": z  (column position)\n";
  movieStream << "# column " << column_nr++ << ": V  (velocity)\n";
  movieStream << "# column " << column_nr++ << ": P (total pressure)\n";
  movieStream << "# column " << column_nr++ << ": T  (Gas Temperature)\n";
  movieStream << "# column " << column_nr++ << ": DPdt (derivative P with t)\n";
  movieStream << "# column " << column_nr++ << ": DTdt (derivative T with t)\n";

  for (size_t j = 0; j < Ncomp; ++j)
  {
    movieStream << "# column " << column_nr++ << ": component " << j << " Q     (loading) \n";
    movieStream << "# column " << column_nr++ << ": component " << j << " Qeq   (equilibrium loading)\n";
    movieStream << "# column " << column_nr++ << ": component " << j << " y     (mole fraction)\n";
    movieStream << "# column " << column_nr++ << ": component " << j << " ynorm (normalized mole fraction)\n";
    movieStream << "# column " << column_nr++ << ": component " << j << " Dydt  (derivative y with t)\n";
    movieStream << "# column " << column_nr++ << ": component " << j << " Dqdt  (derivative Q with tn\n";
  }

  for (size_t step = 0; (step < Nsteps || autoSteps); ++step)
  {
    // compute new step
    computeStep(step);

    double t = static_cast<double>(step) * dt;

    if (step % writeEvery == 0)
    {
      // write breakthrough output to files
      // column 1: dimensionless time
      // column 2: time [minutes]
      // column 3: normalized mole fraction
      for (size_t j = 0; j < Ncomp; ++j)
      {
        // streams[j] << t * v_in / L << " " << t / 60.0 << " "
        //            << P[Ngrid * Ncomp + j] / ((p_total + dptdx * L) * components[j].Yi0) << std::endl;
        // streams[j] << t * v_in / L << " " << t / 60.0 << " "
        //            << y[Ngrid * Ncomp + j] / components[j].Yi0 << std::endl;
        streams[j] << t * L / v_in << " " << t * L / v_in / 60.0 << " "
                   << y[Ngrid * Ncomp + j] / components[j].Yi0 << std::endl;
      }

      for (size_t i = 0; i < Ngrid + 1; ++i)
      {
        movieStream << static_cast<double>(i) * dx*L << " ";
        movieStream << V[i] * v_in << " ";
        movieStream << P[i] * p_total<< " ";
        movieStream << T[i] * T_gas<< " ";
        movieStream << DPdt[i] * p_total*v_in/L<< " ";
        movieStream << DTdt[i] * T_gas*v_in/L<< " ";

        for (size_t j = 0; j < Ncomp; ++j)
        {
          movieStream << Q[i * Ncomp + j] * Q_s0 << " " << Qeq[i * Ncomp + j] << " " << y[i * Ncomp + j] << " "
                      << y[i * Ncomp + j] / (components[j].Yi0) << " " << Dydt[i * Ncomp + j] * v_in/L<< " "
                      << Dqdt[i * Ncomp + j] * Q_s0*v_in/L << " ";
        }
        movieStream << "\n";
      }
      movieStream << "\n\n";
    }

    if (step % printEvery == 0)
    {
      std::cout << "Timestep " + std::to_string(step) + ", time: " + std::to_string(t) + " [s]\n";
      std::cout << "    Average number of mixture-prediction steps: " +
                       std::to_string(static_cast<double>(iastPerformance.first) /
                                      static_cast<double>(iastPerformance.second))
                << std::endl;
    }
  }

  std::cout << "Final timestep " + std::to_string(Nsteps) +
                   ", time: " + std::to_string(dt * static_cast<double>(Nsteps)) + " [s]\n";
}

void Breakthrough::computeStep(size_t step)
{
  double t = static_cast<double>(step) * dt;

  // check if we can set the expected end-time based on 10% longer time than when all
  // adorbed mol-fractions are smaller than 1% of unity
  if (autoSteps)
  {
    double tolerance = 0.0;
    for (size_t j = 0; j < Ncomp; ++j)
    {
      tolerance =
          std::max(tolerance, std::abs((y[Ngrid * Ncomp + j] / (components[j].Yi0)) - 1.0));
    }

    // consider 1% as being visibily indistinguishable from 'converged'
    // use a 10% longer time for display purposes
    if (tolerance < 0.01)
    {
      std::cout << "\nConvergence criteria reached, running 10% longer\n\n" << std::endl;
      Nsteps = static_cast<size_t>(1.1 * static_cast<double>(step));
      autoSteps = false;
    }
  }

  // SSP-RK Step 1
  // ======================================================================

  // calculate the derivatives Dq/dt and Dp/dt based on Qeq, Q, V, and P
  computeFirstDerivatives(Dqdt, DPdt, DTdt, Dydt, Qeq, Q, V, P, T, y);

  // Dqdt and Dpdt are calculated at old time step
  // make estimate for the new loadings and new gas phase partial pressures
  // first iteration is made using the Explicit Euler scheme
  for (size_t i = 0; i < Ngrid + 1; ++i)
  {
    for (size_t j = 0; j < Ncomp; ++j)
    {
      Qnew[i * Ncomp + j] = std::max(Q[i * Ncomp + j] + dt * Dqdt[i * Ncomp + j], 0.0);
    }
    
    // Loop for all components except Carrier gas
    double sum_y = 0.0;
    for (size_t j = 1; j < Ncomp; ++j)
    {
      ynew[i * Ncomp + j] = std::max(y[i * Ncomp + j] + dt * Dydt[i * Ncomp + j], 0.0);
      ynew[i * Ncomp + j] = std::min(ynew[i * Ncomp + j], 1.0);
      sum_y += ynew[i * Ncomp + j]; 
    }

    // Carrier gas mole fraction is computed using molefractions of other components
    ynew[i * Ncomp + carrierGasComponent] = std::max((1.0 - sum_y), 0.0);
    
    // Add equation for Tnew and Pnew
    Tnew[i] = std::max(T[i] + dt * DTdt[i], 0.0);
    Pnew[i] = std::max(P[i] + dt * DPdt[i], 0.0);
  }

  computeEquilibriumLoadings();
  for (size_t i; i<Ngrid+1; ++i)
  {
    if (Pnew[i]<=0.0){
      throw std::runtime_error("Error: Pressure becoming zero\n"); 
    }
  }
  computeVelocity();

  // SSP-RK Step 2
  // ======================================================================

  // calculate new derivatives at new (current) timestep
  // calculate the derivatives Dq/dt and Dp/dt based on Qeq, Q, V, and P at new (current) timestep
  computeFirstDerivatives(Dqdtnew, DPdtnew, DTdtnew, Dydtnew,  Qeqnew, Qnew, Vnew, Pnew, Tnew, ynew);

  for (size_t i = 0; i < Ngrid + 1; ++i)
  {
    for (size_t j = 0; j < Ncomp; ++j)
    {
      Qnew[i * Ncomp + j] = std::max(0.75 * Q[i * Ncomp + j] + 0.25 * Qnew[i * Ncomp + j] + 0.25 * dt * Dqdtnew[i * Ncomp + j], 0.0);
    }

    double sum_y = 0.0;
    for (size_t j = 1; j < Ncomp; ++j)
    {
      ynew[i * Ncomp + j] = std::max(0.75 * y[i * Ncomp + j] + 0.25 * ynew[i * Ncomp + j] + 0.25 * dt * Dydtnew[i * Ncomp + j], 0.0);
      ynew[i * Ncomp + j] = std::min(ynew[i * Ncomp + j], 1.0);
      sum_y += ynew[i * Ncomp + j];
    }
    // Carrier gas mole fraction is computed using molefractions of other components
    ynew[i * Ncomp + carrierGasComponent] = std::max((1.0 - sum_y), 0.0);

    // Add equation for Tnew and Pnew
    Tnew[i] = std::max(0.75 * T[i] + 0.25 * Tnew[i] + 0.25 * dt * DTdtnew[i], 0.0);
    Pnew[i] = std::max(0.75 * P[i] + 0.25 * Pnew[i] + 0.25 * dt * DPdtnew[i], 0.0);
  }

  computeEquilibriumLoadings();

  for (size_t i; i<Ngrid+1; ++i)
  {
    if (Pnew[i]<=0.0){
      throw std::runtime_error("Error: Pressure becoming zero\n"); 
    }
  }
  computeVelocity();

  // SSP-RK Step 3
  // ======================================================================

  // calculate new derivatives at new (current) timestep
  // calculate the derivatives Dq/dt and Dp/dt based on Qeq, Q, V, and P at new (current) timestep
  computeFirstDerivatives(Dqdtnew, DPdtnew, DTdtnew, Dydtnew,  Qeqnew, Qnew, Vnew, Pnew, Tnew, ynew);

  for (size_t i = 0; i < Ngrid + 1; ++i)
  {
    for (size_t j = 0; j < Ncomp; ++j)
    {
      Qnew[i * Ncomp + j] = std::max((1.0 / 3.0) * Q[i * Ncomp + j] + (2.0 / 3.0) * Qnew[i * Ncomp + j] +
                            (2.0 / 3.0) * dt * Dqdtnew[i * Ncomp + j], 0.0);
    }

    double sum_y = 0.0;
    for (size_t j = 1; j < Ncomp; ++j)
    {
      ynew[i * Ncomp + j] = std::max((1.0 / 3.0) * y[i * Ncomp + j] + (2.0 / 3.0) * ynew[i * Ncomp + j] +
                            (2.0 / 3.0) * dt * Dydtnew[i * Ncomp + j], 0.0);
      ynew[i * Ncomp + j] = std::min(ynew[i * Ncomp + j], 1.0);
      sum_y += ynew[i * Ncomp + j];
    }
    
    // Carrier gas mole fraction is computed using molefractions of other components
    ynew[i * Ncomp + carrierGasComponent] = std::max((1.0 - sum_y), 0.0);

    Tnew[i] = std::max((1.0 / 3.0) * T[i] + (2.0 / 3.0) * Tnew[i] + (2.0 / 3.0) * dt * DTdtnew[i], 0.0);
    Pnew[i] = std::max((1.0 / 3.0) * P[i] + (2.0 / 3.0) * Pnew[i] + (2.0 / 3.0) * dt * DPdtnew[i], 0.0);
  }

  computeEquilibriumLoadings();

  for (size_t i; i<Ngrid+1; ++i)
  {
    if (Pnew[i]<=0.0){
      throw std::runtime_error("Error: Pressure becoming zero \n"); 
    }
  }
  computeVelocity();
  
  
  // update to the new time step
  std::copy(Qnew.begin(), Qnew.end(), Q.begin());
  std::copy(Pnew.begin(), Pnew.end(), P.begin());
  std::copy(Tnew.begin(), Tnew.end(), T.begin());
  std::copy(Qeqnew.begin(), Qeqnew.end(), Qeq.begin());
  std::copy(ynew.begin(), ynew.end(), y.begin());
  std::copy(Vnew.begin(), Vnew.end(), V.begin());

  // pulse boundary condition
  if (pulse == true)
  {
    if (t*L/v_in > tpulse)
    {
      for (size_t j = 0; j < Ncomp; ++j)
      {
        if (j == carrierGasComponent)
        {
          y[0 * Ncomp + j] = p_total/p_total;
        }
        else
        {
          y[0 * Ncomp + j] = 0.0;
        }
      }
    }
  }
}

void Breakthrough::computeEquilibriumLoadings()
{
  // Hassan modification
  // extended lagnmuir model instead of IAST
  std::vector <double> b(Ncomp);
  std::vector <double> d(Ncomp);
  double den_b = 1.0;
  double den_d = 1.0;

  // calculate new equilibrium loadings Qeqnew corresponding to the new timestep
  for(size_t i = 0; i < Ngrid + 1; ++i)
  {
    // estimation of total pressure Pt at each grid point from partial pressures
    // Hassan: Pt is a dummy variable to be used only for mixture prediction
    // Pt[i] = 0.0;
    Pt[i] = std::max(0.0, Pnew[i]);

    // for(size_t j = 0; j < Ncomp; ++j)
    // {
    //   Pt[i] += std::max(0.0, Pnew[i * Ncomp + j]);
    // }

    // compute gas-phase mol-fractions
    // force the gas-phase mol-fractions to be positive and normalized
    // double sum = 0.0;
    // for(size_t j = 0; j < Ncomp; ++j)
    // {
    //   Yi[j] = std::max(Pnew[i * Ncomp + j], 0.0);
    //   sum += Yi[j];
    // }

    // for(size_t j = 0; j < Ncomp; ++j)
    // {
    //   // Yi[j] /= sum;
    //   Yi[j] = std::max(ynew[i * Ncomp + j], 0.0);
    // }
    den_b = 1.0;
    den_d = 1.0;
    for(size_t j = 0; j < Ncomp; ++j)
    {
      // Yi[j] /= sum;
      b[j] = b_0[j] * std::exp(-del_H[j]/R * (1.0/Tnew[i]/T_gas - 1.0/T_ref));
      d[j] = d_0[j] * std::exp(-del_H[j]/R * (1.0/Tnew[i]/T_gas - 1.0/T_ref));

      den_b += (b[j]*ynew[i * Ncomp + j]*Pt[i]*p_total);
      den_d += (d[j]*ynew[i * Ncomp + j]*Pt[i]*p_total);
    }


    // use Yi and Pt[i] to compute the loadings in the adsorption mixture via mixture prediction
    // iastPerformance += mixture.predictMixture(Yi, Pt[i]*p_total, Xi, Ni, 
    //     &cachedP0[i * Ncomp * maxIsothermTerms], &cachedPsi[i * maxIsothermTerms]);

    for(size_t j = 0; j < Ncomp; ++j)
    {
      // Qeqnew[i * Ncomp + j] = Ni[j];
      Qeqnew[i * Ncomp + j] = sat_q_b[j]*b[j] * ynew[i * Ncomp + j]*Pt[i]*p_total / den_b
                            + sat_q_d[j]*d[j] * ynew[i * Ncomp + j]*Pt[i]*p_total / den_d;
    }
  }

  // check the total pressure at the outlet, it should not be negative
  if (Pt[Pt.size()-1] < 0.0)
  {
    throw std::runtime_error("Error: pressure gradient is too large/ Or some other problem (negative outlet pressure)\n");
  }
}


// calculate the derivatives Dq/dt and Dp/dt along the column
void Breakthrough::computeFirstDerivatives(std::vector<double> &dqdt,
                                           std::vector<double> &dpdt,
                                           std::vector<double> &dTdt,
                                           std::vector<double> &dydt,
                                           const std::vector<double> &q_eq,
                                           const std::vector<double> &q,
                                           const std::vector<double> &v,
                                           const std::vector<double> &p,
                                          //  const std::vector<double> &T_arg,
                                          //  const std::vector<double> &y_arg)
                                          const std::vector<double> &Temp,
                                           const std::vector<double> &y_vec)
{
  double idx = 1.0 / dx;
  double idx2 = 1.0 / (dx * dx);

  // The following variables have already been initialized in the initialize function
  // double K_z = 0.09;                        // Thermal conductivity of gas [J/mol/K]
  // double C_ps = 750.0;                      // Heat capacity of adsorbent [J/kg/K]
  // double C_pg = 35.8;                       // Heat capacity of gas [J/mol/K]
  // double C_pa = 35.8;                       // Heat capacity of adsorbate [J/mol/K]
  // double mu = 1.13e-05;                     // Viscoisty of gas [Pa.s]
  // double r_p = 5.0e-03;                     // Radius of adsorbent particles [m]  
  
  // %%%%%%%%%%%%%%%%%%%% Variables required for balances equations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  std::vector <double> ro_g(Ngrid+1);
  std::vector <double> sink_term(Ngrid+1);
  std::vector <double> kinetic_term(Ngrid+1);
  std::vector <double> p_dum(Ngrid+1);
  // std::vector <double> Temp(Ngrid+1);
  // std::vector <double> y_vec(Ngrid+1);

  // Variables for finite differences
  std::vector <double> dTdx(Ngrid+1, 0.0);
  std::vector <double> d2Tdx2(Ngrid+1, 0.0);
  std::vector <double> dydx((Ngrid+1)*Ncomp, 0.0);
  std::vector <double> d2ydx2((Ngrid+1)*Ncomp, 0.0);
  std::vector <double> PvT(Ngrid+1, 0.0);
  std::vector <double> Pv(Ngrid+1, 0.0);
  std::vector <double> dPdx(Ngrid+1, 0.0);
  
  double phi = R*rho_p*Q_s0*T_gas*(1-epsilon)/epsilon/p_total;
  double dTdt1, dTdt2, dTdt3;
  double dPdt1, dPdt2, dPdt3;
  double dydt1, dydt2, dydt3;

  std::copy(p.begin(), p.end(), p_dum.begin());
  // std::copy(y_arg.begin(), y_arg.end(), y_vec.begin());
  // std::copy(T_arg.begin(), T_arg.end(), Temp.begin());

  // Calculation of variable properties
  double sum_q; 
  for(size_t i = 0; i < Ngrid+1; i++)
  {
    sum_q = 0.0;
    for(size_t j = 0; j < Ncomp; ++j)
    {
      sum_q += q[i * Ncomp + j];      // Sum of adsorbent loading
    }
    ro_g[i] = (p_dum[i]*p_total) / R / Temp[i] / T_gas;
    sink_term[i] = (1 - epsilon) * (rho_p*C_ps + rho_p*sum_q*C_pa)
                  + (epsilon*ro_g[i]*C_pg);
  }

  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Inlet Boundary Pressure correction %%%%%%%%%%%%%%
  
  double vis_term = 150.0 * mu * std::pow((1-epsilon), 2) 
                    / 4.0 / std::pow(r_p, 2) / std::pow(epsilon, 2);

  double sum = 0.0;
  for(size_t j = 0; j < Ncomp; ++j)
  {
    sum += (MW[j] * y_vec[0 * Ncomp + j]);
  }
  
  // Calculate Inlet pressure using inlet velocity and 2nd grid point pressure
  // based on Ergun's equation
  // p_dum[0] = p_dum[1] + (vis_term*v_in + kinetic_term[1]*std::pow(v_in, 2.0))*dx;
  p_dum[0] = ((vis_term*v_in*dx*L/p_total) + p_dum[1]) 
          / (1.0 - (dx*L/R/Temp[0]/T_gas)*sum*((1.75*(1-epsilon))/2.0/r_p/epsilon)*std::pow(v_in, 2.0));
  if (p_dum[p_dum.size()-1] > p_dum[p_dum.size()-2]){
    p_dum[p_dum.size()-1] = p_dum[p_dum.size()-2];
  }

  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BC Check %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // BCs check (This check ensures that variables at the boundary, also implemented in MATLAB code)
  // In order to do that I need modifiable copy of vectors y_vec and Temp. Unfotunately, if I am initializing those copies within this function,...
  // I get the C++ error : free(): invalid pointer ...
  // I am a new to C++ and trying to address this issue. 
  // Thanks

  // for(size_t j = 0; j < Ncomp; ++j)
  // {
  //   y_vec[Ngrid * Ncomp + j] = y_vec[Ngrid-1 * Ncomp + j];
  // }
  // Temp[Ngrid] = Temp[Ngrid-1];
  p_dum[Ngrid] = p_dum[Ngrid-1];

  if (p_dum[Ngrid-1] >= 1.0){
    p_dum[Ngrid] = 1.0;
  }else{
    p_dum[Ngrid] = p_dum[Ngrid-1];
  }

  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Finite Difference Calculation %%%%%%%%%%%%%%%%%%%
  // First order differences
  for(size_t i = 1; i < Ngrid+1; i++)
  {
    for(size_t j = 0; j < Ncomp; ++j)
    {
      dydx[i * Ncomp + j] = (y_vec[i * Ncomp + j] - y_vec[(i-1) * Ncomp + j]) * idx;
    }
    dTdx[i] = (Temp[i] - Temp[i-1]) * idx;
    dPdx[i] = (p_dum[i] - p_dum[i-1]) * idx;

    Pv[i] = (p_dum[i]*v[i] - p_dum[i-1]*v[i-1]) * idx;
    PvT[i] = (p_dum[i]*v[i]/Temp[i] - p_dum[i-1]*v[i-1]/Temp[i-1]) * idx;
  }

  // Second order differences
  for(size_t i = 1; i < Ngrid; i++)     // For middle nodes
  {
    for(size_t j = 0; j < Ncomp; ++j)
    {
      d2ydx2[i * Ncomp + j] = (y_vec[(i+1) * Ncomp + j] - 2.0*y_vec[i * Ncomp + j] + y_vec[(i-1) * Ncomp + j]) * idx2;
    }
    
    d2Tdx2[i] = (Temp[i+1] - 2.0*Temp[i] + Temp[i-1]) * idx2;
  }

  // For last node
  for(size_t j = 0; j < Ncomp; ++j)
  {
    d2ydx2[Ngrid * Ncomp + j] = (y_vec[(Ngrid-1) * Ncomp + j] - y_vec[Ngrid * Ncomp + j]) * idx2;
  }
  d2Tdx2[Ngrid] = (Temp[Ngrid-1] - Temp[Ngrid]) * idx2;
  
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Balance equations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // first/inlet  gridpoint, the conditions remains same
  dpdt[0] = 0.0;
  dTdt[0] = 0.0;
  for(size_t j = 0; j < Ncomp; ++j)
  {
    dydt[0 * Ncomp + j] = 0.0;
    dqdt[0 * Ncomp + j] = (components[j].Kl * L/v_in) * (q_eq[0 * Ncomp + j]/Q_s0 - q[0 * Ncomp + j]);
  }
  
  // middle gridpoints 
  for(size_t i = 1; i < Ngrid; i++)
  {
    // %%%%%%%%%%%%%%%%%%%%%%% Solid mass balance %%%%%%%%%%%%%%%%%%%%%%%%%%
    for(size_t j = 0; j < Ncomp; ++j)
    {
      dqdt[i * Ncomp + j] = (components[j].Kl * L/v_in) * (q_eq[i * Ncomp + j]/Q_s0 - q[i * Ncomp + j]);
    }
    
    // %%%%%%%%%%%%%%%%%%%%%%%%% Temperature balance %%%%%%%%%%%%%%%%%%%%%%%%
    dTdt1 = K_z/v_in/L * d2Tdx2[i] / sink_term[i];
    dTdt2 = -(epsilon*C_pg*p_total/R/T_gas) * (Pv[i] - Temp[i]*PvT[i]) / sink_term[i];
    
    dTdt3 = 0.0; 
    for(size_t j = 0; j < Ncomp; ++j)
    {
      dTdt3 += -(1-epsilon) * rho_p * Q_s0 * dqdt[i * Ncomp + j] * del_H[j] / T_gas / sink_term[i];
    }
    
    // dTdt[i] = 0.0 * (dTdt1 + dTdt2 + dTdt3);
    dTdt[i] = (dTdt1 + dTdt2 + dTdt3);

    // %%%%%%%%%%%%%%%%%%%%%%%%% Total Pressure  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    dPdt1 = -Temp[i] * PvT[i];
    dPdt2 = p_dum[i] * dTdt[i] / Temp[i];

    dPdt3 = 0.0; 
    for(size_t j = 0; j < Ncomp; ++j)
    {
      dPdt3 += -phi * Temp[i] * dqdt[i * Ncomp + j];
    }

    dpdt[i] = dPdt1 + dPdt2 + dPdt3;
    
    // %%%%%%%%%%%%%%%%%%%%%%%%%% Mole balance %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // For all components except carrier gas
    for(size_t j = 1; j < Ncomp; ++j)
    {
      double dyPVT;
      dyPVT = ((y_vec[i * Ncomp + j]*p_dum[i]*v[i]/Temp[i]) 
              - (y_vec[(i-1) * Ncomp + j]*p_dum[i-1]*v[i-1]/Temp[i-1])) * idx; 
      dydt1 = - (Temp[i]/p_dum[i]) * (dyPVT - y_vec[i * Ncomp + j]*PvT[i]);
      
      dydt2 = (components[j].D/L/v_in) * (d2ydx2[i * Ncomp + j] 
                                + (dPdx[i]*dydx[i * Ncomp + j])/p_dum[i] 
                                - (dTdx[i]*dydx[i * Ncomp + j])/Temp[i]);
      double res_dqdt = 0.0;
      for(size_t k = 0; k < Ncomp; ++k)
      {
        if (k != j)
        {
          res_dqdt += dqdt[i * Ncomp + k];
        }
      }

      dydt3 = phi * Temp[i] / p_dum[i]
              * ((y_vec[i * Ncomp + j] - 1)*dqdt[i * Ncomp + j]
              + y_vec[i * Ncomp + j] * res_dqdt);
      
      dydt[i * Ncomp + j] = dydt1 + dydt2 + dydt3;
    }
  }
  // Outlet gridpoints 
  // %%%%%%%%%%%%%%%%%%%%%%% Solid mass balance %%%%%%%%%%%%%%%%%%%%%%%%%%
  for(size_t j = 0; j < Ncomp; ++j)
  {
    dqdt[Ngrid * Ncomp + j] = (components[j].Kl * L/v_in) * (q_eq[Ngrid * Ncomp + j]/Q_s0 - q[Ngrid * Ncomp + j]);
  }
  
  // %%%%%%%%%%%%%%%%%%%%%%%%% Temperature balance %%%%%%%%%%%%%%%%%%%%%%%%
  dTdt[Ngrid] = dTdt[Ngrid-1];

  // %%%%%%%%%%%%%%%%%%%%%%%%% Total Pressure  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  dpdt[Ngrid] = dpdt[Ngrid-1];
  
  // %%%%%%%%%%%%%%%%%%%%%%%%%% Mole balance %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // For all components except carrier gas
  for(size_t j = 1; j < Ncomp; ++j)
  {
    dydt[Ngrid * Ncomp + j] = dydt[(Ngrid-1) * Ncomp + j];
  }
}

// calculate new velocity Vnew from Qnew, Qeqnew, Pnew, Pt
void Breakthrough::computeVelocity()
{
  double idx = 1.0 / dx;
  // double mu = 1.13e-05;
  // double r_p = 5.0e-03;

  double vis_term = 150.0 * mu * std::pow((1-epsilon), 2) 
                    / 4.0 / std::pow(r_p, 2) / std::pow(epsilon, 2);
  
  std::vector <double> ro_g(Ngrid+1);
  std::vector <double> kinetic_term(Ngrid+1);
  std::vector <double> dPdx(Ngrid+1);
  std::vector <double> p_dum(Ngrid+1); 
  std::copy(Pnew.begin(), Pnew.end(), p_dum.begin());

  for (size_t i = 0; i < Ngrid+1; ++i)
  {
    double sum = 0.0;
    for(size_t j = 0; j < Ncomp; ++j)
    {
      sum += (MW[j] * ynew[i * Ncomp + j]);
    }
    if (p_dum[i]<=0.0){
      throw std::runtime_error("Error: Pressure becoming zero\n"); 
    }
    ro_g[i] = p_dum[i]*p_total / R / Tnew[i] / T_gas;
    kinetic_term[i] = (ro_g[i] * sum) * (1.75*(1-epsilon)) / 2.0 / r_p / epsilon;
  }

  // Calculate Inlet pressure using inlet velocity and 2nd grid point pressure
  // based on Ergun's equation
  double sum = 0.0;
  for(size_t j = 0; j < Ncomp; ++j)
  {
    sum += (MW[j] * ynew[0 * Ncomp + j]);
  }
 
  // p_dum[0] = p_dum[1] + (vis_term*v_in + kinetic_term[1]*std::pow(v_in, 2.0))*dx;
  p_dum[0] = ((vis_term*v_in*dx*L/p_total) + p_dum[1]) 
          / (1.0 - (dx*L/R/Tnew[0]/T_gas)*sum*((1.75*(1-epsilon))/2.0/r_p/epsilon)*std::pow(v_in, 2.0));
  if (p_dum[p_dum.size()-1] > p_dum[p_dum.size()-2]){
    p_dum[p_dum.size()-1] = p_dum[p_dum.size()-2];
  }

  // Pressure gradient
  dPdx[0] = 0.0;      // Inlet grid point
  for(size_t i = 1; i < Ngrid+1; ++i)
  {
    dPdx[i] = (p_dum[i] - p_dum[i-1]) * idx;
  }

  // Calculate Velocity based on Ergun's equation
  Vnew[0] = v_in/v_in;     // first grid point
  
  int v_sign;
  for(size_t i = 1; i < Ngrid+1; ++i)
  {  
    if (dPdx[i] <= 0.0)
    {
      v_sign = 1.0;
    }
    else
    {
      v_sign = -1.0;
    }
    
    Vnew[i] = v_sign * (-vis_term + std::pow(
              (std::abs(std::pow(vis_term, 2) + 4.0 * kinetic_term[i] * std::abs(dPdx[i]*p_total/L))), 0.5))
              / 2.0 / kinetic_term[i] / v_in;
    if (std::isnan(Vnew[i]))
    {
      std::cout << "Nan encountered" << std::endl;
    } 
  }
  
  // middle gridpoints
  // for(size_t i = 1; i < Ngrid; ++i)  
  // {
    
    // // sum = derivative at the actual gridpoint i
    // double sum = 0.0;
    // for(size_t j = 0; j < Ncomp; ++j)
    // {
    //   sum = sum - prefactor[j] * (Qeqnew[i * Ncomp + j] - Qnew[i * Ncomp + j]) +
    //         components[j].D * (Pnew[(i - 1) * Ncomp + j] - 2.0 * Pnew[i * Ncomp + j] + Pnew[(i + 1) * Ncomp + j]) * idx2;
    // }
  
    // // explicit version
    // Vnew[i] = Vnew[i - 1] + dx * (sum - Vnew[i - 1] * dptdx) / Pt[i];
  // }
  
  // last grid point
  // double sum = 0.0;
  // for(size_t j = 0; j < Ncomp; ++j)
  // {
  //   sum = sum - prefactor[j] * (Qeqnew[Ngrid * Ncomp + j] - Qnew[Ngrid * Ncomp + j]) +
  //         components[j].D * (Pnew[(Ngrid - 1) * Ncomp + j] - Pnew[Ngrid * Ncomp + j]) * idx2;
  // }
  
  // // explicit version
  // Vnew[Ngrid] = Vnew[Ngrid-1] + dx * (sum - Vnew[Ngrid - 1] * dptdx) / Pt[Ngrid];
}

std::string Breakthrough::repr() const
{
  std::string s;
  s += "Column properties\n";
  s += "=======================================================\n";
  s += "Display-name:                          " + displayName + "\n";
  s += "Temperature:                           " + std::to_string(T_gas) + " [K]\n";
  s += "Column length:                         " + std::to_string(L) + " [m]\n";
  s += "Column void-fraction:                  " + std::to_string(epsilon) + " [-]\n";
  s += "Particle density:                      " + std::to_string(rho_p) + " [kg/m^3]\n";
  s += "Total pressure:                        " + std::to_string(p_total) + " [Pa]\n";
  s += "Pressure gradient:                     " + std::to_string(dptdx) + " [Pa/m]\n";
  s += "Column entrance interstitial velocity: " + std::to_string(v_in) + " [m/s]\n";
  s += "\n\n";

  s += "Breakthrough settings\n";
  s += "=======================================================\n";
  s += "Number of time steps:          " + std::to_string(Nsteps) + "\n";
  s += "Print every step:              " + std::to_string(printEvery) + "\n";
  s += "Write data every step:         " + std::to_string(writeEvery) + "\n";
  s += "\n\n";

  s += "Integration details\n";
  s += "=======================================================\n";
  s += "Time step:                     " + std::to_string(dt) + " [s]\n";
  s += "Number of column grid points:  " + std::to_string(Ngrid) + "\n";
  s += "Column spacing:                " + std::to_string(dx) + " [m]\n";
  s += "\n\n";

  s += "Component data\n";
  s += "=======================================================\n";
  s += "maximum isotherm terms:        " + std::to_string(maxIsothermTerms) + "\n";
  for(size_t i = 0; i < Ncomp; ++i)
  {
    s += components[i].repr() + "\n";
  }
  return s;
}

void Breakthrough::createPlotScript()
{
  #if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
    std::ofstream stream_graphs("make_graphs.bat");
    stream_graphs << "set PATH=%PATH%;C:\\Program Files\\gnuplot\\bin;C:\\Program Files\\ffmpeg-master-latest-win64-gpl\\bin;C:\\Program Files\\ffmpeg\\bin\n";
    stream_graphs << "gnuplot.exe plot_breakthrough\n";
  #else
    std::ofstream stream_graphs("make_graphs");
    stream_graphs << "#!/bin/sh\n";
    stream_graphs << "export LC_ALL='en_US.UTF-8'\n";
    stream_graphs << "cd -- \"$(dirname \"$0\")\"\n";
    stream_graphs << "gnuplot plot_breakthrough\n";
  #endif

  #if (__cplusplus >= 201703L)
    std::filesystem::path path{"make_graphs"};
    std::filesystem::permissions(path, std::filesystem::perms::owner_exec, std::filesystem::perm_options::add);
  #else 
    chmod("make_graphs", S_IRWXU);
  #endif

  std::ofstream stream("plot_breakthrough");
  stream << "set encoding utf8\n";
  #if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
    stream << "set xlabel 'Dimensionless time, {/Arial-Italic τ}={/Arial-Italic tv/L} / [-]' font \"Arial,14\"\n";
    stream << "set ylabel 'Concentration exit gas, {/Arial-Italic c}_i/{/Arial-Italic c}_{i,0} / [-]' offset 0.0,0 font \"Arial,14\"\n";
    stream << "set key outside top center horizontal samplen 2.5 height 0.5 spacing 1.5 font 'Arial, 10'\n";
  #else
    stream << "set xlabel 'Dimensionless time, {/Helvetica-Italic τ}={/Helvetica-Italic tv/L} / [-]' font \"Helvetica,18\"\n";
    stream << "set ylabel 'Concentration exit gas, {/Helvetica-Italic c}_i/{/Helvetica-Italic c}_{i,0} / [-]' offset 0.0,0 font \"Helvetica,18\"\n";
    stream << "set key outside top center horizontal samplen 2.5 height 0.5 spacing 1.5 font 'Helvetica, 10'\n";
  #endif
  stream << "set bmargin 4\n";
  stream << "set yrange[0:]\n";

  stream << "set key title '" << displayName << " {/:Italic T}=" << T_gas << " K, {/:Italic p_t}=" << p_total*1e-3 << " kPa'\n";

  stream << "set output 'breakthrough_dimensionless.pdf'\n";
  stream << "set term pdf color solid\n";

  // colorscheme from book 'gnuplot in action', listing 12.7
  stream << "set linetype 1 pt 5 ps 1 lw 4 lc rgb '0xee0000'\n";
  stream << "set linetype 2 pt 7 ps 1 lw 4 lc rgb '0x008b00'\n";
  stream << "set linetype 3 pt 9 ps 1 lw 4 lc rgb '0x0000cd'\n";
  stream << "set linetype 4 pt 11 ps 1 lw 4 lc rgb '0xff3fb3'\n";
  stream << "set linetype 5 pt 13 ps 1 lw 4 lc rgb '0x00cdcd'\n";
  stream << "set linetype 6 pt 15 ps 1 lw 4 lc rgb '0xcd9b1d'\n";
  stream << "set linetype 7 pt  4 ps 1 lw 4 lc rgb '0x8968ed'\n";
  stream << "set linetype 8 pt  6 ps 1 lw 4 lc rgb '0x8b8b83'\n";
  stream << "set linetype 9 pt  8 ps 1 lw 4 lc rgb '0x00bb00'\n";
  stream << "set linetype 10 pt 10 ps 1 lw 4 lc rgb '0x1e90ff'\n";
  stream << "set linetype 11 pt 12 ps 1 lw 4 lc rgb '0x8b2500'\n";
  stream << "set linetype 12 pt 14 ps 1 lw 4 lc rgb '0x000000'\n";

  stream << "ev=1\n";
  stream << "plot \\\n";
  for (size_t i = 0; i < Ncomp; i++)
  {
    std::string fileName = "component_" + std::to_string(i) + "_" + components[i].name + ".data";
    stream << "    " << "\"" << fileName << "\"" << " us ($1):($3) every ev" << " title \""
           << components[i].name << " (y_i=" << components[i].Yi0 << ")\""
           << " with li lt " << i+1 << (i < Ncomp - 1 ? ",\\" : "") << "\n";
  }
  stream << "set output 'breakthrough.pdf'\n";
  #if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
     stream << "set xlabel 'Time, {/Arial-Italic t} / [min.]' font \"Arial,14\"\n";
  #else
     stream << "set xlabel 'Time, {/Helvetica-Italic t} / [min.]' font \"Helvetica,18\"\n";
  #endif
  stream << "plot \\\n";
  for (size_t i = 0; i < Ncomp; i++)
  {
    std::string fileName = "component_" + std::to_string(i) + "_" + components[i].name + ".data";
    stream << "    " << "\"" << fileName << "\"" << " us ($2):($3) every ev" << " title \""
           << components[i].name << " (y_i=" << components[i].Yi0 << ")\""
           << " with li lt " << i+1 << (i < Ncomp - 1 ? ",\\" : "") << "\n";
  }
}

void Breakthrough::createMovieScripts()
{
  #if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
    std::ofstream makeMovieStream("make_movies.bat");
    makeMovieStream << "CALL make_movie_V.bat %1 %2 %3 %4\n";
    makeMovieStream << "CALL make_movie_Pt.bat %1 %2 %3 %4\n";
    makeMovieStream << "CALL make_movie_Q.bat %1 %2 %3 %4\n";
    makeMovieStream << "CALL make_movie_Qeq.bat %1 %2 %3 %4\n";
    makeMovieStream << "CALL make_movie_P.bat %1 %2 %3 %4\n";
    makeMovieStream << "CALL make_movie_Pnorm.bat %1 %2 %3 %4\n";
    makeMovieStream << "CALL make_movie_Dpdt.bat %1 %2 %3 %4\n";
    makeMovieStream << "CALL make_movie_Dqdt.bat %1 %2 %3 %4\n";
  #else
    std::ofstream makeMovieStream("make_movies");
    makeMovieStream << "#!/bin/sh\n";
    makeMovieStream << "cd -- \"$(dirname \"$0\")\"\n";
    makeMovieStream << "./make_movie_V \"$@\"\n";
    makeMovieStream << "./make_movie_Pt \"$@\"\n";
    makeMovieStream << "./make_movie_Q \"$@\"\n";
    makeMovieStream << "./make_movie_Qeq \"$@\"\n";
    makeMovieStream << "./make_movie_P \"$@\"\n";
    makeMovieStream << "./make_movie_Pnorm \"$@\"\n";
    makeMovieStream << "./make_movie_Dpdt \"$@\"\n";
    makeMovieStream << "./make_movie_Dqdt \"$@\"\n";
  #endif

  #if (__cplusplus >= 201703L)
    std::filesystem::path path{"make_movies"};
    std::filesystem::permissions(path, std::filesystem::perms::owner_exec, std::filesystem::perm_options::add);
  #else 
    chmod("make_movies", S_IRWXU);
  #endif

  createMovieScriptColumnV();
  createMovieScriptColumnPt();
  createMovieScriptColumnQ();
  createMovieScriptColumnQeq();
  createMovieScriptColumnP();
  createMovieScriptColumnDpdt();
  createMovieScriptColumnDqdt();
  createMovieScriptColumnPnormalized();
}

// -crf 18: the range of the CRF scale is 0–51, where 0 is lossless, 23 is the default, 
//          and 51 is worst quality possible; 18 is visually lossless or nearly so.
// -pix_fmt yuv420p: needed on apple devices
std::string movieScriptTemplate(std::string s)
{
  std::ostringstream stream;

  #if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
    stream << "del column_movie_" << s << ".mp4\n";
    stream << "set /A argVec[1]=1\n";
    stream << "set /A argVec[2]=1200\n";
    stream << "set /A argVec[3]=800\n";
    stream << "set /A argVec[4]=18\n";
    stream << "setlocal enabledelayedexpansion\n";
    stream << "set argCount=0\n";
    stream << "for %%x in (%*) do (\n";
    stream << "   set /A argCount+=1\n";
    stream << "   set \"argVec[!argCount!]=%%~x\"'n";
    stream << ")\n";
    stream << "set PATH=%PATH%;C:\\Program Files\\gnuplot\\bin;C:\\Program Files\\ffmpeg-master-latest-win64-gpl\\bin;C:\\Program Files\\ffmpeg\\bin\n";
    stream << "gnuplot.exe -c plot_column_" << s << " %argVec[1]% %argVec[2]% %argVec[3]% | ffmpeg.exe -f png_pipe -s:v \"%argVec[2]%,%argVec[3]%\" -i pipe: -c:v libx264 -pix_fmt yuv420p -crf %argVec[4]% -c:a aac column_movie_" << s + ".mp4\n";
  #else
    stream << "rm -f " << "column_movie_" << s << ".mp4\n";
    stream << "every=1\n";
    stream << "format=\"-c:v libx265 -tag:v hvc1\"\n";
    stream << "width=1200\n";
    stream << "height=800\n";
    stream << "quality=18\n";
    stream << "while getopts e:w:h:q:l flag\n";
    stream << "do\n";
    stream << "    case \"${flag}\" in\n";
    stream << "        e) every=${OPTARG};;\n";
    stream << "        w) width=${OPTARG};;\n";
    stream << "        h) height=${OPTARG};;\n";
    stream << "        q) quality=${OPTARG};;\n";
    stream << "        l) format=\"-c:v libx264\";;\n";
    stream << "    esac\n";
    stream << "done\n";
    stream << "gnuplot -c plot_column_" << s << " $every $width $height | ffmpeg -f png_pipe -s:v \"${width},${height}\" -i pipe: $format -pix_fmt yuv420p -crf $quality -c:a aac column_movie_" << s + ".mp4\n";
  #endif
  return stream.str();
}

void Breakthrough::createMovieScriptColumnV()
{
  #if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
    std::ofstream makeMovieStream("make_movie_V.bat");
  #else
    std::ofstream makeMovieStream("make_movie_V");
    makeMovieStream << "#!/bin/sh\n";
    makeMovieStream << "cd -- \"$(dirname \"$0\")\"\n";
  #endif
  makeMovieStream << movieScriptTemplate("V");

  #if (__cplusplus >= 201703L)
    std::filesystem::path path{"make_movie_V"};
    std::filesystem::permissions(path, std::filesystem::perms::owner_exec, std::filesystem::perm_options::add);
  #else
    chmod("make_movie_V", S_IRWXU);
  #endif

  std::ofstream stream("plot_column_V");

  stream << "set encoding utf8\n";
  #if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
    stream << "set terminal pngcairo size ARG2,ARG3 enhanced font 'Arial,10'\n";
    stream << "set xlabel 'Adsorber position / [m]' font 'Arial,14'\n";
    stream << "set ylabel 'Interstitial velocity, {/Arial-Italic v} / [m/s]' offset 0.0,0 font 'Arial,14'\n";
    stream << "set key outside top center horizontal samplen 2.5 height 0.5 spacing 1.5 font 'Arial, 10'\n";
  #else
    stream << "set terminal pngcairo size ARG2,ARG3 enhanced font 'Helvetica,10'\n";
    stream << "set xlabel 'Adsorber position / [m]' font 'Helvetica,18'\n";
    stream << "set ylabel 'Interstitial velocity, {/Helvetica-Italic v} / [m/s]' offset 0.0,0 font 'Helvetica,18'\n";
    stream << "set key outside top center horizontal samplen 2.5 height 0.5 spacing 1.5 font 'Helvetica, 10'\n";
  #endif

  // colorscheme from book 'gnuplot in action', listing 12.7
  stream << "set linetype 1 pt 5 ps 1 lw 4 lc rgb '0xee0000'\n";
  stream << "set linetype 2 pt 7 ps 1 lw 4 lc rgb '0x008b00'\n";
  stream << "set linetype 3 pt 9 ps 1 lw 4 lc rgb '0x0000cd'\n";
  stream << "set linetype 4 pt 11 ps 1 lw 4 lc rgb '0xff3fb3'\n";
  stream << "set linetype 5 pt 13 ps 1 lw 4 lc rgb '0x00cdcd'\n";
  stream << "set linetype 6 pt 15 ps 1 lw 4 lc rgb '0xcd9b1d'\n";
  stream << "set linetype 7 pt  4 ps 1 lw 4 lc rgb '0x8968ed'\n";
  stream << "set linetype 8 pt  6 ps 1 lw 4 lc rgb '0x8b8b83'\n";
  stream << "set linetype 9 pt  8 ps 1 lw 4 lc rgb '0x00bb00'\n";
  stream << "set linetype 10 pt 10 ps 1 lw 4 lc rgb '0x1e90ff'\n";
  stream << "set linetype 11 pt 12 ps 1 lw 4 lc rgb '0x8b2500'\n";
  stream << "set linetype 12 pt 14 ps 1 lw 4 lc rgb '0x000000'\n";

  stream << "set bmargin 4\n";
  stream << "set title '" << displayName << " {/:Italic T}=" << T_gas << " K, {/:Italic p_t}=" << p_total*1e-3 << " kPa'\n";
  stream << "stats 'column.data' us 2 nooutput\n";
  stream << "max=STATS_max\n";
  stream << "stats 'column.data' us 1 nooutput\n";
  stream << "set xrange[0:STATS_max]\n";
  stream << "set yrange[0:1.1*max]\n";
  stream << "ev=int(ARG1)\n";
  stream << "do for [i=0:int((STATS_blocks-2)/ev)] {\n";
  stream << "  plot \\\n";
  stream << "    " << "'column.data'" << " us 1:2 index ev*i notitle with li lt 1,\\\n";
  stream << "    " << "'column.data'" << " us 1:2 index ev*i notitle with po lt 1\n";
  stream << "}\n";
}


void Breakthrough::createMovieScriptColumnPt()
{
  #if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
    std::ofstream makeMovieStream("make_movie_Pt.bat");
  #else
    std::ofstream makeMovieStream("make_movie_Pt");
    makeMovieStream << "#!/bin/sh\n";
    makeMovieStream << "cd -- \"$(dirname \"$0\")\"\n";
  #endif
  makeMovieStream << movieScriptTemplate("Pt");

  #if (__cplusplus >= 201703L)
    std::filesystem::path path{"make_movie_Pt"};
    std::filesystem::permissions(path, std::filesystem::perms::owner_exec, std::filesystem::perm_options::add);
  #else
    chmod("make_movie_Pt", S_IRWXU);
  #endif

  std::ofstream stream("plot_column_Pt");

  stream << "set encoding utf8\n";
  #if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
    stream << "set terminal pngcairo size ARG2,ARG3 enhanced font 'Arial,10'\n";
    stream << "set xlabel 'Adsorber position / [m]' font 'Arial,14'\n";
    stream << "set ylabel 'Total Pressure, {/Arial-Italic p_t} / [Pa]' offset 0.0,0 font 'Arial,14'\n";
    stream << "set key outside top center horizontal samplen 2.5 height 0.5 spacing 1.5 font 'Arial, 10'\n";
  #else 
    stream << "set terminal pngcairo size ARG2,ARG3 enhanced font 'Helvetica,10'\n";
    stream << "set xlabel 'Adsorber position / [m]' font 'Helvetica,18'\n";
    stream << "set ylabel 'Total Pressure, {/Helvetica-Italic p_t} / [Pa]' offset 0.0,0 font 'Helvetica,18'\n";
    stream << "set key outside top center horizontal samplen 2.5 height 0.5 spacing 1.5 font 'Helvetica, 10'\n";
  #endif

  // colorscheme from book 'gnuplot in action', listing 12.7
  stream << "set linetype 1 pt 5 ps 1 lw 4 lc rgb '0xee0000'\n";
  stream << "set linetype 2 pt 7 ps 1 lw 4 lc rgb '0x008b00'\n";
  stream << "set linetype 3 pt 9 ps 1 lw 4 lc rgb '0x0000cd'\n";
  stream << "set linetype 4 pt 11 ps 1 lw 4 lc rgb '0xff3fb3'\n";
  stream << "set linetype 5 pt 13 ps 1 lw 4 lc rgb '0x00cdcd'\n";
  stream << "set linetype 6 pt 15 ps 1 lw 4 lc rgb '0xcd9b1d'\n";
  stream << "set linetype 7 pt  4 ps 1 lw 4 lc rgb '0x8968ed'\n";
  stream << "set linetype 8 pt  6 ps 1 lw 4 lc rgb '0x8b8b83'\n";
  stream << "set linetype 9 pt  8 ps 1 lw 4 lc rgb '0x00bb00'\n";
  stream << "set linetype 10 pt 10 ps 1 lw 4 lc rgb '0x1e90ff'\n";
  stream << "set linetype 11 pt 12 ps 1 lw 4 lc rgb '0x8b2500'\n";
  stream << "set linetype 12 pt 14 ps 1 lw 4 lc rgb '0x000000'\n";

  stream << "set bmargin 4\n";
  stream << "set title '" << displayName << " {/:Italic T}=" << T_gas << " K, {/:Italic p_t}=" << p_total*1e-3 << " kPa'\n";
  stream << "stats 'column.data' us 3 nooutput\n";
  stream << "max=STATS_max\n";
  stream << "stats 'column.data' us 1 nooutput\n";
  stream << "set xrange[0:STATS_max]\n";
  stream << "set yrange[0:1.1*max]\n";
  stream << "ev=int(ARG1)\n";
  stream << "do for [i=0:int((STATS_blocks-2)/ev)] {\n";
  stream << "  plot \\\n";
  stream << "    " << "'column.data'" << " us 1:3 index ev*i notitle with li lt 1,\\\n";
  stream << "    " << "'column.data'" << " us 1:3 index ev*i notitle with po lt 1\n";
  stream << "}\n";
}

void Breakthrough::createMovieScriptColumnQ()
{
  #if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
    std::ofstream makeMovieStream("make_movie_Q.bat");
  #else
    std::ofstream makeMovieStream("make_movie_Q");
    makeMovieStream << "#!/bin/sh\n";
    makeMovieStream << "cd -- \"$(dirname \"$0\")\"\n";
  #endif
  makeMovieStream << movieScriptTemplate("Q");

  #if (__cplusplus >= 201703L)
    std::filesystem::path path{"make_movie_Q"};
    std::filesystem::permissions(path, std::filesystem::perms::owner_exec, std::filesystem::perm_options::add);
  #else
    chmod("make_movie_Q", S_IRWXU);
  #endif

  std::ofstream stream("plot_column_Q");

  stream << "set encoding utf8\n";
  #if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
    stream << "set terminal pngcairo size ARG2,ARG3 enhanced font 'Arial,10'\n";
    stream << "set xlabel 'Adsorber position / [m]' font 'Arial,14'\n";
    stream << "set ylabel 'Concentration, {/Arial-Italic c}_i / [mol/kg]' offset 0.0,0 font 'Arial,14'\n";
    stream << "set key outside top center horizontal samplen 2.5 height 0.5 spacing 1.5 font 'Arial, 10'\n";
  #else
    stream << "set terminal pngcairo size ARG2,ARG3 enhanced font 'Helvetica,10'\n";
    stream << "set xlabel 'Adsorber position / [m]' font 'Helvetica,18'\n";
    stream << "set ylabel 'Concentration, {/Helvetica-Italic c}_i / [mol/kg]' offset 0.0,0 font 'Helvetica,18'\n";
    stream << "set key outside top center horizontal samplen 2.5 height 0.5 spacing 1.5 font 'Helvetica, 10'\n";
  #endif

  // colorscheme from book 'gnuplot in action', listing 12.7
  stream << "set linetype 1 pt 5 ps 1 lw 4 lc rgb '0xee0000'\n";
  stream << "set linetype 2 pt 7 ps 1 lw 4 lc rgb '0x008b00'\n";
  stream << "set linetype 3 pt 9 ps 1 lw 4 lc rgb '0x0000cd'\n";
  stream << "set linetype 4 pt 11 ps 1 lw 4 lc rgb '0xff3fb3'\n";
  stream << "set linetype 5 pt 13 ps 1 lw 4 lc rgb '0x00cdcd'\n";
  stream << "set linetype 6 pt 15 ps 1 lw 4 lc rgb '0xcd9b1d'\n";
  stream << "set linetype 7 pt  4 ps 1 lw 4 lc rgb '0x8968ed'\n";
  stream << "set linetype 8 pt  6 ps 1 lw 4 lc rgb '0x8b8b83'\n";
  stream << "set linetype 9 pt  8 ps 1 lw 4 lc rgb '0x00bb00'\n";
  stream << "set linetype 10 pt 10 ps 1 lw 4 lc rgb '0x1e90ff'\n";
  stream << "set linetype 11 pt 12 ps 1 lw 4 lc rgb '0x8b2500'\n";
  stream << "set linetype 12 pt 14 ps 1 lw 4 lc rgb '0x000000'\n";

  stream << "set bmargin 4\n";
  stream << "set key title '" << displayName << " {/:Italic T}=" << T_gas << " K, {/:Italic p_t}=" << p_total*1e-3 << " kPa'\n";
  stream << "stats 'column.data' nooutput\n";
  stream << "max = 0.0;\n";
  stream << "do for [i=4:STATS_columns:6] {\n";
  stream << "  stats 'column.data' us i nooutput\n";
  stream << "  if (max<STATS_max) {\n";
  stream << "    max=STATS_max\n";
  stream << "  }\n";
  stream << "}\n";
  stream << "stats 'column.data' us 1 nooutput\n";
  stream << "set xrange[0:STATS_max]\n";
  stream << "set yrange[0:1.1*max]\n";
  stream << "ev=int(ARG1)\n";
  stream << "do for [i=0:int((STATS_blocks-2)/ev)] {\n";
  stream << "  plot \\\n";
  for (size_t i = 0; i < Ncomp; i++)
  {
    stream << "    " << "'column.data'" << " us 1:" << std::to_string(4 + i * 6) << " index ev*i notitle " 
           << " with li lt " << i+1 << ",\\\n";
  }
  for (size_t i = 0; i < Ncomp; i++)
  {
    stream << "    " << "'column.data'" << " us 1:" << std::to_string(4 + i * 6) << " index ev*i title '"
           << components[i].name << " (y_i=" << components[i].Yi0 << ")'"
           << " with po lt " << i+1 << (i < Ncomp - 1 ? ",\\" : "") << "\n";
  }
  stream << "}\n";
}

void Breakthrough::createMovieScriptColumnQeq()
{
  #if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
    std::ofstream makeMovieStream("make_movie_Qeq.bat");
  #else
    std::ofstream makeMovieStream("make_movie_Qeq");
    makeMovieStream << "#!/bin/sh\n";
    makeMovieStream << "cd -- \"$(dirname \"$0\")\"\n";
  #endif
  makeMovieStream << movieScriptTemplate("Qeq");

  #if (__cplusplus >= 201703L)
    std::filesystem::path path{"make_movie_Qeq"};
    std::filesystem::permissions(path, std::filesystem::perms::owner_exec, std::filesystem::perm_options::add);
  #else
    chmod("make_movie_Qeq", S_IRWXU);
  #endif

  std::ofstream stream("plot_column_Qeq");

  stream << "set encoding utf8\n";
  #if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
    stream << "set terminal pngcairo size ARG2,ARG3 enhanced font 'Arial,10'\n";
    stream << "set xlabel 'Adsorber position / [m]' font 'Arial,14'\n";
    stream << "set ylabel 'Concentration, {/Arial-Italic c}_i / [mol/kg]' offset 0.0,0 font 'Arial,14'\n";
    stream << "set key outside top center horizontal samplen 2.5 height 0.5 spacing 1.5 font 'Arial, 10'\n";
  #else
    stream << "set terminal pngcairo size ARG2,ARG3 enhanced font 'Helvetica,10'\n";
    stream << "set xlabel 'Adsorber position / [m]' font 'Helvetica,18'\n";
    stream << "set ylabel 'Concentration, {/Helvetica-Italic c}_i / [mol/kg]' offset 0.0,0 font 'Helvetica,18'\n";
    stream << "set key outside top center horizontal samplen 2.5 height 0.5 spacing 1.5 font 'Helvetica, 10'\n";
  #endif
  
  // colorscheme from book 'gnuplot in action', listing 12.7
  stream << "set linetype 1 pt 5 ps 1 lw 4 lc rgb '0xee0000'\n";
  stream << "set linetype 2 pt 7 ps 1 lw 4 lc rgb '0x008b00'\n";
  stream << "set linetype 3 pt 9 ps 1 lw 4 lc rgb '0x0000cd'\n";
  stream << "set linetype 4 pt 11 ps 1 lw 4 lc rgb '0xff3fb3'\n";
  stream << "set linetype 5 pt 13 ps 1 lw 4 lc rgb '0x00cdcd'\n";
  stream << "set linetype 6 pt 15 ps 1 lw 4 lc rgb '0xcd9b1d'\n";
  stream << "set linetype 7 pt  4 ps 1 lw 4 lc rgb '0x8968ed'\n";
  stream << "set linetype 8 pt  6 ps 1 lw 4 lc rgb '0x8b8b83'\n";
  stream << "set linetype 9 pt  8 ps 1 lw 4 lc rgb '0x00bb00'\n";
  stream << "set linetype 10 pt 10 ps 1 lw 4 lc rgb '0x1e90ff'\n";
  stream << "set linetype 11 pt 12 ps 1 lw 4 lc rgb '0x8b2500'\n";
  stream << "set linetype 12 pt 14 ps 1 lw 4 lc rgb '0x000000'\n";

  stream << "set bmargin 4\n";
  stream << "set key title '" << displayName << " {/:Italic T}=" << T_gas << " K, {/:Italic p_t}=" << p_total*1e-3 << " kPa'\n";
  stream << "stats 'column.data' nooutput\n";
  stream << "max = 0.0;\n";
  stream << "do for [i=5:STATS_columns:6] {\n";
  stream << "  stats 'column.data' us i nooutput\n";
  stream << "  if (max<STATS_max) {\n";
  stream << "    max=STATS_max\n";
  stream << "  }\n";
  stream << "}\n";
  stream << "stats 'column.data' us 1 nooutput\n";
  stream << "set xrange[0:STATS_max]\n";
  stream << "set yrange[0:1.1*max]\n";
  stream << "ev=int(ARG1)\n";
  stream << "do for [i=0:int((STATS_blocks-2)/ev)] {\n";
  stream << "  plot \\\n";
  for (size_t i = 0; i < Ncomp; i++)
  {
    stream << "    " << "'column.data'" << " us 1:" << std::to_string(5 + i * 6) << " index ev*i notitle " 
           << " with li lt " << i+1 << ",\\\n";
  }
  for (size_t i = 0; i < Ncomp; i++)
  {
    stream << "    " << "'column.data'" << " us 1:" << std::to_string(5 + i * 6) << " index ev*i title '"
           << components[i].name << " (y_i=" << components[i].Yi0 << ")'"
           << " with po lt " << i+1 << (i < Ncomp - 1 ? ",\\" : "") << "\n";
  }
  stream << "}\n";
}

void Breakthrough::createMovieScriptColumnP()
{
  #if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
    std::ofstream makeMovieStream("make_movie_P.bat");
  #else
    std::ofstream makeMovieStream("make_movie_P");
    makeMovieStream << "#!/bin/sh\n";
    makeMovieStream << "cd -- \"$(dirname \"$0\")\"\n";
  #endif
  makeMovieStream << movieScriptTemplate("P");

  #if (__cplusplus >= 201703L)
    std::filesystem::path path{"make_movie_P"};
    std::filesystem::permissions(path, std::filesystem::perms::owner_exec, std::filesystem::perm_options::add);
  #else
    chmod("make_movie_P", S_IRWXU);
  #endif

  std::ofstream stream("plot_column_P");

  stream << "set encoding utf8\n";
  #if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
    stream << "set terminal pngcairo size ARG2,ARG3 enhanced font 'Arial,10'\n";
    stream << "set xlabel 'Adsorber position / [m]' font 'Arial,14'\n";
    stream << "set ylabel 'Partial pressure, {/Arial-Italic p}_i / [Pa]' offset 0.0,0 font 'Arial,14'\n";
    stream << "set key outside top center horizontal samplen 2.5 height 0.5 spacing 1.5 font 'Arial, 10'\n";
  #else
    stream << "set terminal pngcairo size ARG2,ARG3 enhanced font 'Helvetica,10'\n";
    stream << "set xlabel 'Adsorber position / [m]' font 'Helvetica,18'\n";
    stream << "set ylabel 'Partial pressure, {/Helvetica-Italic p}_i / [Pa]' offset 0.0,0 font 'Helvetica,18'\n";
    stream << "set key outside top center horizontal samplen 2.5 height 0.5 spacing 1.5 font 'Helvetica, 10'\n";
  #endif

  // colorscheme from book 'gnuplot in action', listing 12.7
  stream << "set linetype 1 pt 5 ps 1 lw 4 lc rgb '0xee0000'\n";
  stream << "set linetype 2 pt 7 ps 1 lw 4 lc rgb '0x008b00'\n";
  stream << "set linetype 3 pt 9 ps 1 lw 4 lc rgb '0x0000cd'\n";
  stream << "set linetype 4 pt 11 ps 1 lw 4 lc rgb '0xff3fb3'\n";
  stream << "set linetype 5 pt 13 ps 1 lw 4 lc rgb '0x00cdcd'\n";
  stream << "set linetype 6 pt 15 ps 1 lw 4 lc rgb '0xcd9b1d'\n";
  stream << "set linetype 7 pt  4 ps 1 lw 4 lc rgb '0x8968ed'\n";
  stream << "set linetype 8 pt  6 ps 1 lw 4 lc rgb '0x8b8b83'\n";
  stream << "set linetype 9 pt  8 ps 1 lw 4 lc rgb '0x00bb00'\n";
  stream << "set linetype 10 pt 10 ps 1 lw 4 lc rgb '0x1e90ff'\n";
  stream << "set linetype 11 pt 12 ps 1 lw 4 lc rgb '0x8b2500'\n";
  stream << "set linetype 12 pt 14 ps 1 lw 4 lc rgb '0x000000'\n";

  stream << "set bmargin 4\n";
  stream << "set key title '" << displayName << " {/:Italic T}=" << T_gas << " K, {/:Italic p_t}=" << p_total*1e-3 << " kPa'\n";
  stream << "stats 'column.data' nooutput\n";
  stream << "max = 0.0;\n";
  stream << "do for [i=6:STATS_columns:6] {\n";
  stream << "  stats 'column.data' us i nooutput\n";
  stream << "  if (max<STATS_max) {\n";
  stream << "    max=STATS_max\n";
  stream << "  }\n";
  stream << "}\n";
  stream << "stats 'column.data' us 1 nooutput\n";
  stream << "set xrange[0:STATS_max]\n";
  stream << "set yrange[0:1.1*max]\n";
  stream << "ev=int(ARG1)\n";
  stream << "do for [i=0:int((STATS_blocks-2)/ev)] {\n";
  stream << "  plot \\\n";
  for (size_t i = 0; i < Ncomp; i++)
  {
    stream << "    " << "'column.data'" << " us 1:" << std::to_string(6 + i * 6) << " index ev*i notitle " 
           << " with li lt " << i+1 << ",\\\n";
  }
  for (size_t i = 0; i < Ncomp; i++)
  {
    stream << "    " << "'column.data'" << " us 1:" << std::to_string(6 + i * 6) << " index ev*i title '"
           << components[i].name << " (y_i=" << components[i].Yi0 << ")'"
           << " with po lt " << i+1 << (i < Ncomp - 1 ? ",\\" : "") << "\n";
  }
  stream << "}\n";
}

void Breakthrough::createMovieScriptColumnPnormalized()
{
  #if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
    std::ofstream makeMovieStream("make_movie_Pnorm.bat");
  #else
    std::ofstream makeMovieStream("make_movie_Pnorm");
    makeMovieStream << "#!/bin/sh\n";
    makeMovieStream << "cd -- \"$(dirname \"$0\")\"\n";
  #endif
  makeMovieStream << movieScriptTemplate("Pnorm");

  #if (__cplusplus >= 201703L)
    std::filesystem::path path{"make_movie_Pnorm"};
    std::filesystem::permissions(path, std::filesystem::perms::owner_exec, std::filesystem::perm_options::add);
  #else
    chmod("make_movie_Pnorm", S_IRWXU);
  #endif

  std::ofstream stream("plot_column_Pnorm");

  stream << "set encoding utf8\n";
  #if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
    stream << "set terminal pngcairo size ARG2,ARG3 enhanced font 'Arial,10'\n";
    stream << "set xlabel 'Adsorber position / [m]' font 'Arial,14'\n";
    stream << "set ylabel 'Partial pressure, {/Arial-Italic p}_i / [-]' offset 0.0,0 font 'Arial,14'\n";
    stream << "set key outside top center horizontal samplen 2.5 height 0.5 spacing 1.5 font 'Arial, 10'\n";
  #else
    stream << "set terminal pngcairo size ARG2,ARG3 enhanced font 'Helvetica,10'\n";
    stream << "set xlabel 'Adsorber position / [m]' font 'Helvetica,18'\n";
    stream << "set ylabel 'Partial pressure, {/Helvetica-Italic p}_i / [-]' offset 0.0,0 font 'Helvetica,18'\n";
    stream << "set key outside top center horizontal samplen 2.5 height 0.5 spacing 1.5 font 'Helvetica, 10'\n";
  #endif

  // colorscheme from book 'gnuplot in action', listing 12.7
  stream << "set linetype 1 pt 5 ps 1 lw 4 lc rgb '0xee0000'\n";
  stream << "set linetype 2 pt 7 ps 1 lw 4 lc rgb '0x008b00'\n";
  stream << "set linetype 3 pt 9 ps 1 lw 4 lc rgb '0x0000cd'\n";
  stream << "set linetype 4 pt 11 ps 1 lw 4 lc rgb '0xff3fb3'\n";
  stream << "set linetype 5 pt 13 ps 1 lw 4 lc rgb '0x00cdcd'\n";
  stream << "set linetype 6 pt 15 ps 1 lw 4 lc rgb '0xcd9b1d'\n";
  stream << "set linetype 7 pt  4 ps 1 lw 4 lc rgb '0x8968ed'\n";
  stream << "set linetype 8 pt  6 ps 1 lw 4 lc rgb '0x8b8b83'\n";
  stream << "set linetype 9 pt  8 ps 1 lw 4 lc rgb '0x00bb00'\n";
  stream << "set linetype 10 pt 10 ps 1 lw 4 lc rgb '0x1e90ff'\n";
  stream << "set linetype 11 pt 12 ps 1 lw 4 lc rgb '0x8b2500'\n";
  stream << "set linetype 12 pt 14 ps 1 lw 4 lc rgb '0x000000'\n";

  stream << "set bmargin 4\n";
  stream << "set key title '" << displayName << " {/:Italic T}=" << T_gas << " K, {/:Italic p_t}=" << p_total*1e-3 << " kPa'\n";
  stream << "stats 'column.data' nooutput\n";
  stream << "max = 0.0;\n";
  stream << "do for [i=7:STATS_columns:6] {\n";
  stream << "  stats 'column.data' us i nooutput\n";
  stream << "  if (max<STATS_max) {\n";
  stream << "    max=STATS_max\n";
  stream << "  }\n";
  stream << "}\n";
  stream << "stats 'column.data' us 1 nooutput\n";
  stream << "set xrange[0:STATS_max]\n";
  stream << "set yrange[0:1.1*max]\n";
  stream << "ev=int(ARG1)\n";
  stream << "do for [i=0:int((STATS_blocks-2)/ev)] {\n";
  stream << "  plot \\\n";
  for (size_t i = 0; i < Ncomp; i++)
  {
    stream << "    " << "'column.data'" << " us 1:" << std::to_string(7 + i * 6) << " index ev*i notitle " 
           << " with li lt " << i+1 << ",\\\n";
  }
  for (size_t i = 0; i < Ncomp; i++)
  {
    stream << "    " << "'column.data'" << " us 1:" << std::to_string(7 + i * 6) << " index ev*i title '"
           << components[i].name << " (y_i=" << components[i].Yi0 << ")'"
           << " with po lt " << i+1 << (i < Ncomp - 1 ? ",\\" : "") << "\n";
  }
  stream << "}\n";
}

void Breakthrough::createMovieScriptColumnDpdt()
{
  #if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
    std::ofstream makeMovieStream("make_movie_Dpdt.bat");
  #else
    std::ofstream makeMovieStream("make_movie_Dpdt");
    makeMovieStream << "#!/bin/sh\n";
    makeMovieStream << "cd -- \"$(dirname \"$0\")\"\n";
  #endif
  makeMovieStream << movieScriptTemplate("Dpdt");

  #if (__cplusplus >= 201703L)
    std::filesystem::path path{"make_movie_Dpdt"};
    std::filesystem::permissions(path, std::filesystem::perms::owner_exec, std::filesystem::perm_options::add);
  #else
    chmod("make_movie_Dpdt", S_IRWXU);
  #endif

  std::ofstream stream("plot_column_Dpdt");

  stream << "set encoding utf8\n";
  #if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
    stream << "set terminal pngcairo size ARG2,ARG3 enhanced font 'Arial,10'\n";
    stream << "set xlabel 'Adsorber position / [m]' font 'Arial,14'\n";
    stream << "set ylabel 'Pressure derivative, {/Arial-Italic dp_/dt} / [Pa/s]' offset 0.0,0 font 'Arial,14'\n";
    stream << "set key outside top center horizontal samplen 2.5 height 0.5 spacing 1.5 font 'Arial, 10'\n";
  #else
    stream << "set terminal pngcairo size ARG2,ARG3 enhanced font 'Helvetica,10'\n";
    stream << "set xlabel 'Adsorber position / [m]' font 'Helvetica,18'\n";
    stream << "set ylabel 'Pressure derivative, {/Helvetica-Italic dp_/dt} / [Pa/s]' offset 0.0,0 font 'Helvetica,18'\n";
    stream << "set key outside top center horizontal samplen 2.5 height 0.5 spacing 1.5 font 'Helvetica, 10'\n";
  #endif

  // colorscheme from book 'gnuplot in action', listing 12.7
  stream << "set linetype 1 pt 5 ps 1 lw 4 lc rgb '0xee0000'\n";
  stream << "set linetype 2 pt 7 ps 1 lw 4 lc rgb '0x008b00'\n";
  stream << "set linetype 3 pt 9 ps 1 lw 4 lc rgb '0x0000cd'\n";
  stream << "set linetype 4 pt 11 ps 1 lw 4 lc rgb '0xff3fb3'\n";
  stream << "set linetype 5 pt 13 ps 1 lw 4 lc rgb '0x00cdcd'\n";
  stream << "set linetype 6 pt 15 ps 1 lw 4 lc rgb '0xcd9b1d'\n";
  stream << "set linetype 7 pt  4 ps 1 lw 4 lc rgb '0x8968ed'\n";
  stream << "set linetype 8 pt  6 ps 1 lw 4 lc rgb '0x8b8b83'\n";
  stream << "set linetype 9 pt  8 ps 1 lw 4 lc rgb '0x00bb00'\n";
  stream << "set linetype 10 pt 10 ps 1 lw 4 lc rgb '0x1e90ff'\n";
  stream << "set linetype 11 pt 12 ps 1 lw 4 lc rgb '0x8b2500'\n";
  stream << "set linetype 12 pt 14 ps 1 lw 4 lc rgb '0x000000'\n";

  stream << "set bmargin 4\n";
  stream << "set key title '" << displayName << " {/:Italic T}=" << T_gas << " K, {/:Italic p_t}=" << p_total*1e-3 << " kPa'\n";
  stream << "stats 'column.data' nooutput\n";
  stream << "max = -1e10;\n";
  stream << "min = 1e10;\n";
  stream << "do for [i=8:STATS_columns:6] {\n";
  stream << "  stats 'column.data' us i nooutput\n";
  stream << "  if (STATS_max>max) {\n";
  stream << "    max=STATS_max\n";
  stream << "  }\n";
  stream << "  if (STATS_min<min) {\n";
  stream << "    min=STATS_min\n";
  stream << "  }\n";
  stream << "}\n";
  stream << "stats 'column.data' us 1 nooutput\n";
  stream << "set xrange[0:STATS_max]\n";
  stream << "set yrange[1.1*min:1.1*max]\n";
  stream << "ev=int(ARG1)\n";
  stream << "do for [i=0:int((STATS_blocks-2)/ev)] {\n";
  stream << "  plot \\\n";
  for (size_t i = 0; i < Ncomp; i++)
  {
    stream << "    " << "'column.data'" << " us 1:" << std::to_string(8 + i * 6) << " index ev*i notitle " 
           << " with li lt " << i+1 << ",\\\n";
  }
  for (size_t i = 0; i < Ncomp; i++)
  {
    stream << "    " << "'column.data'" << " us 1:" << std::to_string(8 + i * 6) << " index ev*i title '"
           << components[i].name << " (y_i=" << components[i].Yi0 << ")'"
           << " with po lt " << i+1 << (i < Ncomp - 1 ? ",\\" : "") << "\n";
  }
  stream << "}\n";
}

void Breakthrough::createMovieScriptColumnDqdt()
{
  #if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
    std::ofstream makeMovieStream("make_movie_Dqdt.bat");
  #else
    std::ofstream makeMovieStream("make_movie_Dqdt");
    makeMovieStream << "#!/bin/sh\n";
    makeMovieStream << "cd -- \"$(dirname \"$0\")\"\n";
  #endif
  makeMovieStream << movieScriptTemplate("Dqdt");

  #if (__cplusplus >= 201703L)
    std::filesystem::path path{"make_movie_Dqdt"};
    std::filesystem::permissions(path, std::filesystem::perms::owner_exec, std::filesystem::perm_options::add);
  #else
    chmod("make_movie_Dqdt", S_IRWXU);
  #endif

  std::ofstream stream("plot_column_Dqdt");

  stream << "set encoding utf8\n";
  #if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
    stream << "set terminal pngcairo size ARG2,ARG3 enhanced font 'Arial,10'\n";
    stream << "set xlabel 'Adsorber position / [m]' font 'Arial,14'\n";
    stream << "set ylabel 'Loading derivative, {/Arial-Italic dq_i/dt} / [mol/kg/s]' offset 0.0,0 font 'Arial,14'\n";
    stream << "set key outside top center horizontal samplen 2.5 height 0.5 spacing 1.5 font 'Arial, 10'\n";
  #else
    stream << "set terminal pngcairo size ARG2,ARG3 enhanced font 'Helvetica,10'\n";
    stream << "set xlabel 'Adsorber position / [m]' font 'Helvetica,18'\n";
    stream << "set ylabel 'Loading derivative, {/Helvetica-Italic dq_i/dt} / [mol/kg/s]' offset 0.0,0 font 'Helvetica,18'\n";
    stream << "set key outside top center horizontal samplen 2.5 height 0.5 spacing 1.5 font 'Helvetica, 10'\n";
  #endif

  // colorscheme from book 'gnuplot in action', listing 12.7
  stream << "set linetype 1 pt 5 ps 1 lw 4 lc rgb '0xee0000'\n";
  stream << "set linetype 2 pt 7 ps 1 lw 4 lc rgb '0x008b00'\n";
  stream << "set linetype 3 pt 9 ps 1 lw 4 lc rgb '0x0000cd'\n";
  stream << "set linetype 4 pt 11 ps 1 lw 4 lc rgb '0xff3fb3'\n";
  stream << "set linetype 5 pt 13 ps 1 lw 4 lc rgb '0x00cdcd'\n";
  stream << "set linetype 6 pt 15 ps 1 lw 4 lc rgb '0xcd9b1d'\n";
  stream << "set linetype 7 pt  4 ps 1 lw 4 lc rgb '0x8968ed'\n";
  stream << "set linetype 8 pt  6 ps 1 lw 4 lc rgb '0x8b8b83'\n";
  stream << "set linetype 9 pt  8 ps 1 lw 4 lc rgb '0x00bb00'\n";
  stream << "set linetype 10 pt 10 ps 1 lw 4 lc rgb '0x1e90ff'\n";
  stream << "set linetype 11 pt 12 ps 1 lw 4 lc rgb '0x8b2500'\n";
  stream << "set linetype 12 pt 14 ps 1 lw 4 lc rgb '0x000000'\n";

  stream << "set bmargin 4\n";
  stream << "set key title '" << displayName << " {/:Italic T}=" << T_gas << " K, {/:Italic p_t}=" << p_total*1e-3 << " kPa'\n";
  stream << "stats 'column.data' nooutput\n";
  stream << "max = -1e10;\n";
  stream << "min = 1e10;\n";
  stream << "min = 10000000000000.0;\n";
  stream << "do for [i=9:STATS_columns:6] {\n";
  stream << "  stats 'column.data' us i nooutput\n";
  stream << "  if (STATS_max>max) {\n";
  stream << "    max=STATS_max\n";
  stream << "  }\n";
  stream << "  if (STATS_min<min) {\n";
  stream << "    min=STATS_min\n";
  stream << "  }\n";
  stream << "}\n";
  stream << "stats 'column.data' us 1 nooutput\n";
  stream << "set xrange[0:STATS_max]\n";
  stream << "set yrange[1.1*min:1.1*max]\n";
  stream << "ev=int(ARG1)\n";
  stream << "do for [i=0:int((STATS_blocks-2)/ev)] {\n";
  stream << "  plot \\\n";
  for (size_t i = 0; i < Ncomp; i++)
  {
    stream << "    " << "'column.data'" << " us 1:" << std::to_string(9 + i * 6) << " index ev*i notitle " 
           << " with li lt " << i+1 << ",\\\n";
  }
  for (size_t i = 0; i < Ncomp; i++)
  {
    stream << "    " << "'column.data'" << " us 1:" << std::to_string(9 + i * 6) << " index ev*i title '"
           << components[i].name << " (y_i=" << components[i].Yi0 << ")'"
           << " with po lt " << i+1 << (i < Ncomp - 1 ? ",\\" : "") << "\n";
  }
  stream << "}\n";
}

#ifdef PYBUILD
py::array_t<double> Breakthrough::compute()
{
  size_t colsize = 6 * Ncomp + 5;
  std::vector<std::vector<std::vector<double>>> brk;

  // loop can quit early if autoSteps
  for (size_t step = 0; (step < Nsteps || autoSteps); ++step)
  {
    // check for error from python side (keyboard interrupt)
    if (PyErr_CheckSignals() != 0)
    {
      throw py::error_already_set();
    }

    computeStep(step);
    double t = static_cast<double>(step) * dt;
    if (step % writeEvery == 0)
    {
      std::vector<std::vector<double>> t_brk(Ngrid + 1, std::vector<double>(colsize));
      for (size_t i = 0; i < Ngrid + 1; ++i)
      {
        t_brk[i][0] = t * v_in / L;
        t_brk[i][1] = t / 60.0;
        t_brk[i][2] = static_cast<double>(i) * dx;
        t_brk[i][3] = V[i];
        t_brk[i][4] = Pt[i];

        for (size_t j = 0; j < Ncomp; ++j)
        {
          t_brk[i][5 + 6 * j] = Q[i * Ncomp + j];
          t_brk[i][6 + 6 * j] = Qeq[i * Ncomp + j];
          t_brk[i][7 + 6 * j] = P[i * Ncomp + j];
          t_brk[i][8 + 6 * j] = P[i * Ncomp + j] / (Pt[i] * components[j].Yi0);
          t_brk[i][9 + 6 * j] = Dpdt[i * Ncomp + j];
          t_brk[i][10 + 6 * j] = Dqdt[i * Ncomp + j];
        }
      }
      brk.push_back(t_brk);
    }
    if (step % printEvery == 0)
    {
      std::cout << "Timestep " + std::to_string(step) + ", time: " + std::to_string(t) + " [s]\n";
      std::cout << "    Average number of mixture-prediction steps: " +
                       std::to_string(static_cast<double>(iastPerformance.first) /
                                      static_cast<double>(iastPerformance.second))
                << "\n";
    }
  }
  std::cout << "Final timestep " + std::to_string(Nsteps) +
                   ", time: " + std::to_string(dt * static_cast<double>(Nsteps)) + " [s]\n";

  std::vector<double> buffer;
  buffer.reserve(brk.size() * (Ngrid + 1) * colsize);
  for (const auto &vec1 : brk)
  {
    for (const auto &vec2 : vec1)
    {
      buffer.insert(buffer.end(), vec2.begin(), vec2.end());
    }
  }
  std::array<size_t, 3> shape{{brk.size(), Ngrid + 1, colsize}};
  py::array_t<double> py_breakthrough(shape, buffer.data());
  py_breakthrough.resize(shape);

  return py_breakthrough;
}
#endif  // PYBUILD
