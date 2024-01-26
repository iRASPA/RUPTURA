#include <exception>

#include "special_functions.h"

#include "inputreader.h"
#include "breakthrough.h"
#include "mixture_prediction.h"
#include "fitting.h"

int main(void)
{
  try 
  {
    InputReader reader("simulation.input");

    switch(reader.simulationType)
    {
      case InputReader::SimulationType::Breakthrough:
      default:
      {
        Breakthrough breakthrough(reader);

        std::cout << breakthrough.repr();
        breakthrough.initialize();
        breakthrough.createPlotScript();
        breakthrough.createMovieScripts();
        breakthrough.run();
        break;
      }
      case InputReader::SimulationType::MixturePrediction:
      {
        MixturePrediction mixture(reader);

        std::cout << mixture.repr();
        mixture.run();
        mixture.createPureComponentsPlotScript();
        mixture.createMixturePlotScript();
        mixture.createMixtureAdsorbedMolFractionPlotScript();
        mixture.createPlotScript();
        std::cout << mixture.repr();
        break;
      }
      case InputReader::SimulationType::Fitting:
      {
        Fitting fitting(reader);

        fitting.run();
        break;
      }
    }
  }
  catch (std::exception const& e)
  {
    std::cerr << e.what();
    exit(-1);
  }

  return 0;
}
