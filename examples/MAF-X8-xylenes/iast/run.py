import ruptura
import matplotlib.pyplot as plt

components = ruptura.Components([{
    "MoleculeName": "o-xylene",
    "GasPhaseMolFraction": 0.25,
    "isotherms": [["Langmuir-Freundlich", 1.7, 7.74e-1, 0.4], ["Langmuir-Freundlich", 1.8, 5.53e6, 1.1]]
}, {
    "MoleculeName": "m-xylene",
    "GasPhaseMolFraction": 0.25,
    "isotherms": [["Langmuir-Freundlich", 1.4, 1.22e1, 0.57], ["Langmuir-Freundlich", 1.8, 2.68e6, 1.07]]
}, {
    "MoleculeName": "p-xylene",
    "GasPhaseMolFraction": 0.25,
    "isotherms": [["Langmuir-Freundlich", 1.6, 1.07e1, 0.54], ["Langmuir-Freundlich", 1.8, 4.97e7, 1.2]]
}, {
    "MoleculeName": "ethylbenzene",
    "GasPhaseMolFraction": 0.25,
    "isotherms": [["Langmuir-Freundlich", 1.5, 7.53, 0.73], ["Langmuir", 1.3, 3.76e4]]
}])

mix = ruptura.MixturePrediction(components=components,
                                DisplayName="MAF-X8",
                                Temperature=443.0,
                                PressureStart=1e-8,
                                PressureEnd=1e4,
                                NumberOfPressurePoints=100,
                                PressureScale="log")

data = mix.compute()
fig, ax = plt.subplots(1, 2, figsize=(8, 6))
mix.plot(ax[0], "pure")
mix.plot(ax[1], "mixture")
plt.show()
