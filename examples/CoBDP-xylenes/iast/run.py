import ruptura
import matplotlib.pyplot as plt

components = ruptura.Components([{
    "MoleculeName": "o-xylene",
    "GasPhaseMolFraction": 0.25,
    "isotherms": [["Langmuir-Freundlich", 4.8, 1.28e-5, 1.4], ["Langmuir-Freundlich", 1.5, 1.6e-4, 0.7]]
}, {
    "MoleculeName": "m-xylene",
    "GasPhaseMolFraction": 0.25,
    "isotherms": [["Langmuir-Freundlich", 4.8, 3.19e-5, 1.27], ["Langmuir-Freundlich", 1.2, 4.72e-5, 0.7]]
}, {
    "MoleculeName": "p-xylene",
    "GasPhaseMolFraction": 0.25,
    "isotherms": [["Langmuir-Freundlich", 4.5, 2.3e-6, 1.7], ["Langmuir-Freundlich", 1.6, 1.46e-4, 0.7]]
}, {
    "MoleculeName": "ethylbenzene",
    "GasPhaseMolFraction": 0.25,
    "isotherms": [["Langmuir-Freundlich", 4.5, 7.38e-6, 1.5], ["Langmuir-Freundlich", 1.4, 1.6e-4, 0.7]]
}])

mix = ruptura.MixturePrediction(components=components,
                                DisplayName="CoBDP",
                                Temperature=443.0,
                                PressureStart=1e2,
                                PressureEnd=1e6,
                                NumberOfPressurePoints=100,
                                PressureScale="log")

data = mix.compute()
fig, ax = plt.subplots(1, 2, figsize=(8, 6))
mix.plot(ax[0], "pure")
mix.plot(ax[1], "mixture")
plt.show()
