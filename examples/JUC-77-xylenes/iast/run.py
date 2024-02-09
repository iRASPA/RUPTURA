import ruptura
import matplotlib.pyplot as plt

components = ruptura.Components([{
    "MoleculeName": "o-xylene",
    "GasPhaseMolFraction": 0.25,
    "isotherms": [["Langmuir-Freundlich", 0.4, 2.4e-3, 0.8], ["Langmuir", 0.9, 2.01e-2]]
}, {
    "MoleculeName": "m-xylene",
    "GasPhaseMolFraction": 0.25,
    "isotherms": [["Langmuir-Freundlich", 0.53, 1.24e-2, 0.822], ["Langmuir-Freundlich", 0.77, 3.1e-2, 1.1]]
}, {
    "MoleculeName": "p-xylene",
    "GasPhaseMolFraction": 0.25,
    "isotherms": [["Langmuir-Freundlich", 0.6, 1.24e-2, 0.9], ["Langmuir-Freundlich", 0.7, 3.14e-2, 1.05]]
}, {
    "MoleculeName": "ethylbenzene",
    "GasPhaseMolFraction": 0.25,
    "isotherms": [["Langmuir-Freundlich", 0.4, 7.02e-4, 0.9], ["Langmuir-Freundlich", 0.9, 2.98e-3, 1.2]]
}])

mix = ruptura.MixturePrediction(components=components,
                                DisplayName="JUC-77",
                                Temperature=443.0,
                                PressureStart=1e0,
                                PressureEnd=1e6,
                                NumberOfPressurePoints=100,
                                PressureScale="log")

data = mix.compute()
fig, ax = plt.subplots(1, 2, figsize=(8, 6))
mix.plot(ax[0], "pure")
mix.plot(ax[1], "mixture")
plt.show()
