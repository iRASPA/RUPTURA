import ruptura
import matplotlib.pyplot as plt

components = ruptura.Components([{
    "MoleculeName": "o-xylene",
    "GasPhaseMolFraction": 0.25,
    "isotherms": [["Langmuir-Freundlich", 2.8, 5.03e-3, 1.08], ["Langmuir-Freundlich", 0.8, 3.78e-4, 0.8]]
}, {
    "MoleculeName": "m-xylene",
    "GasPhaseMolFraction": 0.25,
    "isotherms": [["Langmuir-Freundlich", 2.6, 2.48e-3, 1.1], ["Langmuir-Freundlich", 0.8, 3.86e-5, 0.9]]
}, {
    "MoleculeName": "p-xylene",
    "GasPhaseMolFraction": 0.25,
    "isotherms": [["Langmuir-Freundlich", 2.8, 3.22e-3, 1.07], ["Langmuir-Freundlich", 0.8, 3.23e-5, 0.9]]
}, {
    "MoleculeName": "ethylbenzene",
    "GasPhaseMolFraction": 0.25,
    "isotherms": [["Langmuir-Freundlich", 0.9, 3.94e-4, 0.62], ["Langmuir-Freundlich", 2.2, 4.73e-3, 1.05]]
}])

mix = ruptura.MixturePrediction(components=components,
                                DisplayName="MIL-47",
                                Temperature=443.0,
                                PressureStart=1e0,
                                PressureEnd=1e7,
                                NumberOfPressurePoints=100,
                                PressureScale="log")

data = mix.compute()
fig, ax = plt.subplots(1, 2, figsize=(8, 6))
mix.plot(ax[0], "pure")
mix.plot(ax[1], "mixture")
plt.show()
