import ruptura
import matplotlib.pyplot as plt

components = ruptura.Components([{
    "name": "o-xylene",
    "gasPhaseMolFraction": 0.25,
    "isotherms": [["Langmuir-Freundlich", 2.8, 5.03e-3, 1.08], ["Langmuir-Freundlich", 0.8, 3.78e-4, 0.8]]
}, {
    "name": "m-xylene",
    "gasPhaseMolFraction": 0.25,
    "isotherms": [["Langmuir-Freundlich", 2.6, 2.48e-3, 1.1], ["Langmuir-Freundlich", 0.8, 3.86e-5, 0.9]]
}, {
    "name": "p-xylene",
    "gasPhaseMolFraction": 0.25,
    "isotherms": [["Langmuir-Freundlich", 2.8, 3.22e-3, 1.07], ["Langmuir-Freundlich", 0.8, 3.23e-5, 0.9]]
}, {
    "name": "ethylbenzene",
    "gasPhaseMolFraction": 0.25,
    "isotherms": [["Langmuir-Freundlich", 0.9, 3.94e-4, 0.62], ["Langmuir-Freundlich", 2.2, 4.73e-3, 1.05]]
}])

mix = ruptura.MixturePrediction(components=components,
                                displayName="MIL-47",
                                temperature=443.0,
                                pressureStart=1e0,
                                pressureEnd=1e7,
                                numberOfPressurePoints=100,
                                pressureScale="log")

data = mix.compute()
fig, ax = plt.subplots(1, 2, figsize=(8, 6))
mix.plot(ax[0], "pure")
mix.plot(ax[1], "mixture")
plt.show()
