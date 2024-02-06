import ruptura
import matplotlib.pyplot as plt

components = ruptura.Components([{
    "name": "o-xylene",
    "gasPhaseMolFraction": 0.25,
    "isotherms": [["Langmuir-Freundlich", 1.7, 7.74e-1, 0.4], ["Langmuir-Freundlich", 1.8, 5.53e6, 1.1]]
}, {
    "name": "m-xylene",
    "gasPhaseMolFraction": 0.25,
    "isotherms": [["Langmuir-Freundlich", 1.4, 1.22e1, 0.57], ["Langmuir-Freundlich", 1.8, 2.68e6, 1.07]]
}, {
    "name": "p-xylene",
    "gasPhaseMolFraction": 0.25,
    "isotherms": [["Langmuir-Freundlich", 1.6, 1.07e1, 0.54], ["Langmuir-Freundlich", 1.8, 4.97e7, 1.2]]
}, {
    "name": "ethylbenzene",
    "gasPhaseMolFraction": 0.25,
    "isotherms": [["Langmuir-Freundlich", 1.5, 7.53, 0.73], ["Langmuir", 1.3, 3.76e4]]
}])

mix = ruptura.MixturePrediction(components=components,
                                displayName="MAF-X8",
                                temperature=443.0,
                                pressureStart=1e-8,
                                pressureEnd=1e4,
                                numberOfPressurePoints=100,
                                pressureScale="log")

data = mix.compute()
fig, ax = plt.subplots(1, 2, figsize=(8, 6))
mix.plot(ax[0], "pure")
mix.plot(ax[1], "mixture")
plt.show()
