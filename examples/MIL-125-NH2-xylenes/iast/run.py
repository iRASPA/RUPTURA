import ruptura
import matplotlib.pyplot as plt

components = ruptura.Components([{
    "name": "o-xylene",
    "gasPhaseMolFraction": 0.25,
    "isotherms": [["Langmuir-Freundlich", 1.8, 1.59e-5, 0.62], ["Langmuir-Freundlich", 2.9, 1.5e-4, 1.06]]
}, {
    "name": "m-xylene",
    "gasPhaseMolFraction": 0.25,
    "isotherms": [["Langmuir-Freundlich", 1.1, 1.69e-5, 0.68], ["Langmuir-Freundlich", 3.0, 9.98e-5, 1.1]]
}, {
    "name": "p-xylene",
    "gasPhaseMolFraction": 0.25,
    "isotherms": [["Langmuir-Freundlich", 1.2, 3.61e-6, 0.87], ["Langmuir-Freundlich", 3.3, 8.33e-5, 1.15]]
}, {
    "name": "ethylbenzene",
    "gasPhaseMolFraction": 0.25,
    "isotherms": [["Langmuir-Freundlich", 1.3, 2.61e-6, 0.9], ["Langmuir-Freundlich", 3.1, 8.02e-5, 1.1]]
}])

mix = ruptura.MixturePrediction(components=components,
                                displayName="MIL-125(Ti)-NH_2",
                                temperature=443.0,
                                pressureStart=1e2,
                                pressureEnd=1e7,
                                numberOfPressurePoints=100,
                                pressureScale="log")

data = mix.compute()
fig, ax = plt.subplots(1, 2, figsize=(8, 6))
mix.plot(ax[0], "pure")
mix.plot(ax[1], "mixture")
plt.show()
