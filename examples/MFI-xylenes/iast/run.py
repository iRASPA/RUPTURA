import ruptura
import matplotlib.pyplot as plt

components = ruptura.Components([{
    "name": "o-xylene",
    "gasPhaseMolFraction": 0.25,
    "isotherms": [["Langmuir", 4.0, 20.0]]
}, {
    "name": "m-xylene",
    "gasPhaseMolFraction": 0.25,
    "isotherms": [["Langmuir", 4.0, 87.5], ["Langmuir", 1.8, 2.3e-5]]
}, {
    "name": "p-xylene",
    "gasPhaseMolFraction": 0.25,
    "isotherms": [["Langmuir-Freundlich", 4.0, 29.6, 1.04], ["Langmuir-Freundlich", 2.2, 1.05e-3, 0.75]]
}, {
    "name": "ethylbenzene",
    "gasPhaseMolFraction": 0.25,
    "isotherms": [["Langmuir", 4.0, 103.0], ["Langmuir", 2.0, 2.99e-7]]
}])

mix = ruptura.MixturePrediction(components=components,
                                displayName="MFI",
                                temperature=443.0,
                                pressureStart=1e-4,
                                pressureEnd=1e7,
                                numberOfPressurePoints=100,
                                pressureScale="log")

data = mix.compute()
fig, ax = plt.subplots(1, 2, figsize=(8, 6))
mix.plot(ax[0], "pure")
mix.plot(ax[1], "mixture")
plt.show()
