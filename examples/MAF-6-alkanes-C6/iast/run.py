import ruptura
import matplotlib.pyplot as plt

components = ruptura.Components([{
    "name":
        "22DMB",
    "gasPhaseMolFraction":
        0.5,
    "isotherms": [["Langmuir-Freundlich", 2.88752, 1.14774e-07, 2.25341],
                  ["Langmuir-Freundlich", 0.938131, 0.0266853, 1.13388]]
}, {
    "name":
        "23DMB",
    "gasPhaseMolFraction":
        0.5,
    "isotherms": [["Langmuir-Freundlich", 3.13648, 1.37996e-08, 3.11688],
                  ["Langmuir-Freundlich", 0.00385344, 4.51394, 10.1938],
                  ["Langmuir-Freundlich", 0.875397, 0.105484, 1.24258]]
}])

mix = ruptura.MixturePrediction(components=components,
                                displayName="MAF-6",
                                temperature=298.0,
                                pressureStart=1e-1,
                                pressureEnd=1e6,
                                numberOfPressurePoints=256,
                                pressureScale="log")

data = mix.compute()
fig, ax = plt.subplots(1, 2, figsize=(8, 6))
mix.plot(ax[0], "pure")
mix.plot(ax[1], "mixture")
plt.show()
