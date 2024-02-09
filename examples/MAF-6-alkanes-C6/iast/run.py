import ruptura
import matplotlib.pyplot as plt

components = ruptura.Components([{
    "MoleculeName":
        "22DMB",
    "GasPhaseMolFraction":
        0.5,
    "isotherms": [["Langmuir-Freundlich", 2.88752, 1.14774e-07, 2.25341],
                  ["Langmuir-Freundlich", 0.938131, 0.0266853, 1.13388]]
}, {
    "MoleculeName":
        "23DMB",
    "GasPhaseMolFraction":
        0.5,
    "isotherms": [["Langmuir-Freundlich", 3.13648, 1.37996e-08, 3.11688],
                  ["Langmuir-Freundlich", 0.00385344, 4.51394, 10.1938],
                  ["Langmuir-Freundlich", 0.875397, 0.105484, 1.24258]]
}])

mix = ruptura.MixturePrediction(components=components,
                                DisplayName="MAF-6",
                                Temperature=298.0,
                                PressureStart=1e-1,
                                PressureEnd=1e6,
                                NumberOfPressurePoints=256,
                                PressureScale="log")

data = mix.compute()
fig, ax = plt.subplots(1, 2, figsize=(8, 6))
mix.plot(ax[0], "pure")
mix.plot(ax[1], "mixture")
plt.show()
