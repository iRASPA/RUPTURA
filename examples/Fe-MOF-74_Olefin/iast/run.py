import ruptura
import matplotlib.pyplot as plt

components = ruptura.Components([{
    "name":
        "1-3-butadiene",
    "gasPhaseMolFraction":
        0.25,
    "isotherms": [["Langmuir-Freundlich", 1.77359, 0.964548, 0.891275],
                  ["Langmuir-Freundlich", 4.71979, 4.193e+06, 1.76902],
                  ["Langmuir-Freundlich", 0.874237, 2857.35, 1.65435]]
}, {
    "name":
        "1-butene",
    "gasPhaseMolFraction":
        0.25,
    "isotherms": [["Langmuir-Freundlich", 4.12029, 260699.0, 2.19499],
                  ["Langmuir-Freundlich", 2.33644, 0.953314, 0.431999]]
}, {
    "name":
        "2-cis-butene",
    "gasPhaseMolFraction":
        0.25,
    "isotherms": [["Langmuir-Freundlich", 5.25517, 3.14573e+06, 2.10887],
                  ["Langmuir-Freundlich", 1.28055, 7.58487, 0.75982]]
}, {
    "name":
        "2-trans-butene",
    "gasPhaseMolFraction":
        0.25,
    "isotherms": [["Langmuir-Freundlich", 3.81848, 2.12583e+08, 3.11845],
                  ["Langmuir-Freundlich", 2.95468, 3.58644, 0.531051]]
}])

mix = ruptura.MixturePrediction(components=components,
                                displayName="Fe-MOF-74",
                                temperature=298.0,
                                pressureStart=1e-6,
                                pressureEnd=1e4,
                                numberOfPressurePoints=512,
                                pressureScale="log")

data = mix.compute()
fig, ax = plt.subplots(1, 2, figsize=(8, 6))
mix.plot(ax[0], "pure")
mix.plot(ax[1], "mixture")
plt.show()
