import ruptura
import matplotlib.pyplot as plt

components = ruptura.Components([{
    "MoleculeName": "o-xylene",
    "GasPhaseMolFraction": 0.25,
    "isotherms": [["Langmuir", 4.0, 20.0]]
}, {
    "MoleculeName": "m-xylene",
    "GasPhaseMolFraction": 0.25,
    "isotherms": [["Langmuir", 4.0, 87.5], ["Langmuir", 1.8, 2.3e-5]]
}, {
    "MoleculeName": "p-xylene",
    "GasPhaseMolFraction": 0.25,
    "isotherms": [["Langmuir-Freundlich", 4.0, 29.6, 1.04], ["Langmuir-Freundlich", 2.2, 1.05e-3, 0.75]]
}, {
    "MoleculeName": "ethylbenzene",
    "GasPhaseMolFraction": 0.25,
    "isotherms": [["Langmuir", 4.0, 103.0], ["Langmuir", 2.0, 2.99e-7]]
}])

mix = ruptura.MixturePrediction(components=components,
                                DisplayName="MFI",
                                Temperature=443.0,
                                PressureStart=1e-4,
                                PressureEnd=1e7,
                                NumberOfPressurePoints=100,
                                PressureScale="log")

data = mix.compute()
fig, ax = plt.subplots(1, 2, figsize=(8, 6))
mix.plot(ax[0], "pure")
mix.plot(ax[1], "mixture")
plt.show()
