import ruptura
import matplotlib.pyplot as plt

components = ruptura.Components([{
    "MoleculeName": "nC6",
    "GasPhaseMolFraction": 0.2,
    "isotherms": [["Langmuir-Freundlich", 1.1, 2.14e-4, 0.6], ["Langmuir-Freundlich", 4.8, 2.39e-5, 1.35]]
}, {
    "MoleculeName": "2MP",
    "GasPhaseMolFraction": 0.2,
    "isotherms": [["Langmuir-Freundlich", 1.4, 1.68e-4, 0.6], ["Langmuir-Freundlich", 4.8, 3.3e-5, 1.26]]
}, {
    "MoleculeName": "3MP",
    "GasPhaseMolFraction": 0.2,
    "isotherms": [["Langmuir-Freundlich", 1.75, 1.14e-4, 0.62], ["Langmuir-Freundlich", 4.8, 3.22e-5, 1.24]]
}, {
    "MoleculeName": "23DMB",
    "GasPhaseMolFraction": 0.2,
    "isotherms": [["Langmuir-Freundlich", 1.5, 8.31e-5, 0.63], ["Langmuir-Freundlich", 4.8, 5.46e-5, 1.13]]
}, {
    "MoleculeName": "22DMB",
    "GasPhaseMolFraction": 0.2,
    "isotherms": [["Langmuir-Freundlich", 1.74, 5.65e-6, 0.72], ["Langmuir-Freundlich", 4.8, 4.3e-5, 1.04]]
}])

mix = ruptura.MixturePrediction(components=components,
                                DisplayName="CoBDP",
                                Temperature=443.0,
                                PressureStart=1e2,
                                PressureEnd=1e6,
                                NumberOfPressurePoints=100,
                                PressureScale="log")

data = mix.compute()
fig, ax = plt.subplots(1, 2, figsize=(8, 6))
mix.plot(ax[0], "pure")
mix.plot(ax[1], "mixture")
plt.show()
