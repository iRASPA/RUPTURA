import ruptura
import matplotlib.pyplot as plt

components = ruptura.Components([{
    "MoleculeName": "CO2",
    "GasPhaseMolFraction": 0.5,
    "isotherms": [["Bingel-Walton", 2.77, 9.695e-6, -10.e-6]]
}, {
    "MoleculeName": "O2",
    "GasPhaseMolFraction": 0.5,
    "isotherms": [["Bingel-Walton", 1.85, 8.094e-11, 2.214e-6]]
}])

mix = ruptura.MixturePrediction(components=components,
                                DisplayName="Cu(dhbc)2(4,4'bpy)",
                                Temperature=298.0,
                                PressureStart=1e1,
                                PressureEnd=8e6,
                                NumberOfPressurePoints=50,
                                PressureScale="linear")

data = mix.compute()
fig, ax = plt.subplots(1, 2, figsize=(8, 6))
mix.plot(ax[0], "pure")
mix.plot(ax[1], "mixture")
plt.show()
