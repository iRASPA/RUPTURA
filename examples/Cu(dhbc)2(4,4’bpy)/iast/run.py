import ruptura
import matplotlib.pyplot as plt

components = ruptura.Components([{
    "name": "CO2",
    "gasPhaseMolFraction": 0.5,
    "isotherms": [["Bingel-Walton", 2.77, 9.695e-6, -10.e-6]]
}, {
    "name": "O2",
    "gasPhaseMolFraction": 0.5,
    "isotherms": [["Bingel-Walton", 1.85, 8.094e-11, 2.214e-6]]
}])

mix = ruptura.MixturePrediction(components=components,
                                displayName="Cu(dhbc)2(4,4'bpy)",
                                temperature=298.0,
                                pressureStart=1e1,
                                pressureEnd=8e6,
                                numberOfPressurePoints=50,
                                pressureScale="linear")

data = mix.compute()
fig, ax = plt.subplots(1, 2, figsize=(8, 6))
mix.plot(ax[0], "pure")
mix.plot(ax[1], "mixture")
plt.show()
