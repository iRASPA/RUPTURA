import ruptura
import matplotlib.pyplot as plt

components = ruptura.Components([{
    "name": "Helium",
    "isCarrierGas": True,
    "gasPhaseMolFraction": 0.0
}, {
    "name": "CO2",
    "gasPhaseMolFraction": 0.5,
    "isotherms": [["Langmuir", 4.4, 2.91e-04], ["Langmuir", 10.0, 6.96e-07]]
}, {
    "name": "C3H8",
    "gasPhaseMolFraction": 0.5,
    "isotherms": [["Langmuir", 4.97, 9.30e-09], ["Langmuir", 2.99, 2.79e-04]]
}])

mix = ruptura.MixturePrediction(components=components,
                                displayName="MOR",
                                temperature=300.0,
                                pressureStart=1e2,
                                pressureEnd=1e8,
                                numberOfPressurePoints=100,
                                pressureScale="log",
                                predictionMethod="SIAST")

data = mix.compute()
fig, ax = plt.subplots(1, 2, figsize=(8, 6))
mix.plot(ax[0], "pure")
mix.plot(ax[1], "mixture")
plt.show()
