import ruptura
import matplotlib.pyplot as plt

components = ruptura.Components([{
    "MoleculeName": "Helium",
    "CarrierGas": True,
    "GasPhaseMolFraction": 0.0
}, {
    "MoleculeName": "CO2",
    "GasPhaseMolFraction": 0.5,
    "isotherms": [["Langmuir", 4.4, 2.91e-04], ["Langmuir", 10.0, 6.96e-07]]
}, {
    "MoleculeName": "C3H8",
    "GasPhaseMolFraction": 0.5,
    "isotherms": [["Langmuir", 4.97, 9.30e-09], ["Langmuir", 2.99, 2.79e-04]]
}])

mix = ruptura.MixturePrediction(components=components,
                                DisplayName="MOR",
                                Temperature=300.0,
                                PressureStart=1e2,
                                PressureEnd=1e8,
                                NumberOfPressurePoints=100,
                                PressureScale="log",
                                MixturePredictionMethod="SIAST")

data = mix.compute()
fig, ax = plt.subplots(1, 2, figsize=(8, 6))
mix.plot(ax[0], "pure")
mix.plot(ax[1], "mixture")
plt.show()
