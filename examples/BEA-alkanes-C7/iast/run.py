import ruptura
import matplotlib.pyplot as plt

components = ruptura.Components([{
    "MoleculeName": "Helium",
    "CarrierGas": True,
    "GasPhaseMolFraction": 0.93,
    "isotherms": []
}, {
    "MoleculeName": "nC7",
    "CarrierGas": False,
    "GasPhaseMolFraction": 0.01,
    "isotherms": [["Langmuir", 1.09984, 6.55857e-5], ["Langmuir", 0.19466, 8.90731e-07]]
}, {
    "MoleculeName": "C6m2",
    "CarrierGas": False,
    "GasPhaseMolFraction": 0.01,
    "isotherms": [["Langmuir", 1.22228, 3.90895e-05], ["Langmuir", 0.481726, 9.64046e-08]]
}, {
    "MoleculeName": "C6m3",
    "CarrierGas": False,
    "GasPhaseMolFraction": 0.01,
    "isotherms": [["Langmuir", 1.25189, 2.6174e-05], ["Langmuir", 0.518406, 7.99897e-08]]
}, {
    "MoleculeName": "C5m2m2",
    "CarrierGas": False,
    "GasPhaseMolFraction": 0.01,
    "isotherms": [["Langmuir", 1.45556, 6.99287e-06]]
}, {
    "MoleculeName": "C5m2m3",
    "CarrierGas": False,
    "GasPhaseMolFraction": 0.01,
    "isotherms": [["Langmuir", 1.53957, 1.41992e-05]]
}, {
    "MoleculeName": "C5m2m4",
    "CarrierGas": False,
    "GasPhaseMolFraction": 0.01,
    "isotherms": [["Langmuir", 1.52841, 2.11233e-05]]
}, {
    "MoleculeName": "C5m3m3",
    "CarrierGas": False,
    "GasPhaseMolFraction": 0.01,
    "isotherms": [["Langmuir", 1.54948, 4.398e-06]]
}])

mix = ruptura.MixturePrediction(components=components,
                                DisplayName="BEA",
                                Temperature=552.0,
                                PressureStart=1e3,
                                PressureEnd=1e7,
                                NumberOfPressurePoints=100,
                                PressureScale="log")

data = mix.compute()
fig, ax = plt.subplots(1,2,figsize=(8,6))
mix.plot(ax[0], "pure")
mix.plot(ax[1], "mixture")
plt.show()