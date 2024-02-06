import ruptura
import matplotlib.pyplot as plt

components = ruptura.Components([{
    "name": "Helium",
    "isCarrierGas": True,
    "gasPhaseMolFraction": 0.93,
    "isotherms": []
}, {
    "name": "nC7",
    "isCarrierGas": False,
    "gasPhaseMolFraction": 0.01,
    "isotherms": [["Langmuir", 1.09984, 6.55857e-5], ["Langmuir", 0.19466, 8.90731e-07]]
}, {
    "name": "C6m2",
    "isCarrierGas": False,
    "gasPhaseMolFraction": 0.01,
    "isotherms": [["Langmuir", 1.22228, 3.90895e-05], ["Langmuir", 0.481726, 9.64046e-08]]
}, {
    "name": "C6m3",
    "isCarrierGas": False,
    "gasPhaseMolFraction": 0.01,
    "isotherms": [["Langmuir", 1.25189, 2.6174e-05], ["Langmuir", 0.518406, 7.99897e-08]]
}, {
    "name": "C5m2m2",
    "isCarrierGas": False,
    "gasPhaseMolFraction": 0.01,
    "isotherms": [["Langmuir", 1.45556, 6.99287e-06]]
}, {
    "name": "C5m2m3",
    "isCarrierGas": False,
    "gasPhaseMolFraction": 0.01,
    "isotherms": [["Langmuir", 1.53957, 1.41992e-05]]
}, {
    "name": "C5m2m4",
    "isCarrierGas": False,
    "gasPhaseMolFraction": 0.01,
    "isotherms": [["Langmuir", 1.52841, 2.11233e-05]]
}, {
    "name": "C5m3m3",
    "isCarrierGas": False,
    "gasPhaseMolFraction": 0.01,
    "isotherms": [["Langmuir", 1.54948, 4.398e-06]]
}])

mix = ruptura.MixturePrediction(components=components,
                                displayName="BEA",
                                temperature=552.0,
                                pressureStart=1e3,
                                pressureEnd=1e7,
                                numberOfPressurePoints=100,
                                pressureScale="log")

data = mix.compute()
fig, ax = plt.subplots(1,2,figsize=(8,6))
mix.plot(ax[0], "pure")
mix.plot(ax[1], "mixture")
plt.show()