import ruptura
import matplotlib.pyplot as plt

components = ruptura.Components([{
    "name":
        "nC4",
    "gasPhaseMolFraction":
        0.06666666666667,
    "massTransferCoefficient":
        0.06,
    "axialDispersionCoefficient":
        0.0,
    "isotherms": [["Langmuir-Freundlich", 0.777988, 0.0028668, 0.401129],
                  ["Langmuir-Freundlich", 1.14543, 0.000254606, 0.75763]]
}, {
    "name":
        "2mC3",
    "gasPhaseMolFraction":
        0.06666666666667,
    "massTransferCoefficient":
        0.06,
    "axialDispersionCoefficient":
        0.0,
    "isotherms": [["Langmuir-Freundlich", 1.70797, 9.79017e-05, 0.643565],
                  ["Langmuir-Freundlich", 0.458558, 5.62846e-06, 1.14944]]
}, {
    "name":
        "nC5",
    "gasPhaseMolFraction":
        0.06666666666667,
    "massTransferCoefficient":
        0.06,
    "axialDispersionCoefficient":
        0.0,
    "isotherms": [["Langmuir-Freundlich", 0.937968, 0.000530062, 0.841165],
                  ["Langmuir-Freundlich", 0.634205, 4.36785e-05, 0.742879]]
}, {
    "name":
        "2mC4",
    "gasPhaseMolFraction":
        0.06666666666667,
    "massTransferCoefficient":
        0.06,
    "axialDispersionCoefficient":
        0.0,
    "isotherms": [["Langmuir-Freundlich", 0.714103, 2.69729e-05, 0.671999],
                  ["Langmuir-Freundlich", 0.771593, 4.05931e-05, 0.927118]]
}, {
    "name": "22dmC3",
    "gasPhaseMolFraction": 0.06666666666667,
    "massTransferCoefficient": 0.06,
    "axialDispersionCoefficient": 0.0,
    "isotherms": [["Langmuir-Freundlich", 0.690561, 5.13286e-08, 1.0147]]
}, {
    "name":
        "nC6",
    "gasPhaseMolFraction":
        0.06666666666667,
    "massTransferCoefficient":
        0.06,
    "axialDispersionCoefficient":
        0.0,
    "isotherms": [["Langmuir-Freundlich", 0.869905, 0.00316758, 0.475264],
                  ["Langmuir-Freundlich", 0.578986, 0.000238141, 1.11479]]
}, {
    "name":
        "2MP",
    "gasPhaseMolFraction":
        0.06666666666667,
    "massTransferCoefficient":
        0.06,
    "axialDispersionCoefficient":
        0.0,
    "isotherms": [["Langmuir-Freundlich", 0.960776, 0.000296532, 0.788678],
                  ["Langmuir-Freundlich", 0.590593, 0.000100636, 0.460252]]
}, {
    "name":
        "3MP",
    "gasPhaseMolFraction":
        0.06666666666667,
    "massTransferCoefficient":
        0.06,
    "axialDispersionCoefficient":
        0.0,
    "isotherms": [["Langmuir-Freundlich", 0.532722, 9.61523e-05, 0.405003],
                  ["Langmuir-Freundlich", 1.02932, 0.000949822, 0.58711]]
}, {
    "name":
        "23DMB",
    "gasPhaseMolFraction":
        0.06666666666667,
    "massTransferCoefficient":
        0.06,
    "axialDispersionCoefficient":
        0.0,
    "isotherms": [["Langmuir-Freundlich", 0.692959, 5.59469e-09, 1.3469],
                  ["Langmuir-Freundlich", 0.53868, 0.00056127, 0.423081]]
}, {
    "name":
        "22DMB",
    "gasPhaseMolFraction":
        0.06666666666667,
    "massTransferCoefficient":
        0.06,
    "axialDispersionCoefficient":
        0.0,
    "isotherms": [["Langmuir-Freundlich", 0.644694, 6.16965e-10, 1.21463],
                  ["Langmuir-Freundlich", 0.884571, 5.02998e-09, 0.727642]]
}, {
    "name":
        "nC7",
    "gasPhaseMolFraction":
        0.06666666666667,
    "massTransferCoefficient":
        0.06,
    "axialDispersionCoefficient":
        0.0,
    "isotherms": [["Langmuir-Freundlich", 0.688052, 0.00525971, 0.402162],
                  ["Langmuir-Freundlich", 0.622752, 0.00166728, 0.990279]]
}, {
    "name":
        "2MH",
    "gasPhaseMolFraction":
        0.06666666666667,
    "massTransferCoefficient":
        0.06,
    "axialDispersionCoefficient":
        0.0,
    "isotherms": [["Langmuir-Freundlich", 0.310174, 0.00344507, 0.27902],
                  ["Langmuir-Freundlich", 0.932413, 0.00042831, 0.794345]]
}, {
    "name":
        "3MH",
    "gasPhaseMolFraction":
        0.06666666666667,
    "massTransferCoefficient":
        0.06,
    "axialDispersionCoefficient":
        0.0,
    "isotherms": [["Langmuir-Freundlich", 0.627761, 2.0218e-05, 1.10607],
                  ["Langmuir-Freundlich", 0.492039, 0.00563881, 0.341391]]
}, {
    "name": "23DMP",
    "gasPhaseMolFraction": 0.06666666666667,
    "massTransferCoefficient": 0.06,
    "axialDispersionCoefficient": 0.0,
    "isotherms": [["Langmuir-Freundlich", 0.983143, 1.98894e-05, 0.78231]]
}, {
    "name":
        "22DMP",
    "gasPhaseMolFraction":
        0.06666666666667,
    "massTransferCoefficient":
        0.06,
    "axialDispersionCoefficient":
        0.0,
    "isotherms": [["Langmuir-Freundlich", 0.673428, 4.69463e-10, 1.27279],
                  ["Langmuir-Freundlich", 0.428088, 3.8364e-10, 0.888983]]
}])

mix = ruptura.MixturePrediction(components=components,
                                displayName="ZIF-77",
                                temperature=443.0,
                                pressureStart=1e0,
                                pressureEnd=1e6,
                                numberOfPressurePoints=100,
                                pressureScale="log")

data = mix.compute()
fig, ax = plt.subplots(1, 2, figsize=(8, 6))
mix.plot(ax[0], "pure")
mix.plot(ax[1], "mixture")
plt.show()
