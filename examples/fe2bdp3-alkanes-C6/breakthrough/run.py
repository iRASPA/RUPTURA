import ruptura
import matplotlib.pyplot as plt

components = ruptura.Components([{
    "MoleculeName": "Helium",
    "CarrierGas": True,
    "GasPhaseMolFraction": 0.95,
    "isotherms": []
}, {
    "MoleculeName": "nC6",
    "GasPhaseMolFraction": 0.01,
    "MassTransferCoefficient": 0.06,
    "AxialDispersionCoefficient": 0.0,
    "isotherms": [["Langmuir-Freundlich", 1.1, 2.14e-4, 0.6], ["Langmuir-Freundlich", 4.8, 2.39e-5, 1.35]]
}, {
    "MoleculeName": "2MP",
    "GasPhaseMolFraction": 0.01,
    "MassTransferCoefficient": 0.06,
    "AxialDispersionCoefficient": 0.0,
    "isotherms": [["Langmuir-Freundlich", 1.4, 1.68e-4, 0.6], ["Langmuir-Freundlich", 4.8, 3.3e-5, 1.26]]
}, {
    "MoleculeName": "3MP",
    "GasPhaseMolFraction": 0.01,
    "MassTransferCoefficient": 0.06,
    "AxialDispersionCoefficient": 0.0,
    "isotherms": [["Langmuir-Freundlich", 1.75, 1.14e-4, 0.62], ["Langmuir-Freundlich", 4.8, 3.22e-5, 1.24]]
}, {
    "MoleculeName": "23DMB",
    "GasPhaseMolFraction": 0.01,
    "MassTransferCoefficient": 0.06,
    "AxialDispersionCoefficient": 0.0,
    "isotherms": [["Langmuir-Freundlich", 1.5, 8.31e-5, 0.63], ["Langmuir-Freundlich", 4.8, 5.46e-5, 1.13]]
}, {
    "MoleculeName": "22DMB",
    "GasPhaseMolFraction": 0.01,
    "MassTransferCoefficient": 0.06,
    "AxialDispersionCoefficient": 0.0,
    "isotherms": [["Langmuir-Freundlich", 1.74, 5.65e-6, 0.72], ["Langmuir-Freundlich", 4.8, 4.3e-5, 1.04]]
}])

brk = ruptura.Breakthrough(
    components=components,
    DisplayName="fe2bdp3",
    Temperature=443.0,
    ParticleDensity=1146.041505,
    TotalPressure=2.0e6,
    PressureGradient=0.0,
    ColumnEntranceVelocity=0.1,
    ColumnLength=0.3,
    NumberOfTimeSteps="auto",
    TimeStep=5e-4,
)

data = brk.compute()
fig, ax = plt.subplots(figsize=(8,6))
brk.plot(ax, "breakthrough")
plt.show()