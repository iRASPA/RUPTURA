import ruptura
import matplotlib.pyplot as plt

components = ruptura.Components([{
    "MoleculeName": "Helium",
    "CarrierGas": True,
    "GasPhaseMolFraction": 0.96,
    "isotherms": []
}, {
    "MoleculeName": "o-xylene",
    "GasPhaseMolFraction": 0.01,
    "MassTransferCoefficient": 0.06,
    "AxialDispersionCoefficient": 0.0,
    "isotherms": [["Langmuir-Freundlich", 2.8, 5.03e-3, 1.08], ["Langmuir-Freundlich", 0.8, 3.78e-4, 0.8]]
}, {
    "MoleculeName": "m-xylene",
    "GasPhaseMolFraction": 0.01,
    "MassTransferCoefficient": 0.06,
    "AxialDispersionCoefficient": 0.0,
    "isotherms": [["Langmuir-Freundlich", 2.6, 2.48e-3, 1.1], ["Langmuir-Freundlich", 0.8, 3.86e-5, 0.9]]
}, {
    "MoleculeName": "p-xylene",
    "GasPhaseMolFraction": 0.01,
    "MassTransferCoefficient": 0.06,
    "AxialDispersionCoefficient": 0.0,
    "isotherms": [["Langmuir-Freundlich", 2.8, 3.22e-3, 1.07], ["Langmuir-Freundlich", 0.8, 3.23e-5, 0.9]]
}, {
    "MoleculeName": "ethylbenzene",
    "GasPhaseMolFraction": 0.01,
    "MassTransferCoefficient": 0.06,
    "AxialDispersionCoefficient": 0.0,
    "isotherms": [["Langmuir-Freundlich", 0.9, 3.94e-4, 0.62], ["Langmuir-Freundlich", 2.2, 4.73e-3, 1.05]]
}])

brk = ruptura.Breakthrough(
    components=components,
    DisplayName="MIL-47",
    Temperature=443.0,
    ParticleDensity=1000.36,
    TotalPressure=2.5e6,
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