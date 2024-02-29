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
    "isotherms": [["Langmuir-Freundlich", 1.7, 7.74e-1, 0.4], ["Langmuir-Freundlich", 1.8, 5.53e6, 1.1]]
}, {
    "MoleculeName": "m-xylene",
    "GasPhaseMolFraction": 0.01,
    "MassTransferCoefficient": 0.06,
    "AxialDispersionCoefficient": 0.0,
    "isotherms": [["Langmuir-Freundlich", 1.4, 1.22e1, 0.57], ["Langmuir-Freundlich", 1.8, 2.68e6, 1.07]]
}, {
    "MoleculeName": "p-xylene",
    "GasPhaseMolFraction": 0.01,
    "MassTransferCoefficient": 0.06,
    "AxialDispersionCoefficient": 0.0,
    "isotherms": [["Langmuir-Freundlich", 1.6, 1.07e1, 0.54], ["Langmuir-Freundlich", 1.8, 4.97e7, 1.2]]
}, {
    "MoleculeName": "ethylbenzene",
    "GasPhaseMolFraction": 0.01,
    "MassTransferCoefficient": 0.06,
    "AxialDispersionCoefficient": 0.0,
    "isotherms": [["Langmuir-Freundlich", 1.5, 7.53, 0.73], ["Langmuir", 1.3, 3.76e4]]
}])

brk = ruptura.Breakthrough(
    components=components,
    DisplayName="MAF-X8",
    Temperature=443.0,
    ParticleDensity=954.29,
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