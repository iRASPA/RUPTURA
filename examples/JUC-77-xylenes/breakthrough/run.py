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
    "isotherms": [["Langmuir-Freundlich", 0.4, 2.4e-3, 0.8], ["Langmuir", 0.9, 2.01e-2]]
}, {
    "MoleculeName": "m-xylene",
    "GasPhaseMolFraction": 0.01,
    "MassTransferCoefficient": 0.06,
    "AxialDispersionCoefficient": 0.0,
    "isotherms": [["Langmuir-Freundlich", 0.53, 1.24e-2, 0.822], ["Langmuir-Freundlich", 0.77, 3.1e-2, 1.1]]
}, {
    "MoleculeName": "p-xylene",
    "GasPhaseMolFraction": 0.01,
    "MassTransferCoefficient": 0.06,
    "AxialDispersionCoefficient": 0.0,
    "isotherms": [["Langmuir-Freundlich", 0.6, 1.24e-2, 0.9], ["Langmuir-Freundlich", 0.7, 3.14e-2, 1.05]]
}, {
    "MoleculeName": "ethylbenzene",
    "GasPhaseMolFraction": 0.01,
    "MassTransferCoefficient": 0.06,
    "AxialDispersionCoefficient": 0.0,
    "isotherms": [["Langmuir-Freundlich", 0.4, 7.02e-4, 0.9], ["Langmuir-Freundlich", 0.9, 2.98e-3, 1.2]]
}])

brk = ruptura.Breakthrough(
    components=components,
    DisplayName="JUC-77",
    Temperature=443.0,
    ParticleDensity=1144.03,
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