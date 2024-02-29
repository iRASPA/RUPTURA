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
    "isotherms": [["Langmuir-Freundlich", 1.8, 1.59e-5, 0.62], ["Langmuir-Freundlich", 2.9, 1.5e-4, 1.06]]
}, {
    "MoleculeName": "m-xylene",
    "GasPhaseMolFraction": 0.01,
    "MassTransferCoefficient": 0.06,
    "AxialDispersionCoefficient": 0.0,
    "isotherms": [["Langmuir-Freundlich", 1.1, 1.69e-5, 0.68], ["Langmuir-Freundlich", 3.0, 9.98e-5, 1.1]]
}, {
    "MoleculeName": "p-xylene",
    "GasPhaseMolFraction": 0.01,
    "MassTransferCoefficient": 0.06,
    "AxialDispersionCoefficient": 0.0,
    "isotherms": [["Langmuir-Freundlich", 1.2, 3.61e-6, 0.87], ["Langmuir-Freundlich", 3.3, 8.33e-5, 1.15]]
}, {
    "MoleculeName": "ethylbenzene",
    "GasPhaseMolFraction": 0.01,
    "MassTransferCoefficient": 0.06,
    "AxialDispersionCoefficient": 0.0,
    "isotherms": [["Langmuir-Freundlich", 1.3, 2.61e-6, 0.9], ["Langmuir-Freundlich", 3.1, 8.02e-5, 1.1]]
}])

brk = ruptura.Breakthrough(
    components=components,
    DisplayName="MIL-125(Ti)-NH_2",
    Temperature=443.0,
    ParticleDensity=861.6161,
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
