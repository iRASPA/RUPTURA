import ruptura
import matplotlib.pyplot as plt

components = ruptura.Components([{
    "MoleculeName": "Helium",
    "CarrierGas": True,
    "GasPhaseMolFraction": 0.93,
    "isotherms": []
}, {
    "MoleculeName": "o-xylene",
    "GasPhaseMolFraction": 0.01,
    "MassTransferCoefficient": 0.06,
    "AxialDispersionCoefficient": 0.0,
    "isotherms": [["Langmuir-Freundlich", 4.8, 1.28e-5, 1.4], ["Langmuir-Freundlich", 1.5, 1.6e-4, 0.7]]
}, {
    "MoleculeName": "m-xylene",
    "GasPhaseMolFraction": 0.01,
    "MassTransferCoefficient": 0.06,
    "AxialDispersionCoefficient": 0.0,
    "isotherms": [["Langmuir-Freundlich", 4.8, 3.19e-5, 1.27], ["Langmuir-Freundlich", 1.2, 4.72e-5, 0.7]]
}, {
    "MoleculeName": "p-xylene",
    "GasPhaseMolFraction": 0.01,
    "MassTransferCoefficient": 0.06,
    "AxialDispersionCoefficient": 0.0,
    "isotherms": [["Langmuir-Freundlich", 4.5, 2.3e-6, 1.7], ["Langmuir-Freundlich", 1.6, 1.46e-4, 0.7]]
}, {
    "MoleculeName": "ethylbenzene",
    "GasPhaseMolFraction": 0.01,
    "MassTransferCoefficient": 0.06,
    "AxialDispersionCoefficient": 0.0,
    "isotherms": [["Langmuir-Freundlich", 4.5, 7.38e-6, 1.5], ["Langmuir-Freundlich", 1.4, 1.6e-4, 0.7]]
}])

brk = ruptura.Breakthrough(
    components=components,
    DisplayName="CoBDP",
    Temperature=443.0,
    ParticleDensity=721.88,
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