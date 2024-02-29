import ruptura
import matplotlib.pyplot as plt

components = ruptura.Components([{
    "MoleculeName": "Helium",
    "CarrierGas": True,
    "GasPhaseMolFraction": 0.93,
    "isotherms": []
}, {
    "MoleculeName":
        "1-3-butadiene",
    "GasPhaseMolFraction": 0.01,
    "MassTransferCoefficient": 0.06,
    "AxialDispersionCoefficient": 0.0,
    "isotherms": [["Langmuir-Freundlich", 1.77359, 0.964548, 0.891275],
                  ["Langmuir-Freundlich", 4.71979, 4.193e+06, 1.76902],
                  ["Langmuir-Freundlich", 0.874237, 2857.35, 1.65435]]
}, {
    "MoleculeName":
        "1-butene",
    "GasPhaseMolFraction": 0.01,
    "MassTransferCoefficient": 0.06,
    "AxialDispersionCoefficient": 0.0,
    "isotherms": [["Langmuir-Freundlich", 4.12029, 260699.0, 2.19499],
                  ["Langmuir-Freundlich", 2.33644, 0.953314, 0.431999]]
}, {
    "MoleculeName":
        "2-cis-butene",
    "GasPhaseMolFraction": 0.01,
    "MassTransferCoefficient": 0.06,
    "AxialDispersionCoefficient": 0.0,
    "isotherms": [["Langmuir-Freundlich", 5.25517, 3.14573e+06, 2.10887],
                  ["Langmuir-Freundlich", 1.28055, 7.58487, 0.75982]]
}, {
    "MoleculeName":
        "2-trans-butene",
    "GasPhaseMolFraction": 0.01,
    "MassTransferCoefficient": 0.06,
    "AxialDispersionCoefficient": 0.0,
    "isotherms": [["Langmuir-Freundlich", 3.81848, 2.12583e+08, 3.11845],
                  ["Langmuir-Freundlich", 2.95468, 3.58644, 0.531051]]
}])

brk = ruptura.Breakthrough(
    components=components,
    DisplayName="Fe-MOF-74",
    Temperature=298.0,
    ParticleDensity=1126.7,
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