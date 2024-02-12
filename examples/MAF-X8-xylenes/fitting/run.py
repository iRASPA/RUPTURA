import ruptura
import numpy as np
import matplotlib.pyplot as plt
import sys

components = ruptura.Components()

# select quantities to fit to
select = [2, 3]

loadings = []
if len(sys.argv) < 2:
    print("Usage: python run.py filename1 filename2 ...")
    sys.exit(0)
for fn in sys.argv[1:]:
    name = fn.split("-")[-1]
    components.addComponent(
        MoleculeName = name,
        GasPhaseMolFraction=1.0,
        isotherms = [["Langmuir-Freundlich", 0.0, 0.0, 0.0], ["Langmuir-Freundlich", 0.0, 0.0, 0.0], ["Langmuir-Freundlich", 0.0, 0.0, 0.0]]
    )
    loadings.append([(x,y) for x, y in np.genfromtxt(fn, usecols=list(range(21)), invalid_raise=False)[:, select]])

fitting = ruptura.Fitting(
    components=components,
    DisplayName="MAF-X8",
)

fitting.compute(loadings)
p = np.logspace(-8, 10, 200)

fig, ax = plt.subplots(figsize=(12,6))
fitting.plot(ax, loadings, p, *select)
plt.show()