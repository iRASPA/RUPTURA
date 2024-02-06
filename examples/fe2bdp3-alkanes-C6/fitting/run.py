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
    with open(fn, 'r') as f:
        name = f.readline().split(" ")[3].rstrip("\n")
    components.addComponent(
        name = name,
        gasPhaseMolFraction=1.0,
        isotherms = [["Langmuir-Freundlich", 0.0, 0.0, 0.0], ["Langmuir-Freundlich", 0.0, 0.0, 0.0]]
    )
    loadings.append([(x,y) for x, y in np.genfromtxt(fn)[:, select]])

fitting = ruptura.Fitting(
    components=components,
    displayName="fe2bdp3",
)

fitting.compute(loadings)
p = np.logspace(-2, 10, 200)

fig, ax = plt.subplots(figsize=(12,6))
fitting.plot(ax, loadings, p, *select)
plt.show()