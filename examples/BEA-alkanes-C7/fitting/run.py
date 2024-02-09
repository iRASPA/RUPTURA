import ruptura
import numpy as np
import matplotlib.pyplot as plt
import sys

components = ruptura.Components()
loadings = []

# select quantities to fit to
select = [2, 3]

if len(sys.argv) < 2:
    print("Usage: python run.py filename1 filename2 ...")
    sys.exit(0)
for fn in sys.argv[1:]:
    with open(fn, 'r') as f:
        name = f.readline().split(" ")[1]
    components.addComponent(name=name, GasPhaseMolFraction=1.0, isotherms=[["Toth", 0.0, 0.0, 0.0]])
    loadings.append([(x, y) for x, y in np.genfromtxt(fn)[:, select]])

fitting = ruptura.Fitting(
    components=components,
    DisplayName="BEA",
)

fitting.compute(loadings)
p = np.logspace(0, 8, 100)

fig, ax = plt.subplots(figsize=(12, 6))
fitting.plot(ax, loadings, p, *select)
plt.show()
