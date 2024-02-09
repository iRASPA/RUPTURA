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
    name = fn.rstrip()
    components.addComponent(
        name = name,
        GasPhaseMolFraction=1.0,
        isotherms = [["Bingel-Walton", 0.0, 0.0, 0.0]]
    )
    print(np.genfromtxt(fn).shape)
    loadings.append([(x,y) for x, y in np.genfromtxt(fn)])

fitting = ruptura.Fitting(
    components=components,
    DisplayName="CoBDP",
)

fitting.compute(loadings)
p = np.linspace(0, 7, 200)

fig, ax = plt.subplots(figsize=(12,6))
fitting.plot(ax, loadings, p, *select)
ax.set_xscale("linear")
plt.show()