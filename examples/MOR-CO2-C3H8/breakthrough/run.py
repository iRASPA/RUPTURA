import ruptura
import matplotlib.pyplot as plt

ruptura_objs = ruptura.from_input("simulation.input")
ruptura_objs["Breakthrough"].compute()

fig, ax = plt.subplots()
ruptura_objs["Breakthrough"].plot(ax, "breakthrough")
plt.show()