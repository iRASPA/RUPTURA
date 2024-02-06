import _ruptura
import numpy as np
from typing import Union, Literal
import matplotlib.pyplot as plt
from .utils import *


class RupturaError(Exception):
    """Base class for other exceptions"""
    pass


class RupturaDataError(RupturaError):
    """Raised when the input data structure is not as expected"""

    def __init__(self, message="Input data structure is not as expected"):
        self.message = message
        super().__init__(self.message)


class RupturaCppError(RupturaError):
    """Raised when there's an error made by the pybind11 C++ code"""

    def __init__(self, message="An error occurred in the Pybind11 C++ code"):
        self.message = message
        super().__init__(self.message)


class Components:
    """
    Class to manage a list of components for simulation.

    Each component is represented as a dictionary with
    various attributes like name, molFraction, isotherms, 
    massTransferCoefficient, axialDispersionCoefficient, 
    and whether it's a carrierGas.

    Attributes:
        components (list): A list of all components (cpp objects) added to this instance.
        labels (list): A list of labels for all components in this instance.
        carrierGas (int): The index of the carrier gas in the components list, if any.
    """

    def __init__(self, components: list[dict] = []):
        """
        Initialize an instance of Components class.

        Parameters:
            components (list[dict], optional): A list of dictionaries representing the components. 
                Defaults to an empty list.
        """

        self.components = []
        self.labels = []
        self.carrierGas = None

        # Add each component in the list
        for comp in components:
            self.addComponent(**comp)

    def addComponent(
        self,
        name: str,
        gasPhaseMolFraction: float,
        isotherms: list = [],
        massTransferCoefficient: float = 0.0,
        axialDispersionCoefficient: float = 0.0,
        isCarrierGas: bool = False,
    ):
        """
        Adds a component to the list of components.

        Parameters:
            name (str): The name of the component.
            gasPhaseMolFraction (float): The mol fraction of the component.
            isotherms (list, optional): List of isotherms. Defaults to empty list.
            massTransferCoefficient (float, optional): The mass transfer coefficient of the component. 
                Defaults to 0.0.
            axialDispersionCoefficient (float, optional): The axial dispersion coefficient of the component.
                Defaults to 0.0.
            isCarrierGas (bool, optional): Specifies whether this component is the carrier gas. 
                Defaults to False.
        """

        # get idx from existing components
        idx = len(self.components)

        # create labels for each isotherm in the component
        for site, isotherm in enumerate(isotherms):
            self.labels += [f"c{idx}_s{site}_{label}" for label in isothermMeta[isotherm[0]]['labels']]

        # add isotherm information
        cpp_isotherms = [_ruptura.Isotherm(isotherm[0], isotherm[1:], len(isotherm) - 1) for isotherm in isotherms]

        # create and append the new component
        comp = _ruptura.Component(
            idx,
            name,
            cpp_isotherms,
            gasPhaseMolFraction,
            massTransferCoefficient,
            axialDispersionCoefficient,
            isCarrierGas,
        )

        # if this component is the carrier gas, update the carrierGas attribute
        if isCarrierGas:
            self.carrierGas = idx

        self.components.append(comp)

    def getLabels(self) -> list[str]:
        """
        Returns a list of labels for all components.

        Returns:
            list[str]: List of labels for all components.
        """
        return [f"{comp.name} (y_i={comp.gasPhaseMolFraction})" for comp in self.components]


class Fitting:
    """
    Class to manage fitting of the component data.

    This class handles the fitting of the component data, 
    which can be loaded from the components or from a file.
    Once the data is fitted, it can be computed and plotted.

    Attributes:
        components (Components): Instance of the Components class.
        fullData (list): Full data for fitting.
        data (np.ndarray): Computed data after fitting.
        Fitting (_ruptura.Fitting): Instance of the ruptura Fitting class.
    """
    ax_labels = [
        "Pressure [Pa]", "Fugacity coeffcient", "Fugacity [Pa]", "Loading (Absolute) [molecules/total cell]",
        "Error [molecules/total cell]", "Loading (Absolute) [molecules/unit cell]", "Error [molecules/unit cell]",
        "Loading (Absolute) [mol/kg framework]", "Error [mol/kg framework]",
        "Loading (Absolute) [milligram/gram framework]", "Error [milligram/gram framework]",
        "Loading (Excess) [molecules/total cell]", "Error [molecules/total cell]",
        "Loading (Excess) [molecules/unit cell]", "Error [molecules/unit cell]", "Loading (Excess) [mol/kg framework]",
        "Error [mol/kg framework]", "Loading (Excess) [milligram/gram framework]", "Error [milligram/gram framework]",
        "heat of desorption [K]", "heat of desorption error [K]"
    ]

    def __init__(self,
                 components: Components,
                 pressureScale: Literal["log", "linear"] = "log",
                 displayName: str = "Column"):
        """
        Initialize an instance of Fitting class.

        Parameters:
            components (Components): An instance of the Components class.
            fullData (list): Full data for fitting.
            pressureScale (str, optional): The scale for pressure: "log", "linear", defaults to "log".
            displayName (str, optional): The name to be displayed, defaults to "Column".
        """
        self.components = components
        self.data = None
        self.displayName = displayName

        # create cpp object
        self.Fitting = _ruptura.Fitting(displayName, components.components, pressureScales[pressureScale])

    def compute(self, data):
        """
        Computes the fitted data and returns it.

        Returns:
            np.ndarray: Computed data.
        """
        self.data = self.Fitting.compute(data)
        return self.data

    def evaluate(self, p):
        evaluatedPoints = self.Fitting.evaluate(p)
        return evaluatedPoints

    def plot(self, ax, data, p, ixlabel, iylabel):
        """
        Plots the fitted data.

        Parameters:
            ax (matplotlib.axes.Axes): The axes object to draw the plot onto.

        Raises:
            AssertionError: If the data has not been computed before plotting.
        """

        loadings = self.Fitting.evaluate(p)
        ax.set_xscale("log")
        ax.set_title(self.displayName)
        ax.set_xlabel(self.ax_labels[ixlabel])
        ax.set_ylabel(self.ax_labels[iylabel])

        labels = [comp.name for comp in self.components.components]
        ncomp = len(labels)

        for i in range(ncomp):
            ax.scatter(*np.array(data[i]).T, label=self.components.components[i].name)
            ax.plot(p, loadings[:, i])
        ax.legend()


class MixturePrediction:
    """
    Class that implements prediction of a mixture of components based on various parameters.

    Attributes:
    -----------
    shape: tuple(int, int, int)
        The shape of the mixture data.
    components: Components
        The components involved in the mixture.
    MixturePrediction: _ruptura.MixturePrediction
        An instance of the MixturePrediction class from the ruptura library.
    data: np.ndarray
        The computed data of the mixture prediction.
    """

    def __init__(
        self,
        components: Components,
        displayName: str = "Column",
        temperature: float = 433.0,
        pressureStart: float = -1.0,
        pressureEnd: float = -1.0,
        numberOfPressurePoints: int = 100,
        pressureScale: Literal["log", "linear"] = "log",
        predictionMethod: Literal["IAST", "SIAST", "EI", "SEI"] = "IAST",
        iastMethod: Literal["FastIAST", "NestedLoopBisection"] = "FastIAST",
    ):
        """
        Initialize the MixturePrediction class with various properties.

        Parameters:
        -----------
        components: Components
            The components involved in the mixture.
        displayName: str
            Name to be displayed for the column (default is "Column").
        temperature: float
            Temperature for the prediction (default is 433.0 K).
        pressureStart: float
            Starting pressure point (default is -1.0).
        pressureEnd: float
            Ending pressure point (default is -1.0).
        numberOfPressurePoints: int
            Number of pressure points to consider (default is 100).
        pressureScale: str
            Scale of pressure: 'log' or 'linear' (default is 'log').
        predictionMethod: str
            Method of prediction: 'IAST', 'SIAST', 'EI', 'SEI' (default is 'IAST').
        iastMethod: str
            IAST Method: 'FastIAST' or 'NestedLoopBisection' (default is 'FastIAST').
        """
        # set shape attribute
        self.shape = (numberOfPressurePoints, len(components.components), 6)

        # select method integers (enum)
        pressureScale = pressureScales[pressureScale]
        predictionMethod = {"IAST": 0, "SIAST": 1, "EI": 2, "SEI": 3}[predictionMethod]
        iastMethod = {"FastIAST": 0, "NestedLoopBisection": 1}[iastMethod]
        self.components = components

        # on carriergas: ruptura first checks if the numberOfCarrierGases is 0 or more. more is
        # always 1, as carrierGasComponent is size_t. if the carriergas is notpresent, set
        # component to 0 (default), it is not checked.

        self.MixturePrediction = _ruptura.MixturePrediction(
            displayName,
            components.components,
            1 if components.carrierGas is not None else 0,
            components.carrierGas or 0,
            temperature,
            pressureStart,
            pressureEnd,
            numberOfPressurePoints,
            pressureScale,
            predictionMethod,
            iastMethod,
        )
        self.data = None

    def compute(self):
        """
        Compute the mixture prediction and store the data.
        
        Returns:
        --------
        np.ndarray
            The computed data of the mixture prediction. Returns array of size (Npress, Ncomp, 6). Last dim has info:
            0: p_i
            1: pure loading
            2: mix loading
            3: mol-frac
            4: adsorbed phase mol-frac 
            5: hypothetical pressure
        """
        self.data = self.MixturePrediction.compute()
        return self.data

    def plot(self, ax, plot_type: Literal["pure", "mixture", "mixture_molfrac"]):
        """
        Plot the mixture prediction data.

        Parameters:
        -----------
        ax: matplotlib.axes.Axes
            The axes on which to draw the plot.
        plot_type: str
            The type of plot: 'pure', 'mixture', 'mixture_molfrac'.
        """

        # get columne
        select = {"pure": 1, "mixture": 2, "mixture_molfrac": 4}[plot_type]

        if self.data is None:
            raise ValueError("Data not computed yet")

        # set axes
        ax.set_xlabel("Total bulk fluid phase fugacity, f/Pa")
        ylabel = {
            "pure": "Absolute loading q_i",
            "mixture": "Absolute loading q_i",
            "mixture_molfrac": "Adsorbed mol-fraction Y_i"
        }
        ax.set_ylabel(ylabel[plot_type])
        ax.set_xscale("log")

        # set values for loop
        ncomp = self.data.shape[1]
        labels = self.components.getLabels()

        # plot all components
        for comp in range(ncomp):
            ax.scatter(self.data[:, comp, 0],
                       self.data[:, comp, select],
                       label=labels[comp],
                       marker=markers[comp],
                       s=8.0)
        ax.legend()


class Breakthrough:
    """
    The Breakthrough class represents a prediction model for gas mixture components. It's designed to predict
    the behaviour of a mixture of components over time as they pass through a column.

    Attributes:
        components (Components): An object containing the components of the gas mixture.
        shape (tuple): The shape of the data that will be predicted by the model.
        Breakthrough (_ruptura.Breakthrough): The Breakthrough model from _ruptura.
        data (np.ndarray): The data predicted by the model, which is computed by calling the `compute` method. Before computation, it's None.
    """

    def __init__(
        self,
        components: Components,
        displayName: str = "Column",
        temperature: float = 433.0,
        numberOfTimeSteps: Union[int, Literal["auto"]] = "auto",
        numberOfGridPoints: int = 100,
        printEvery: int = 10000,
        writeEvery: int = 10000,
        totalPressure: float = 1e6,
        columnVoidFraction: float = 0.4,
        pressureGradient: float = 0.0,
        particleDensity: float = 1e3,
        columnEntranceVelocity: float = 0.1,
        columnLength: float = 0.3,
        timeStep: float = 5e-4,
        pulseTime: float = None,
    ):
        """
        Initializes the Breakthrough object.

        Parameters:
            components (Components): An object containing the components of the gas mixture.
            displayName (str, optional): The display name of the column. Defaults to "Column".
            temperature (float, optional): The temperature of the column. Defaults to 433.0.
            numberOfTimeSteps (Union[int, str], optional): The number of time steps for the prediction model. Can be integer or "auto". Defaults to "auto".
            numberOfGridPoints (int, optional): The number of grid points along the column. Defaults to 100.
            printEvery (int, optional): Frequency of prints. Defaults to 10000.
            writeEvery (int, optional): Frequency of writes. Defaults to 10000.
            totalPressure (float, optional): Total pressure in the column. Defaults to 1e6.
            columnVoidFraction (float, optional): Void fraction of the column. Defaults to 0.4.
            pressureGradient (float, optional): Pressure gradient in the column. Defaults to 0.0.
            particleDensity (float, optional): Density of particles in the column. Defaults to 1e3.
            columnEntranceVelocity (float, optional): Velocity of column entrance. Defaults to 0.1.
            columnLength (float, optional): Length of the column. Defaults to 0.3.
            timeStep (float, optional): Time step for the prediction model. Defaults to 5e-4.
            pulseTime (float, optional): Pulse time for the prediction model. If it's None, there is no pulse. Defaults to None.
        """

        # set attributes
        self.components = components
        self.data = None
        self.shape = (
            numberOfTimeSteps // writeEvery,
            numberOfGridPoints,
            len(components.components) + 1,
            6,
        )

        # determine number of timesteps and set autosteps flag
        autoSteps = numberOfTimeSteps == "auto"
        numberOfTimeSteps = 0 if numberOfTimeSteps == "auto" else int(numberOfTimeSteps)

        # determine pulsetime and pulse flag
        pulse = pulseTime is not None
        pulseTime = 0 if pulseTime is None else pulseTime

        # create mixtureprediction object
        carrierGas = components.carrierGas if components.carrierGas else 0
        mix = MixturePrediction(displayName=displayName, temperature=temperature, components=components)

        # create breakthrough cpp object
        self.Breakthrough = _ruptura.Breakthrough(
            displayName,
            components.components,
            carrierGas,
            numberOfGridPoints,
            printEvery,
            writeEvery,
            temperature,
            totalPressure,
            columnVoidFraction,
            pressureGradient,
            particleDensity,
            columnEntranceVelocity,
            columnLength,
            timeStep,
            numberOfTimeSteps,
            autoSteps,
            pulse,
            pulseTime,
            mix.MixturePrediction,
        )

    def compute(self):
        """
        Computes the data for the breakthrough model. The computed data is stored in the class attribute 'data'.
        
        Returns:
            np.ndarray: The computed data.
        """
        self.data = self.Breakthrough.compute()
        return self.data

    def plot(self, ax, plot_type: Literal["breakthrough", "Dpdt", "Dqdt", "P", "Pnorm", "Pt", "Q", "Qeq", "V"]):
        """
        Plots the data for the breakthrough model. If data has not yet been computed, raises a ValueError.

        Parameters:
            ax (matplotlib.axes.Axes): The axes on which to plot the data.
            plot_type (str): The type of plot to produce. Can be one of "breakthrough", "Dpdt", "Dqdt", "P",
                             "Pnorm", "Pt", "Q", "Qeq", "V".
        
        Raises:
            ValueError: If the data has not yet been computed.
        """

        if self.data is None:
            raise ValueError("Data not computed yet")

        # set values for loop
        labels = self.components.getLabels()
        ncomp = len(labels)

        if plot_type == "breakthrough":

            x = self.data[:, -1, 0]
            ax.set_xlim(x.min(), x.max())
            ax.set_xticks(np.linspace(x.min(), x.max(), 7))
            ax.set_xlabel("Dimensionless time")
            ax.set_ylabel("Concentration exit gas")

            x2 = self.data[:, -1, 1]
            ax2 = ax.twiny()
            ax2.set_xticks(np.linspace(x2.min(), x2.max(), 7))
            ax2.set_xlabel("Time (min)")

            for comp in range(ncomp):
                ax.scatter(x, self.data[:, -1, 8 + comp * 6], label=labels[comp], marker=getMarker(comp), s=8.0)

        ax.legend()
