import _ruptura
import numpy as np
from typing import Union, Literal
import matplotlib.pyplot as plt
from .utils import *

def from_input(file_path):
    return from_config(input_to_config(file_path))

def from_config(config):
    ruptura_objects = {"Components": Components(config["components"])}
    if "Fitting" in config.keys():
        ruptura_objects["Fitting"] = Fitting(components=ruptura_objects["Components"], **config["Fitting"])
    if "Breakthrough" in config.keys():
        ruptura_objects["Breakthrough"] = Breakthrough(components=ruptura_objects["Components"], **config["Breakthrough"])
    if "MixturePrediction" in config.keys():
        ruptura_objects["MixturePrediction"] = MixturePrediction(components=ruptura_objects["Components"], **config["MixturePrediction"])
    return ruptura_objects


class Components:
    """
    Class to manage a list of components for simulation.

    Each component is represented as a dictionary with
    various attributes like name, molFraction, isotherms, 
    MassTransferCoefficient, AxialDispersionCoefficient, 
    and whether it's a CarrierGas.

    Attributes:
        components (list): A list of all components (cpp objects) added to this instance.
        labels (list): A list of labels for all components in this instance.
        CarrierGas (int): The index of the carrier gas in the components list, if any.
    """

    def __init__(self, components: list[dict] = []):
        """
        Initialize an instance of Components class.

        Parameters:
            components (list[dict], optional): A list of dictionaries representing the components. 
                Defaults to an empty list.
        """

        self._components = []
        self.labels = []
        self.CarrierGas = None

        # Add each component in the list
        for comp in components:
            self.addComponent(**comp)

    def addComponent(
        self,
        MoleculeName: str,
        GasPhaseMolFraction: float,
        isotherms: list = [],
        MassTransferCoefficient: float = 0.0,
        AxialDispersionCoefficient: float = 0.0,
        CarrierGas: bool = False,
    ):
        """
        Adds a component to the list of components.

        Parameters:
            MoleculeName (str): The name of the component.
            GasPhaseMolFraction (float): The mol fraction of the component.
            isotherms (list, optional): List of isotherms. Defaults to empty list.
            MassTransferCoefficient (float, optional): The mass transfer coefficient of the component. 
                Defaults to 0.0.
            AxialDispersionCoefficient (float, optional): The axial dispersion coefficient of the component.
                Defaults to 0.0.
            CarrierGas (bool, optional): Specifies whether this component is the carrier gas. 
                Defaults to False.
        """

        # check valid molfrac
        if GasPhaseMolFraction < 0.0:
            raise ValueError("Mol fraction can not be negative!")

        # get idx from existing components
        idx = len(self._components)

        # create labels for each isotherm in the component
        for site, isotherm in enumerate(isotherms):
            self.labels += [f"c{idx}_s{site}_{label}" for label in isothermMeta[isotherm[0]]['labels']]

        # add isotherm information
        _cpp_isotherms = [_ruptura.Isotherm(isotherm[0], isotherm[1:], len(isotherm) - 1) for isotherm in isotherms]

        # create and append the new component
        _comp = _ruptura.Component(
            idx,
            MoleculeName,
            _cpp_isotherms,
            GasPhaseMolFraction,
            MassTransferCoefficient,
            AxialDispersionCoefficient,
            CarrierGas,
        )

        # if this component is the carrier gas, update the CarrierGas attribute
        if CarrierGas:
            self.CarrierGas = idx

        self._components.append(_comp)

    def getLabels(self) -> list[str]:
        """
        Returns a list of labels for all components.

        Returns:
            list[str]: List of labels for all components.
        """
        return [f"{_comp.MoleculeName} (y_i={_comp.GasPhaseMolFraction:3.2e})" for _comp in self._components]


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
                 PressureScale: Literal["log", "linear"] = "log",
                 DisplayName: str = "Column"):
        """
        Initialize an instance of Fitting class.

        Parameters:
            components (Components): An instance of the Components class.
            fullData (list): Full data for fitting.
            PressureScale (str, optional): The scale for pressure: "log", "linear", defaults to "log".
            DisplayName (str, optional): The name to be displayed, defaults to "Column".
        """
        self.components = components
        self.data = None
        self.DisplayName = DisplayName

        # create cpp object
        self._Fitting = _ruptura.Fitting(DisplayName, components._components, pressureScales[PressureScale])

    def compute(self, data: list[list[tuple]]) -> np.ndarray:
        """
        Computes the fitted data and returns it.

        Returns:
            np.ndarray: Computed data.
        """
        self.data = self._Fitting.compute(data)
        return self.data

    def evaluate(self, p: np.ndarray) -> np.ndarray:
        evaluatedPoints = self._Fitting.evaluate(p)
        return evaluatedPoints

    def plot(self, ax, data, p, ixlabel, iylabel):
        """
        Plots the fitted data.

        Parameters:
            ax (matplotlib.axes.Axes): The axes object to draw the plot onto.

        Raises:
            AssertionError: If the data has not been computed before plotting.
        """

        loadings = self._Fitting.evaluate(p)
        ax.set_xscale("log")
        ax.set_title(self.DisplayName)
        ax.set_xlabel(self.ax_labels[ixlabel])
        ax.set_ylabel(self.ax_labels[iylabel])

        labels = self.components.getLabels()
        ncomp = len(labels)

        for i in range(ncomp):
            ax.scatter(*np.array(data[i]).T, label=labels[i])
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
        DisplayName: str = "Column",
        Temperature: float = 433.0,
        PressureStart: float = -1.0,
        PressureEnd: float = -1.0,
        NumberOfPressurePoints: int = 100,
        PressureScale: Literal["log", "linear"] = "log",
        MixturePredictionMethod: Literal["IAST", "SIAST", "EI", "SEI"] = "IAST",
        IASTMethod: Literal["FastIAST", "NestedLoopBisection"] = "FastIAST",
    ):
        """
        Initialize the MixturePrediction class with various properties.

        Parameters:
        -----------
        components: Components
            The components involved in the mixture.
        DisplayName: str
            Name to be displayed for the column (default is "Column").
        Temperature: float
            Temperature for the prediction (default is 433.0 K).
        PressureStart: float
            Starting pressure point (default is -1.0).
        PressureEnd: float
            Ending pressure point (default is -1.0).
        NumberOfPressurePoints: int
            Number of pressure points to consider (default is 100).
        PressureScale: str
            Scale of pressure: 'log' or 'linear' (default is 'log').
        MixturePredictionMethod: str
            Method of prediction: 'IAST', 'SIAST', 'EI', 'SEI' (default is 'IAST').
        IASTMethod: str
            IAST Method: 'FastIAST' or 'NestedLoopBisection' (default is 'FastIAST').
        """
        # set shape attribute
        self.shape = (NumberOfPressurePoints, len(components._components), 6)
        self.DisplayName=DisplayName

        # select method integers (enum)
        PressureScale = pressureScales[PressureScale]
        MixturePredictionMethod = {"IAST": 0, "SIAST": 1, "EI": 2, "SEI": 3}[MixturePredictionMethod]
        IASTMethod = {"FastIAST": 0, "NestedLoopBisection": 1}[IASTMethod]
        self.components = components

        # on carriergas: ruptura first checks if the numberOfCarrierGases is 0 or more. more is
        # always 1, as CarrierGasComponent is size_t. if the carriergas is notpresent, set
        # component to 0 (default), it is not checked.

        self._MixturePrediction = _ruptura.MixturePrediction(
            DisplayName,
            components._components,
            1 if components.CarrierGas is not None else 0,
            components.CarrierGas or 0,
            Temperature,
            PressureStart,
            PressureEnd,
            NumberOfPressurePoints,
            PressureScale,
            MixturePredictionMethod,
            IASTMethod,
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
        self.data = self._MixturePrediction.compute()
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
        ax.set_title(self.DisplayName + f" ({plot_type})")
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
                       marker=getMarker(comp),
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
        DisplayName: str = "Column",
        Temperature: float = 433.0,
        NumberOfTimeSteps: Union[int, Literal["auto"]] = "auto",
        NumberOfGridPoints: int = 100,
        PrintEvery: int = 10000,
        WriteEvery: int = 10000,
        TotalPressure: float = 1e6,
        ColumnVoidFraction: float = 0.4,
        PressureGradient: float = 0.0,
        ParticleDensity: float = 1e3,
        ColumnEntranceVelocity: float = 0.1,
        ColumnLength: float = 0.3,
        TimeStep: float = 5e-4,
        PulseTime: float = None,
        MixturePredictionMethod: Literal["IAST", "SIAST", "EI", "SEI"] = "IAST",
    ):
        """
        Initializes the Breakthrough object.

        Parameters:
            components (Components): An object containing the components of the gas mixture.
            DisplayName (str, optional): The display name of the column. Defaults to "Column".
            Temperature (float, optional): The Temperature of the column. Defaults to 433.0.
            NumberOfTimeSteps (Union[int, str], optional): The number of time steps for the prediction model. Can be integer or "auto". Defaults to "auto".
            NumberOfGridPoints (int, optional): The number of grid points along the column. Defaults to 100.
            PrintEvery (int, optional): Frequency of prints. Defaults to 10000.
            WriteEvery (int, optional): Frequency of writes. Defaults to 10000.
            TotalPressure (float, optional): Total pressure in the column. Defaults to 1e6.
            ColumnVoidFraction (float, optional): Void fraction of the column. Defaults to 0.4.
            PressureGradient (float, optional): Pressure gradient in the column. Defaults to 0.0.
            ParticleDensity (float, optional): Density of particles in the column. Defaults to 1e3.
            ColumnEntranceVelocity (float, optional): Velocity of column entrance. Defaults to 0.1.
            ColumnLength (float, optional): Length of the column. Defaults to 0.3.
            TimeStep (float, optional): Time step for the prediction model. Defaults to 5e-4.
            PulseTime (float, optional): Pulse time for the prediction model. If it's None, there is no pulse. Defaults to None.
            MixturePredictionMethod (str): Method of prediction: 'IAST', 'SIAST', 'EI', 'SEI' (default is 'IAST').
        """

        # set attributes
        self.components = components
        self.DisplayName = DisplayName
        self.data = None
        
        # determine number of timesteps and set autosteps flag
        autoSteps = NumberOfTimeSteps == "auto"
        NumberOfTimeSteps = 0 if NumberOfTimeSteps == "auto" else int(NumberOfTimeSteps)

        # determine pulsetime and pulse flag
        pulse = PulseTime is not None
        PulseTime = 0 if PulseTime is None else PulseTime

        # create mixtureprediction object
        CarrierGas = components.CarrierGas if components.CarrierGas else 0
        mix = MixturePrediction(DisplayName=DisplayName, Temperature=Temperature, components=components, MixturePredictionMethod=MixturePredictionMethod)

        # create breakthrough cpp object
        self._Breakthrough = _ruptura.Breakthrough(
            DisplayName,
            components._components,
            CarrierGas,
            NumberOfGridPoints,
            PrintEvery,
            WriteEvery,
            Temperature,
            TotalPressure,
            ColumnVoidFraction,
            PressureGradient,
            ParticleDensity,
            ColumnEntranceVelocity,
            ColumnLength,
            TimeStep,
            NumberOfTimeSteps,
            autoSteps,
            pulse,
            PulseTime,
            mix._MixturePrediction,
        )

    def compute(self):
        """
        Computes the data for the breakthrough model. The computed data is stored in the class attribute 'data'.
        
        Returns:
            np.ndarray: The computed data.
        """
        self.data = self._Breakthrough.compute()
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
        ax.set_title(f"{self.DisplayName} ({plot_type})")
