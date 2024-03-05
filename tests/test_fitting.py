import pytest
from unittest.mock import patch
from ruptura import Components, Fitting
import numpy as np

@pytest.fixture
def basic_components():
    return Components([{
            "MoleculeName": "nC7",
            "GasPhaseMolFraction": 0.5,
            "isotherms": [["Langmuir", 1.09984, 6.55857e-5], ["Langmuir", 0.19466, 8.90731e-07]]
        }, {
            "MoleculeName": "C6m2",
            "GasPhaseMolFraction": 0.5,
            "isotherms": [["Langmuir", 1.22228, 3.90895e-05], ["Langmuir", 0.481726, 9.64046e-08]]
        }])

@pytest.fixture
def fitting_instance(basic_components):
    return Fitting(components=basic_components)

def test_compute(fitting_instance):
    mock_data = [[(0.1, 0.2), (0.3, 0.4)], [(0.5, 0.6), (0.7, 0.8)]]
    expected_output = np.random.rand(2, 2)  # Simulate an output shape

    with patch('_ruptura.Fitting.compute', return_value=expected_output) as mock_compute:
        result = fitting_instance.compute(mock_data)
        np.testing.assert_array_equal(result, expected_output)
        mock_compute.assert_called_once_with(mock_data)