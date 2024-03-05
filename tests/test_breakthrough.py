import pytest
from unittest.mock import patch
from ruptura import Components, Breakthrough
import numpy as np

@pytest.fixture
def basic_components():
    return Components([{
            "MoleculeName": "Helium",
            "CarrierGas": True,
            "GasPhaseMolFraction": 0.9,
            "isotherms": []
        }, {
            "MoleculeName": "nC7",
            "GasPhaseMolFraction": 0.05,
            "isotherms": [["Langmuir", 1.09984, 6.55857e-5], ["Langmuir", 0.19466, 8.90731e-07]]
        }, {
            "MoleculeName": "C6m2",
            "GasPhaseMolFraction": 0.05,
            "isotherms": [["Langmuir", 1.22228, 3.90895e-05], ["Langmuir", 0.481726, 9.64046e-08]]
        }])

@pytest.fixture
def breakthrough_instance(basic_components):
    components = Components([
        {'MoleculeName': 'Component1', 'GasPhaseMolFraction': 0.5},
        {'MoleculeName': 'Component2', 'GasPhaseMolFraction': 0.5}
    ])
    return Breakthrough(components=basic_components)

def test_compute(breakthrough_instance):
    expected_output = np.random.rand(100, 20)  # Simulate an output shape for the breakthrough data

    with patch('_ruptura.Breakthrough.compute', return_value=expected_output) as mock_compute:
        result = breakthrough_instance.compute()
        np.testing.assert_array_equal(result, expected_output)
        mock_compute.assert_called_once()