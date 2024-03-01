import pytest
from unittest.mock import patch
from ruptura import Components, MixturePrediction
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
def mixture_prediction(basic_components):
    return MixturePrediction(components=basic_components)

def test_initialization(mixture_prediction):
    assert mixture_prediction.shape == (100, 2, 6)
    assert mixture_prediction.DisplayName == "Column"
    assert mixture_prediction._MixturePrediction is not None
    assert mixture_prediction.data is None

def test_compute(mixture_prediction):
    with patch('_ruptura.MixturePrediction.compute') as mock_compute:
        mock_data = np.random.rand(100, 2, 6)  # Mocking the computation result
        mock_compute.return_value = mock_data

        result = mixture_prediction.compute()
        np.testing.assert_array_equal(result, mock_data)
        mock_compute.assert_called_once()