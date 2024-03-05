import pytest
from ruptura import Components

def test_initialization_with_no_components():
    components = Components()
    assert len(components._components) == 0
    assert components.labels == []
    assert components.CarrierGas is None

def test_initialization_with_single_component():
    component = [{'MoleculeName': 'Nitrogen', 'GasPhaseMolFraction': 0.78}]
    components = Components(component)
    assert len(components._components) == 1
    assert components._components[0].MoleculeName == 'Nitrogen'
    assert components.CarrierGas is None

def test_initialization_with_multiple_components():
    component_list = [
        {'MoleculeName': 'Nitrogen', 'GasPhaseMolFraction': 0.78},
        {'MoleculeName': 'Oxygen', 'GasPhaseMolFraction': 0.21, 'CarrierGas': True}
    ]
    components = Components(component_list)
    assert len(components._components) == 2
    assert components.CarrierGas == 1

def test_add_component():
    components = Components()
    components.addComponent('Argon', 0.01)
    assert len(components._components) == 1
    assert components._components[0].MoleculeName == 'Argon'

def test_add_carrier_gas():
    components = Components()
    components.addComponent('Helium', 0.0005, CarrierGas=True)
    assert components.CarrierGas == 0

def test_get_labels():
    component_list = [
        {'MoleculeName': 'Nitrogen', 'GasPhaseMolFraction': 0.78},
        {'MoleculeName': 'Oxygen', 'GasPhaseMolFraction': 0.21}
    ]
    components = Components(component_list)
    labels = components.getLabels()
    assert len(labels) == 2
    assert 'Nitrogen' in labels[0]
    assert 'Oxygen' in labels[1]


def test_add_component_with_negative_mol_fraction():
    components = Components()
    with pytest.raises(ValueError):
        components.addComponent('Xenon', -0.1)
