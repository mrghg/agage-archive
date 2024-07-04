import pytest

from agage_archive.definitions import get_instrument_type, get_instrument_number


def test_get_instrument_type():
    
    # Test for single instrument number
    instrument_number = 1
    instrument_type = get_instrument_type(instrument_number)
    assert instrument_type == "GAGE"

    # Test for list of instrument numbers
    instrument_numbers = [1, 2]
    instrument_types = get_instrument_type(instrument_numbers)
    assert instrument_types == ["GAGE", "GCMD"]
    
    # Test for invalid input
    instrument_number = "1"
    with pytest.raises(ValueError):
        instrument_type = get_instrument_type(instrument_number)


def test_get_instrument_number():
    '''Test get_instrument_number function'''

    # Test for single instrument type
    instrument = "GAGE"
    instrument_number = get_instrument_number(instrument)
    assert instrument_number == 1

    # Test for invalid input
    instrument = 1
    with pytest.raises(ValueError):
        instrument_number = get_instrument_number(instrument)

    # Test for partial match
    instrument = "Picarro-1"
    instrument_number = get_instrument_number(instrument)
    assert instrument_number == 7

    # Test for no match
    instrument = "Invalid"
    with pytest.raises(KeyError):
        instrument_number = get_instrument_number(instrument)

