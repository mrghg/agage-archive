from agage_archive.util import parse_fortran_format


def test_parse_fortran_format():

    format_string = "(F10.5, 2I4,I6, 2I4,I6,1X,70(F12.3,a1))"

    column_specs, column_types = parse_fortran_format(format_string)

    assert len(column_specs) == 7+2*70
    assert column_specs[0] == (0, 10)
    assert column_specs[1] == (10, 14)
    assert len(column_types) == 7+2*70
    assert column_types[0] == float

    format_string = "I4"
    column_specs, column_types = parse_fortran_format(format_string)
    assert len(column_specs) == 1
    assert column_specs[0] == (0, 4)
    assert len(column_types) == 1
    assert column_types[0] == int

    format_string = "10(F12.3)"
    column_specs, column_types = parse_fortran_format(format_string)
    assert len(column_specs) == 10
    assert column_specs[0] == (0, 12)
    assert len(column_types) == 10
    assert column_types[0] == float

    format_string = "(a1, 2I4, 2F10.5)"
    column_specs, column_types = parse_fortran_format(format_string)
    assert len(column_specs) == 5
    assert column_specs[0] == (0, 1)
    assert len(column_types) == 5
    assert column_types[0] == str
    assert column_types[1] == int
    assert column_types[2] == int
    assert column_types[3] == float

