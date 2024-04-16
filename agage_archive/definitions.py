

nc4_types = {"f4": "float32",
            "f8": "float64",
            "i4": "int32",
            "i2": "int16",
            "i8": "int64",
            "i1": "int8",
            "u1": "uint8",
            "u2": "uint16",
            "u4": "uint32",
            "u8": "uint64",
            "S1": "S1"}

scale_translator = {"TU1987": "TU-87"}

unit_translator = {"ppm": "1e-6",
                    "ppb": "1e-9",
                    "ppt": "1e-12",
                    "ppq": "1e-15",
                    "nmol/mol": "1e-9",
                    "nmol mol-1": "1e-9",
                    "pmol/mol": "1e-12",
                    "pmol mol-1": "1e-12",
                    }

species_translator = {"pfc-116": "c2f6",
                      "pfc-218": "c3f8",
                      "pfc-318": "c4f8"}

minimum_averaging_period = {"Picarro": "1H"}


def instrument_type_definition():
    '''Define instrument numbers for each instrument type

    Returns:
        str: Instrument type definition
    '''

    instrument_number = {"UNDEFINED": -1,
                        "ALE": 0,
                        "GAGE": 1,
                        "GCMD": 2,
                        "GCMS-ADS": 3,
                        "GCMS-Medusa": 4,
                        "GCECD": 5,
                        "GCTOFMS": 6,
                        "Picarro": 7,
                        "LGR": 8,
                        "GCMS-MteCimone": 9}
    
    # Create string from dictionary defining instrument numbers
    instrument_number_string = ", ".join([f"{k}={v}" for k, v in instrument_number.items()])

    return instrument_number, instrument_number_string
