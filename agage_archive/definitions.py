

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

species_translator_flask = {"c2f6": "PFC-116",
                            "c3f8": "PFC-218",
                            "c4f8": "PFC-318",
                            "benzene": "benzene",
                            "hfc-134a": "HFC-134a",
                            "hfc-152a": "HFC-152a",
                            "hfc-143a": "HFC-143a",
                            "hfc-227ea": "HFC-227ea",
                            "hfc-236fa": "HFC-236fa",
                            "hfc-245fa": "HFC-245fa",
                            "hfc-365mfc": "HFC-365mfc",
                            "hfc-4310mee": "HFC-4310mee",
                            "hcfc-22": "HCFC-22",
                            "hcfc-141b": "HCFC-141b",
                            "hcfc-142b": "HCFC-142b",
                            "hcfc-132b": "HCFC-132b",
                            "hcfc-133a": "HCFC-133a",
                            "ch3cl": "CH3Cl",
                            "ch3br": "CH3Br",
                            "ch2cl2": "CH2Cl2",
                            "chcl3": "CHCl3",
                            "ch3ccl3": "CH3CCl3",
                            "ccl4": "CCl4"
                            }

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
