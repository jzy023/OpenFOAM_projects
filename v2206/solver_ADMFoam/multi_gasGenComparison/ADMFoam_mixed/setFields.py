def write_setFieldsDict(filename: str, regions: list):
    header = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v2212                                 |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      setFieldsDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

"""

    default_field_values = """defaultFieldValues
(
    volScalarFieldValue alpha.liquid   1.0
    volScalarFieldValue alpha.gas      0.0
    volScalarFieldValue Gh2            0.0
    volScalarFieldValue Gch4           0.0
    volScalarFieldValue Gco2           0.0
);

"""

    regions_start = "regions\n(\n"

    regions_blocks = ""
    for region in regions:
        box = region.get("box", "(0 -0.013 0) (1 0.00 1)")
        states = region.get("states", {})

        block = f"""    boxToCell
    {{
        box {box};
        fieldValues
        (
"""

        block += "".join(
            f"            volScalarFieldValue {key:<15} {value}\n" for key, value in states.items()
        )

        block += """        );
    }
"""

        regions_blocks += block + "\n"

    regions_end = ");\n\n\n// ************************************************************************* //\n"

    with open(filename, "w") as f:
        f.write(header)
        f.write(default_field_values)
        f.write(regions_start)
        f.write(regions_blocks)
        f.write(regions_end)


def mix_states(states_a: dict, states_b: dict, weight_a: float, weight_b: float, normalize: bool = True) -> dict:
    """Return a weighted-average combination of two state dictionaries.

    - states_a, states_b: dictionaries mapping field name -> float
    - weight_a, weight_b: weights (e.g., layer thicknesses) for each state
    - normalize: if True, normalize weights so they sum to 1

    Keys present in only one of the dictionaries are treated as 0.0 in the other.
    The output key order follows states_a first, then any remaining keys from states_b.
    """
    wa, wb = float(weight_a), float(weight_b)
    s = wa + wb
    if normalize and s > 0:
        wa, wb = wa / s, wb / s

    # Preserve order from states_a, then append keys unique to states_b
    ordered_keys = list(states_a.keys()) + [k for k in states_b.keys() if k not in states_a]

    mixed = {}
    for k in ordered_keys:
        v1 = float(states_a.get(k, 0.0))
        v2 = float(states_b.get(k, 0.0))
        mixed[k] = v1 * wa + v2 * wb
    return mixed


if __name__ == "__main__":
    # import argparse
    statesGas = {
        'alpha.liquid': 0.0,
        'alpha.gas':    1.0,
        'Saa': 0.00,
        'Ssu': 0.00,
        'Sfa': 0.00,
        'Sva': 0.00,
        'Sbu': 0.00,
        'Spro': 0.00,
        'Sac': 0.00,
        'Sh2': 0.00,
        'Sch4': 0.00,
        'SIC': 0.00,
        'SIN': 0.00,
        'SI': 0.00,
        'Xc': 0.00,
        'Xch': 0.00,
        'Xpr': 0.00,
        'Xli': 0.00,
        'Xsu': 0.00,
        'Xaa': 0.00,
        'Xfa': 0.00,
        'Xc4': 0.00,
        'Xpro': 0.00,
        'Xac': 0.00,
        'Xh2': 0.00,
        'XI': 0.00,
        'Snh3': 0.00,
        'Gh2': 1.1032e-2,
        'Gch4': 1.6535,
        'Gco2': 13.500,
    }

    statesH00 = {
        'alpha.liquid': 1.0,
        'alpha.gas':    0.0,
        'Ssu': 0.012394,
        'Saa': 0.005500,
        'Sfa': 0.107400,
        'Sva': 0.012300,
        'Sbu': 0.014000,
        'Spro': 0.01760,
        'Sac': 0.089300,
        'Sh2': 2.5055e-7,
        'Sch4': 0.05550,
        'SIC': 95.100,
        'SIN': 94.500,
        'SI': 0.130900,
        'Xc': 0.107900,
        'Xch': 0.020500,
        'Xpr': 0.084200,
        'Xli': 0.043600,
        'Xsu': 0.312200,
        'Xaa': 0.93170,
        'Xfa': 0.338400,
        'Xc4': 0.32580,
        'Xpro': 0.10110,
        'Xac': 0.677200,
        'Xh2': 0.284800,
        'XI': 17.21620,
        'Snh3': 0.001884,
        'Gh2':  0.00,
        'Gch4': 0.00,
        'Gco2': 0.00,
    }

    ## Benchmark case without inflows -----------------------------------------
    statesH01 = {
        'alpha.liquid': 1.0,
        'alpha.gas':    0.0,
        'Ssu': 0.0111176,
        'Saa': 0.00389562,
        'Sfa': 0.104667,
        'Sva': 0.0116361,
        'Sbu': 0.0132247,
        'Spro': 0.0171919,
        'Sac': 0.0874941,
        'Sh2': 2.31122e-07,
        'Shp': 5.39349e-05,
        'Sch4': 0.0554709,
        'SIC': 95.2933,
        'SIN': 94.6716,
        'SI': 0.131125,
        'Xc': 0.108128,
        'Xch': 0.0138819,
        'Xpr': 0.0558735,
        'Xli': 0.0292937,
        'Xsu': 0.312849,
        'Xaa': 0.933353,
        'Xfa': 0.339136,
        'Xc4': 0.3265,
        'Xpro': 0.10132,
        'Xac': 0.678664,
        'Xh2': 0.285395,
        'XI': 17.2167,
        'Sco2': 9.38502,
        'Snh3': 1.90957,
        'Gh2':  0.00, # 1.06015e-05
        'Gch4': 0.00, # 1.65787
        'Gco2': 0.00, # 13.4228
    }

    statesH03 = {
        'alpha.liquid': 1.0,
        'alpha.gas':    0.0,
        'Ssu': 0.00656388,
        'Saa': 0.00171624,
        'Sfa': 0.0888668,
        'Sva': 0.00739683,
        'Sbu': 0.00834974,
        'Spro': 0.0137591,
        'Sac': 0.0719942,
        'Sh2': 1.72227e-07,
        'Shp': 5.28208e-05,
        'Sch4': 0.0544192,
        'SIC': 95.605,
        'SIN': 94.8489,
        'SI': 0.131577,
        'Xc': 0.108583,
        'Xch': 0.00664534,
        'Xpr': 0.0248931,
        'Xli': 0.0136491,
        'Xsu': 0.313681,
        'Xaa': 0.934521,
        'Xfa': 0.340486,
        'Xc4': 0.327518,
        'Xpro': 0.101715,
        'Xac': 0.681336,
        'Xh2': 0.286322,
        'XI': 17.2176,
        'Sco2': 9.24001,
        'Snh3': 1.95267,
        'Gh2':  0.00, # 8.38073e-06
        'Gch4': 0.00, # 1.6528
        'Gco2': 0.00, # 13.2264
    }

    statesH06 = {
        'alpha.liquid': 1.0,
        'alpha.gas':    0.0,
        'Ssu': 0.00254834,
        'Saa': 0.000539913,
        'Sfa': 0.0600253,
        'Sva': 0.0027227,
        'Sbu': 0.00305923,
        'Spro': 0.00684772,
        'Sac': 0.0403097,
        'Sh2': 9.6478e-08,
        'Shp': 5.08592e-05,
        'Sch4': 0.0520997,
        'SIC': 95.9491,
        'SIN': 94.9373,
        'SI': 0.132257,
        'Xc': 0.109247,
        'Xch': 0.0026813,
        'Xpr': 0.0079086,
        'Xli': 0.00507643,
        'Xsu': 0.313884,
        'Xaa': 0.933748,
        'Xfa': 0.341969,
        'Xc4': 0.327894,
        'Xpro': 0.102047,
        'Xac': 0.683983,
        'Xh2': 0.286993,
        'XI': 17.2189,
        'Sco2': 8.96106,
        'Snh3': 2.02826,
        'Gh2':  0.00, # 5.06568e-06
        'Gch4': 0.00, # 1.63892
        'Gco2': 0.00, # 12.8475
    }

    statesH12 = {
        'alpha.liquid': 1.0,
        'alpha.gas':    0.0,
        'Ssu': 0.000794575,
        'Saa': 0.000109709,
        'Sfa': 0.0225361,
        'Sva': 0.000394737,
        'Sbu': 0.00047003,
        'Spro': 0.00112412,
        'Sac': 0.00942798,
        'Sh2': 2.98875e-08,
        'Shp': 4.9022e-05,
        'Sch4': 0.0489446,
        'SIC': 96.2163,
        'SIN': 94.991,
        'SI': 0.133631,
        'Xc': 0.110459,
        'Xch': 0.00123065,
        'Xpr': 0.00165961,
        'Xli': 0.00193251,
        'Xsu': 0.312948,
        'Xaa': 0.929842,
        'Xfa': 0.342919,
        'Xc4': 0.326851,
        'Xpro': 0.101933,
        'Xac': 0.684397,
        'Xh2': 0.28675,
        'XI': 17.2217,
        'Sco2': 8.69072,
        'Snh3': 2.10378,
        'Gh2':  0.00, # 1.68704e-06
        'Gch4': 0.00, # 1.61264
        'Gco2': 0.00, # 12.4848
    }

    statesH24 = {
        'alpha.liquid': 1.0,
        'alpha.gas':    0.0,
        'Ssu': 0.000647386,
        'Saa': 7.3294e-05,
        'Sfa': 0.00481113,
        'Sva': 0.000170749,
        'Sbu': 0.000223943,
        'Spro': 0.000359969,
        'Sac': 0.00204735,
        'Sh2': 8.8505e-09,
        'Shp': 4.86248e-05,
        'Sch4': 0.0478863,
        'SIC': 96.3742,
        'SIN': 95.0936,
        'SI': 0.136417,
        'Xc': 0.112335,
        'Xch': 0.0011209,
        'Xpr': 0.00112379,
        'Xli': 0.00168193,
        'Xsu': 0.310457,
        'Xaa': 0.921079,
        'Xfa': 0.341051,
        'Xc4': 0.323836,
        'Xpro': 0.101051,
        'Xac': 0.679233,
        'Xh2': 0.284539,
        'XI': 17.2272,
        'Sco2': 8.64078,
        'Snh3': 2.12287,
        'Gh2':  0.00, # 4.75746e-07
        'Gch4': 0.00, # 1.60118
        'Gco2': 0.00, # 12.4208
    }

    ## Benchmark case with inflows --------------------------------------------
    # statesH01 = {
    #     'alpha.liquid': 1.0,
    #     'alpha.gas':    0.0,
    #     'Ssu': 0.0111176,
    #     'Saa': 0.00389562,
    #     'Sfa': 0.104667,
    #     'Sva': 0.0116361,
    #     'Sbu': 0.0132247,
    #     'Spro': 0.0171919,
    #     'Sac': 0.0874941,
    #     'Sh2': 2.31122e-07,
    #     'Shp': 5.39349e-05,
    #     'Sch4': 0.0554709,
    #     'SIC': 95.2933,
    #     'SIN': 94.6716,
    #     'SI': 0.131125,
    #     'Xc': 0.108128,
    #     'Xch': 0.0138819,
    #     'Xpr': 0.0558735,
    #     'Xli': 0.0292937,
    #     'Xsu': 0.312849,
    #     'Xaa': 0.933353,
    #     'Xfa': 0.339136,
    #     'Xc4': 0.3265,
    #     'Xpro': 0.10132,
    #     'Xac': 0.678664,
    #     'Xh2': 0.285395,
    #     'XI': 17.2167,
    #     'Sco2': 9.38502,
    #     'Snh3': 1.90957,
    #     'Gh2':  0.00, # 1.06015e-05
    #     'Gch4': 0.00, # 1.65787
    #     'Gco2': 0.00, # 13.4228
    # }

    # statesH03 = {
    #     'alpha.liquid': 1.0,
    #     'alpha.gas':    0.0,
    #     'Ssu': 0.00656388,
    #     'Saa': 0.00171624,
    #     'Sfa': 0.0888668,
    #     'Sva': 0.00739683,
    #     'Sbu': 0.00834974,
    #     'Spro': 0.0137591,
    #     'Sac': 0.0719942,
    #     'Sh2': 1.72227e-07,
    #     'Shp': 5.28208e-05,
    #     'Sch4': 0.0544192,
    #     'SIC': 95.605,
    #     'SIN': 94.8489,
    #     'SI': 0.131577,
    #     'Xc': 0.108583,
    #     'Xch': 0.00664534,
    #     'Xpr': 0.0248931,
    #     'Xli': 0.0136491,
    #     'Xsu': 0.313681,
    #     'Xaa': 0.934521,
    #     'Xfa': 0.340486,
    #     'Xc4': 0.327518,
    #     'Xpro': 0.101715,
    #     'Xac': 0.681336,
    #     'Xh2': 0.286322,
    #     'XI': 17.2176,
    #     'Sco2': 9.24001,
    #     'Snh3': 1.95267,
    #     'Gh2':  0.00, # 8.38073e-06
    #     'Gch4': 0.00, # 1.6528
    #     'Gco2': 0.00, # 13.2264
    # }

    # statesH06 = {
    #     'alpha.liquid': 1.0,
    #     'alpha.gas':    0.0,
    #     'Ssu': 0.00254834,
    #     'Saa': 0.000539913,
    #     'Sfa': 0.0600253,
    #     'Sva': 0.0027227,
    #     'Sbu': 0.00305923,
    #     'Spro': 0.00684772,
    #     'Sac': 0.0403097,
    #     'Sh2': 9.6478e-08,
    #     'Shp': 5.08592e-05,
    #     'Sch4': 0.0520997,
    #     'SIC': 95.9491,
    #     'SIN': 94.9373,
    #     'SI': 0.132257,
    #     'Xc': 0.109247,
    #     'Xch': 0.0026813,
    #     'Xpr': 0.0079086,
    #     'Xli': 0.00507643,
    #     'Xsu': 0.313884,
    #     'Xaa': 0.933748,
    #     'Xfa': 0.341969,
    #     'Xc4': 0.327894,
    #     'Xpro': 0.102047,
    #     'Xac': 0.683983,
    #     'Xh2': 0.286993,
    #     'XI': 17.2189,
    #     'Sco2': 8.96106,
    #     'Snh3': 2.02826,
    #     'Gh2':  0.00, # 5.06568e-06
    #     'Gch4': 0.00, # 1.63892
    #     'Gco2': 0.00, # 12.8475
    # }

    # statesH12 = {
    #     'alpha.liquid': 1.0,
    #     'alpha.gas':    0.0,
    #     'Ssu': 0.000794575,
    #     'Saa': 0.000109709,
    #     'Sfa': 0.0225361,
    #     'Sva': 0.000394737,
    #     'Sbu': 0.00047003,
    #     'Spro': 0.00112412,
    #     'Sac': 0.00942798,
    #     'Sh2': 2.98875e-08,
    #     'Shp': 4.9022e-05,
    #     'Sch4': 0.0489446,
    #     'SIC': 96.2163,
    #     'SIN': 94.991,
    #     'SI': 0.133631,
    #     'Xc': 0.110459,
    #     'Xch': 0.00123065,
    #     'Xpr': 0.00165961,
    #     'Xli': 0.00193251,
    #     'Xsu': 0.312948,
    #     'Xaa': 0.929842,
    #     'Xfa': 0.342919,
    #     'Xc4': 0.326851,
    #     'Xpro': 0.101933,
    #     'Xac': 0.684397,
    #     'Xh2': 0.28675,
    #     'XI': 17.2217,
    #     'Sco2': 8.69072,
    #     'Snh3': 2.10378,
    #     'Gh2':  0.00, # 1.68704e-06
    #     'Gch4': 0.00, # 1.61264
    #     'Gco2': 0.00, # 12.4848
    # }

    # statesH24 = {
    #     'alpha.liquid': 1.0,
    #     'alpha.gas':    0.0,
    #     'Ssu': 0.000647386,
    #     'Saa': 7.3294e-05,
    #     'Sfa': 0.00481113,
    #     'Sva': 0.000170749,
    #     'Sbu': 0.000223943,
    #     'Spro': 0.000359969,
    #     'Sac': 0.00204735,
    #     'Sh2': 8.8505e-09,
    #     'Shp': 4.86248e-05,
    #     'Sch4': 0.0478863,
    #     'SIC': 96.3742,
    #     'SIN': 95.0936,
    #     'SI': 0.136417,
    #     'Xc': 0.112335,
    #     'Xch': 0.0011209,
    #     'Xpr': 0.00112379,
    #     'Xli': 0.00168193,
    #     'Xsu': 0.310457,
    #     'Xaa': 0.921079,
    #     'Xfa': 0.341051,
    #     'Xc4': 0.323836,
    #     'Xpro': 0.101051,
    #     'Xac': 0.679233,
    #     'Xh2': 0.284539,
    #     'XI': 17.2272,
    #     'Sco2': 8.64078,
    #     'Snh3': 2.12287,
    #     'Gh2':  0.00, # 4.75746e-07
    #     'Gch4': 0.00, # 1.60118
    #     'Gco2': 0.00, # 12.4208
    # }

    # states1 = {k: 0.0 for k in states0.keys()}  # example second region with zero values

    # # Parse user-provided layer thicknesses to compute weights for mixing
    # parser = argparse.ArgumentParser(description="Generate setFieldsDict with optional mixed initial states.")
    # parser.add_argument('--hH00', type=float, default=1.0, help='Layer thickness (weight) for statesH00')
    # parser.add_argument('--hH12', type=float, default=1.0, help='Layer thickness (weight) for statesH12')
    # args = parser.parse_args()

    # Create mixed state using thickness-proportional weights
    # statesMixed = mix_states(statesH00, statesH12, args.hH00, args.hH12, normalize=True)

    statesWeight1 = (0.08 - 0.013)/(0.16 - 0.013)  # weight for statesH00 in the mixed region
    statesWeight2 = (0.16 - 0.08)/(0.16 - 0.013)  # weight for statesH12 in the mixed region
    statesMixed = mix_states(statesH00, statesH12, statesWeight1, statesWeight2, normalize=True)

    regions = [
        # {"box": "(0 -0.013 0) (1 -0.000 1)", "states": statesGas},
        # {"box": "(0 -0.080 0) (1 -0.013 1)", "states": statesH00},
        # {"box": "(0 -0.160 0) (1 -0.080 1)", "states": statesH00},

        # {"box": "(0 -0.160 0) (1 0.000 1)", "states": statesH00},
        {"box": "(0 -0.160 0) (1 0.000 1)", "states": statesMixed},
    ]

    write_setFieldsDict("./system/setFieldsDict", regions)
