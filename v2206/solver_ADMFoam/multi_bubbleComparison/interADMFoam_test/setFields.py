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
            volScalarFieldValue alpha.liquid    0.0
            volScalarFieldValue alpha.gas       1.0

            // TODO: make a sperate file to initialize
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


if __name__ == "__main__":
    states0 = {
        'Ssu': 0.012394,
        'Saa': 0.005500,
        'Sfa': 0.107400,
        'Sva': 0.012300,
        'Sbu': 0.014000,
        'Spro': 0.01760,
        'Sac': 0.089300,
        'Sh2': 2.5055e-7,
        'Sch4': 0.05550,
        'SIC': 0.095100,
        'SIN': 0.094500,
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
        'Gh2': 1.1032e-2,
        'Gch4': 1653.50,
        'Gco2': 13.500,
    }

    states1 = {
        'Ssu': 0.012394,
        'Saa': 0.005500,
        'Sfa': 0.107400,
        'Sva': 0.012300,
        'Sbu': 0.014000,
        'Spro': 0.01760,
        'Sac': 0.089300,
        'Sh2': 2.5055e-7,
        'Sch4': 0.05550,
        'SIC': 0.095100,
        'SIN': 0.094500,
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
        'Gh2': 1.1032e-2,
        'Gch4': 1653.50,
        'Gco2': 13.500,
    }

    # states1 = {k: 0.0 for k in states0.keys()}  # example second region with zero values

    regions = [
        {"box": "(0 -0.08 0) (1 -0.013 1)", "states": states0},
        {"box": "(0 -0.16 0) (1 -0.080 1)", "states": states1},
    ]

    write_setFieldsDict("./system/setFieldsDict", regions)
