import numpy as np
import pandas as pd

# Initial State ===============================================================================
names = ['Ssu', 'Saa', 'Sfa', 'Sva', 'Sbu', 'Spro', 'Sac', 'Sh2', 'Sch4', 'SIC', 'SIN', 'SI',
         'Xc',  'Xch', 'Xpr', 'Xli', 'Xsu', 'Xaa',  'Xfa', 'Xc4', 'Xpro', 'Xac', 'Xh2', 'XI',
         'Snh3','Gh2', 'Gch4','Gco2']

values = [0.012394, 0.005500, 0.107400, 0.012300, 0.014000, 0.01760, 0.089300, 2.5055e-7, 0.05550, 0.095100, 0.094500, 0.130900,
          0.107900, 0.020500, 0.084200, 0.043600, 0.312200, 0.93170, 0.338400, 0.32580,   0.10110, 0.677200, 0.284800, 17.21620, 
          0.001884, 1.1032e-5, 1.65350, 0.013500]


# Other State   ===============================================================================



# Headers =====================================================================================
header_0 = """
/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\    /   O peration     | Version:  v2206                                 |
|   \\\  /    A nd           | Website:  www.openfoam.com                      |
|    \\\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;\n"""


header_1 = """
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n
dimensions      [1 -3 0 0 0 0 0];\n\n"""


content = """
boundaryField
{
    allBoundary
    {
        type            zeroGradient;
    }

    top
    {
        type            zeroGradient;
    }

    GeoTank
    {
        type            zeroGradient;
    }

    GeoShaft
    {
        type            zeroGradient;
    }
}\n\n
// ************************************************************************* //\n"""




for i, name in enumerate(names):
    header_name = f'    object      {name};\n'
    content_value = f'internalField\tuniform {values[i]};\n'
    with open(f'./0.orig/{name}', 'w') as file:
        file.write(header_0)
        file.write(header_name)
        file.write(header_1)
        file.write(content_value)
        file.write(content)
