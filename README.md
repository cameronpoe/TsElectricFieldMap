# TsElectricFieldMap: A Custom TOPAS Extension for Non-uniform Electric Fields

Created by Cameron Poe (cameronpoe@uchicago.edu) and Kepler Domurat-Sousa (domuratsousa@uchicago.edu)

Enrico Fermi Institute, University of Chicago

3 October 2023

## About the extension

This custom TOPAS extension modifies TsMagneticFieldMap to create non-uniform electric fields from a modified Opera-3d `.TABLE` file containing electric field data.

This package also includes a Python script to turn Ansys `.fld`-formatted files into the modified Opera-3d `.TABLE`-formatted files required by the extension.

## How to use the extension

To install the extension, follow the instructions on the TOPAS installation guide at https://www.topasmc.org/user-guides.

The `efield_table_creator.py` script takes as inputs the path of your Ansys `.fld` file and then creates an Opera-3d `.TABLE` file in the `non-uniform_efield_tables/` directory.

To use the `efield_table_creator.py` script, scroll to the `__main__` section of the code. In the General Parameters section, set `MODE` to `non-uniform` and set `output_file_name` to what you want the `.TABLE` file to be called. Do not include `.TABLE` in that variable; it will be added by the code.

In the Non-uniform Mode Parameters section, set `ansys_file_nonuniform` to the path of the Ansys `.fld` file you wish to convert. This Python script uses mm and kV/mm as the intrinsic units, so set `ansys_position_scale_nonuniform` to the value that turns the `.fld` position units to mm. 

This script can process Ansys simulations that have symmetry in one dimension. If this describes your `.fld` file, set `repetition_axis` to the axis of translational symmetry (e.g. `'x'` if your simulation is in y and z). Set `repetition_depth` to the extent in mm of the object containing your field in the symmetric axis. If your simulation does not have translational symmetry in one dimension, set these values to `''` and `0`, respectively.

The Python script also has the ability to run some simple checks. These are included in the `uniform_tests/`, `linear_tests/`, and `nonuniform_tests/` directories, but you can also re-create these files using the Python script. To do so, set the `MODE` in the General Parameters section appropriately, then modify any parameters in the respective Test Parameters section.

## Extension next steps

This extension can still be improved, but our research group does not have incentive since our final detector design does not use non-uniform electric fields. 

The biggest improvement would be to translate the Python `.fld` to `.TABLE` converter script into the C++ of the extension. This would involve a couple things. First, the extension would need to be modified to parse a `.fld` file instead of a `.TABLE` file. The header information is the biggest difference, so instead of retrieving `fNX`, `fNY`, and `fNZ` in the first line of the `.TABLE` file, the extension would need to calculate those quantities from the grid min, max, and size arrays in the `.fld` file. Second, the indexing behavior could be changed. This is optional. The current extension assumes that rows in the `.TABLE` file are sorted least-to-greatest first in X, then in Y, then in Z. While this is how `.fld` files are sorted, a more general solution would be to calculate `ix`, `iy`, and `iz` using the grid min, max, and size arrays and the position vector in a row of the `.fld` file. These index values would then no longer rely on the data being sorted in a specific order.