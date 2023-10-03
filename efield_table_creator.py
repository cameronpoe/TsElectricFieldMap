import os
import re
import numpy as np

def create_opera_3d_file_uniform(table_path, efield, points):

    # opens the .TABLE output in the directory this script is in
    file_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), table_path)
    f = open(file_path, 'w', newline='\n')

    # writes the Opera-3d header data
    f.write(f'''
    \t\t\t{len(points[0])}\t\t\t{len(points[1])}\t\t\t{len(points[0])}
1 X [MM]
2 Y [MM]
3 Z [MM]
4 EX [KV/MM]
5 EY [KV/MM]
6 EZ [KV/MM]
0
''')
    
    for x_val in points[0]:
        for y_val in points[1]:
            for z_val in points[2]:
                new_line = f'\t{x_val}\t{y_val}\t{z_val}\t{efield[0]}\t{efield[1]}\t{efield[2]}\n'
                f.write(new_line)                

    f.close()

    return

def create_opera_3d_file_linear(table_path, efield_dir, efield_grad, points):

    # opens the .TABLE output in the directory this script is in
    file_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), table_path)
    f = open(file_path, 'w', newline='\n')

    # writes the Opera-3d header data
    f.write(f'''
    \t\t\t{len(points[0])}\t\t\t{len(points[1])}\t\t\t{len(points[0])}
1 X [MM]
2 Y [MM]
3 Z [MM]
4 EX [KV/MM]
5 EY [KV/MM]
6 EZ [KV/MM]
0
''')
    
    direc_length = np.linalg.norm(np.array(efield_dir))
    dir_vec = np.array(efield_dir/direc_length)
    distance_to_spacing_ratios = np.array([points[0][1]-points[0][0], points[1][1]-points[1][0], points[2][1]-points[2][0]])

    for i, x_val in enumerate(points[0]):
        for j, y_val in enumerate(points[1]):
            for k, z_val in enumerate(points[2]):

                dist_from_zero_point = np.array([i,j,k])*distance_to_spacing_ratios
                dist_along_dir_vec = np.dot(dist_from_zero_point, dir_vec)
                this_e = dist_along_dir_vec*dir_vec*efield_grad

                new_line = f'\t{x_val}\t{y_val}\t{z_val}\t{this_e[0]}\t{this_e[1]}\t{this_e[2]}\n'
                f.write(new_line)              
                

    f.close()

    return

def create_opera_3d_file_nonuniform(table_path, ansys_path, repeat_axis, repeat_depth, ansys_pos_scale):

    axis_dict = {
        'x': 0,
        'y': 1,
        'z': 2
    }

    table_path = r'non-uniform_efield_tables/' + table_path + r'.TABLE'
    f = open_file_location(table_path)
    simulation_file = open(ansys_path, 'r')

    # takes the first line ofthe .fld file and puts each bracketed block of information into an element of grid_list
    first_line = simulation_file.readline()
    grid_list = re.findall(r'(\[.+?\])', first_line)

    # initializes a list to serve as the nested list full of formatted information taken from the first line of the .fld file
    grid_list_formatted = []
    grid_list_formatted_count = 0

    # loops through and formats each bracketed piece of information in the first line of the .fld file
    for bracketed_unit in grid_list:

        bracketed_unit = bracketed_unit[1:-1]   # strips the string of brackets []
        components = bracketed_unit.split()     # splits the bracketed information in x, y, and z components
        grid_list_formatted.append([])          # adds a new list to grid_list_formatted forthis bracketed_unit

        for component in components:
            value = re.search(r'[\-\.0-9]+', component).group(0)
            grid_list_formatted[grid_list_formatted_count].append(float(value))

        grid_list_formatted_count += 1
    
    midpoint = np.array([
        0.5*(grid_list_formatted[1][0]+grid_list_formatted[0][0]),
        0.5*(grid_list_formatted[1][1]+grid_list_formatted[0][1]),
        0.5*(grid_list_formatted[1][2]+grid_list_formatted[0][2]),
        ])

    if repeat_axis:
        repeat_index = axis_dict[repeat_axis.lower()]
        grid_list_formatted[0][repeat_index] = -repeat_depth/2
        grid_list_formatted[1][repeat_index] = repeat_depth/2
        grid_list_formatted[2][repeat_index] = repeat_depth

    unique_points = [
        int((grid_list_formatted[1][0] - grid_list_formatted[0][0])/grid_list_formatted[2][0]) + 1,
        int((grid_list_formatted[1][1] - grid_list_formatted[0][1])/grid_list_formatted[2][1]) + 1,
        int((grid_list_formatted[1][2] - grid_list_formatted[0][2])/grid_list_formatted[2][2]) + 1
    ]

    # writes the Opera3D header data
    f.write(f'''
    \t\t\t{unique_points[0]}\t\t\t{unique_points[1]}\t\t\t{unique_points[2]}
1 X [MM]
2 Y [MM]
3 Z [MM]
4 EX [KV/MM]
5 EY [KV/MM]
6 EZ [KV/MM]
7 BMOD/HMOD [Wb] 
0
''')

    simulation_file.readline()

    for line in simulation_file.readlines():

        values = line.split()
        
        if repeat_axis:
            repeat_index = axis_dict[repeat_axis.lower()]
            values[repeat_index] = (-repeat_depth/2)/ansys_pos_scale

        pos_vec = ansys_pos_scale*np.array([float(values[0]), float(values[1]), float(values[2])])
        Ex, Ey, Ez = 1e-6*np.array([float(values[3]), float(values[4]), float(values[5])])

        pos_vec -= midpoint

        new_line = f'\t{pos_vec[0]}\t{pos_vec[1]}\t{pos_vec[2]}\t{Ex}\t{Ey}\t{Ez}\t0\n'
        f.write(new_line)
    simulation_file.close()
    
    if repeat_axis:

        simulation_file_again = open(ansys_path, 'r')
        simulation_file_again.readline()
        simulation_file_again.readline()

        repeat_index = axis_dict[repeat_axis.lower()]

        for line in simulation_file_again.readlines():

            values = line.split()
            
            values[repeat_index] = (repeat_depth/2)/ansys_pos_scale

            pos_vec = ansys_pos_scale*np.array([float(values[0]), float(values[1]), float(values[2])])
            Ex, Ey, Ez = 1e-6*np.array([float(values[3]), float(values[4]), float(values[5])])

            pos_vec -= midpoint

            new_line = f'\t{pos_vec[0]}\t{pos_vec[1]}\t{pos_vec[2]}\t{Ex}\t{Ey}\t{Ez}\t0\n'
            f.write(new_line)
        simulation_file_again.close()

    f.close()

    return

def open_file_location(file_name):

    file_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), file_name)
    f = open(file_path, 'w', newline='\n')

    return f

def write_topas_header(file_object):
    file_object.write(f'''sv:Ph/Default/Modules                   = 1 "g4em-standard_opt4"
b:Ts/PauseBeforeQuit    = "True"
i:Ts/NumberOfThreads    = -1
i:Ts/Seed = 2
b:Ts/UseQT  = "True"

s:Gr/ViewA/Type             = "OpenGL"
i:Gr/ViewA/WindowSizeX      = 1024
i:Gr/ViewA/WindowSizeY      = 768
d:Gr/ViewA/Theta            = 0 deg
d:Gr/ViewA/Phi              = 90 deg
s:Gr/ViewA/Projection       = "Orthogonal"
d:Gr/ViewA/PerspectiveAngle = 30 deg
#u:Gr/ViewA/Zoom             = 10000.
b:Gr/ViewA/HiddenLineRemovalForTrajectories = "True"

s:Gr/ViewB/Type             = "OpenGL"
i:Gr/ViewB/WindowSizeX      = 1024
i:Gr/ViewB/WindowSizeY      = 768
d:Gr/ViewB/Theta            = 75 deg
d:Gr/ViewB/Phi              = 10 deg
s:Gr/ViewB/Projection       = "Perspective"
d:Gr/ViewB/PerspectiveAngle = 30 deg
u:Gr/ViewB/Zoom             = 3.
b:Gr/ViewB/HiddenLineRemovalForTrajectories = "True"
  
''')
    return

def write_topas_world_lines(file_object, world_hls):
    file_object.write(f'''d:Ge/World/HLX          =  {world_hls[0]} mm
d:Ge/World/HLY          = {world_hls[1]} mm
d:Ge/World/HLZ          = {world_hls[2]} mm
s:Ge/World/Material     = "Vacuum"
b:Ge/World/Invisible    = "False"
                      
''')
    return

def write_topas_tsbox(file_object, box_dims, box_place, box_rot):
    file_object.write(f'''s:Ge/EnclosingBox/Type          = "TsBox"
s:Ge/EnclosingBox/Parent        = "World"
d:Ge/EnclosingBox/HLX           = {0.5*box_dims[0]} mm
d:Ge/EnclosingBox/HLY           = {0.5*box_dims[1]} mm
d:Ge/EnclosingBox/HLZ           = {0.5*box_dims[2]} mm
d:Ge/EnclosingBox/TransX        = {box_place[0]} mm
d:Ge/EnclosingBox/TransY        = {box_place[1]} mm
d:Ge/EnclosingBox/TransZ        = {box_place[2]} mm
d:Ge/EnclosingBox/RotX          = {box_rot[0]} deg
d:Ge/EnclosingBox/RotY          = {box_rot[1]} deg
d:Ge/EnclosingBox/RotZ          = {box_rot[2]} deg
s:Ge/EnclosingBox/Material      = "Vacuum"
s:Ge/EnclosingBox/Color             = "Yellow"

''')
    return

def create_topas_files_uniform(writing_to_file_name, box_place, box_rot, box_dims, efield_vec, efield_points):

    rel_directory_uniform = r'topas_run_files/uniform_runs/'
    native_filename_uniform = rel_directory_uniform + writing_to_file_name + r'_native.topas'  
    etable_filename_uniform = rel_directory_uniform + writing_to_file_name + r'_etable.topas' 
    f_native = open_file_location(native_filename_uniform)
    f_etable = open_file_location(etable_filename_uniform)
    world_halflengths_uniform = 1.5*(np.absolute(box_place) + 0.5*(box_dims))

    # Writing the TOPAS file for the native uniform electric field
    write_topas_header(f_native)
    write_topas_world_lines(f_native, world_halflengths_uniform)
    write_topas_tsbox(f_native, box_dims, box_place, box_rot)
    emag = np.linalg.norm(efield_vec)
    edirections = efield_vec/emag
    f_native.write(f'''s:Ge/EnclosingBox/Field = "UniformElectromagnetic"
u:Ge/EnclosingBox/ElectricFieldDirectionX = {edirections[0]}
u:Ge/EnclosingBox/ElectricFieldDirectionY = {edirections[1]}
u:Ge/EnclosingBox/ElectricFieldDirectionZ = {edirections[2]}
d:Ge/EnclosingBox/ElectricFieldStrength   = {emag} kV/mm
u:Ge/EnclosingBox/MagneticFieldDirectionX = 0.0
u:Ge/EnclosingBox/MagneticFieldDirectionY = 0.0
u:Ge/EnclosingBox/MagneticFieldDirectionZ = 0.0
d:Ge/EnclosingBox/MagneticFieldStrength   = 0.0 tesla

''')
    
    # Creating the Opera-3d-formatted table file for the electric field
    efield_table_path = rel_directory_uniform + writing_to_file_name + r'.TABLE'
    x_pos = np.linspace(-box_dims[0]/2, box_dims[0]/2, efield_points[0])
    y_pos = np.linspace(-box_dims[1]/2, box_dims[1]/2, efield_points[1])
    z_pos = np.linspace(-box_dims[2]/2, box_dims[2]/2, efield_points[2])
    efield_points = [x_pos, y_pos, z_pos]
    create_opera_3d_file_uniform(efield_table_path, efield_vec, efield_points)
    
    # Writing the TOPAS file for imported uniform electric field 
    write_topas_header(f_etable)
    write_topas_world_lines(f_etable, world_halflengths_uniform)
    write_topas_tsbox(f_etable, box_dims, box_place, box_rot)
    f_etable.write(f'''s:Ge/EnclosingBox/Field         = "TsElectricFieldMap"
s:Ge/EnclosingBox/ElectricField3DTable  = "{writing_to_file_name + r'.TABLE'}"

''')

    f_native.close()
    f_etable.close()
    return

def create_topas_file_linear(writing_to_file_name, box_place, box_rot, box_dims, efield_dir, efield_grad, efield_points):

    rel_directory_uniform = r'topas_run_files/linear_runs/'
    etable_filename_linear = rel_directory_uniform + writing_to_file_name + r'_etable.topas' 
    f_etable = open_file_location(etable_filename_linear)
    world_halflengths_uniform = 1.5*(np.absolute(box_place) + 0.5*(box_dims))
    
    # Creating the Opera-3d-formatted table file for the electric field
    efield_table_path = rel_directory_uniform + writing_to_file_name + r'.TABLE'
    x_pos = np.linspace(-box_dims[0]/2, box_dims[0]/2, efield_points[0])
    y_pos = np.linspace(-box_dims[1]/2, box_dims[1]/2, efield_points[1])
    z_pos = np.linspace(-box_dims[2]/2, box_dims[2]/2, efield_points[2])
    efield_points = [x_pos, y_pos, z_pos]
    create_opera_3d_file_linear(efield_table_path, efield_dir, efield_grad, efield_points)
    
    # Writing the TOPAS file for imported linear electric field 
    write_topas_header(f_etable)
    write_topas_world_lines(f_etable, world_halflengths_uniform)
    write_topas_tsbox(f_etable, box_dims, box_place, box_rot)
    f_etable.write(f'''s:Ge/EnclosingBox/Field         = "TsElectricFieldMap"
s:Ge/EnclosingBox/ElectricField3DTable  = "{writing_to_file_name + r'.TABLE'}"

''')

    f_etable.close()
    return

def create_topas_file_wire_dynodes_nonuniform(output_name, depth):

    output_path = r'topas_run_files/nonuniform_tests/' + output_name + r'.topas'
    table_path = r'../../non-uniform_efield_tables/' + output_name + r'.TABLE'

    # Calculating some constants
    R = 1                       # mm
    W = 4                       # mm
    H = (np.sqrt(3)/2)*W        # mm
    plane_thickness = 1         # mm
    plane_dist = 2              # mm
    vert_layers = 4
    horizontal_rods = 8
    plane_min_x = -10           # mm
    plane_max_x = (horizontal_rods - 0.5)*W + 2*R + 10                          # mm
    top_plane_y = 2*R + (vert_layers-1)*H + plane_dist + 0.5*plane_thickness    # mm   
    bot_plane_y = -plane_dist - 0.5*plane_thickness                             # mm
    plane_hlx = 0.5*(plane_max_x - plane_min_x)                                 # mm
    diagram_box_pos =  np.array([16, 6, depth/2])
    plane_transx = 0.5*(plane_max_x + plane_min_x) - diagram_box_pos[0]         # mm

    dynode_box_dims = np.array([64, 24, depth])     # mm
    dynode_box_place = np.array([0,0,0])            # mm
    dynode_box_rot = np.array([0,0,0])              # deg
    world_halflengths_dynode = 1.5*(np.absolute(dynode_box_place) + 0.5*(dynode_box_dims))

    f = open_file_location(output_path)
    write_topas_header(f)
    write_topas_world_lines(f, world_halflengths_dynode)
    write_topas_tsbox(f, dynode_box_dims, dynode_box_place, dynode_box_rot)

    f.write(f'''s:Ge/EnclosingBox/Field                  = "TsElectricFieldMap"
s:Ge/EnclosingBox/ElectricField3DTable   = "{table_path}"
              
''')

    f.write(f'''s:Ge/PlaneUpper/Parent            = "EnclosingBox"
s:Ge/PlaneUpper/Type              = "TsBox"
s:Ge/PlaneUpper/Material          = "G4_W" 
d:Ge/PlaneUpper/HLX               = {plane_hlx} mm
d:Ge/PlaneUpper/HLY               = {0.5*plane_thickness} mm
d:Ge/PlaneUpper/HLZ               = {0.5*depth} mm
d:Ge/PlaneUpper/TransX            = {plane_transx} mm
d:Ge/PlaneUpper/TransY            = {top_plane_y-diagram_box_pos[1]} mm
d:Ge/PlaneUpper/TransZ            = 0 mm
s:Ge/PlaneUpper/DrawingStyle = "solid"

s:Ge/PlaneLower/Parent            = "EnclosingBox"
s:Ge/PlaneLower/Type              = "TsBox"
s:Ge/PlaneLower/Material          = "G4_W" 
d:Ge/PlaneLower/HLX               = {plane_hlx} mm
d:Ge/PlaneLower/HLY               = {0.5*plane_thickness} mm
d:Ge/PlaneLower/HLZ               = {0.5*depth} mm
d:Ge/PlaneLower/TransX            = {plane_transx} mm
d:Ge/PlaneLower/TransY            = {bot_plane_y-diagram_box_pos[1]} mm
d:Ge/PlaneLower/TransZ            = 0 mm
s:Ge/PlaneLower/DrawingStyle = "solid"

''')

    i = 0
    while i < vert_layers:
        j = 0
        y_coord = R + i*H
        while j < horizontal_rods:
            if i%2 == 0:
                x_coord = (j+0.5)*W + R
            else:
                x_coord = j*W + R
            f.write(f'''s:Ge/Rod_{i+1}_{j+1}/Parent            = "EnclosingBox"
s:Ge/Rod_{i+1}_{j+1}/Type              = "TsCylinder"
s:Ge/Rod_{i+1}_{j+1}/Material          = "G4_W" 
d:Ge/Rod_{i+1}_{j+1}/Rmax              = {R} mm
d:Ge/Rod_{i+1}_{j+1}/HL                = {0.5*depth} mm
d:Ge/Rod_{i+1}_{j+1}/TransX            = {x_coord-diagram_box_pos[0]} mm
d:Ge/Rod_{i+1}_{j+1}/TransY            = {y_coord-diagram_box_pos[1]} mm
d:Ge/Rod_{i+1}_{j+1}/TransZ            = 0 mm
s:Ge/Rod_{i+1}_{j+1}/DrawingStyle = "solid"

''')

            j += 1
        i += 1
    f.close()
    return

def create_topas_file_sphere_plates_nonuniform(output_name):

    output_path = r'topas_run_files/nonuniform_tests/' + output_name + r'.topas'
    table_path = r'../../non-uniform_efield_tables/' + output_name + r'.TABLE'

    containing_box_dims = np.array([15, 15, 15])        # mm
    containing_box_place = np.array([0,0,0])            # mm
    containing_box_rot = np.array([0,0,0])              # deg
    world_halflengths_sphere_plates = 1.5*(np.absolute(containing_box_place) + 0.5*(containing_box_dims))

    f = open_file_location(output_path)
    write_topas_header(f)
    write_topas_world_lines(f, world_halflengths_sphere_plates)
    write_topas_tsbox(f, containing_box_dims, containing_box_place, containing_box_rot)

    f.write(f'''s:Ge/EnclosingBox/Field                  = "TsElectricFieldMap"
s:Ge/EnclosingBox/ElectricField3DTable   = "{table_path}"
              
''')

    f.write(f'''s:Ge/PlaneUpper/Parent            = "EnclosingBox"
s:Ge/PlaneUpper/Type              = "TsBox"
s:Ge/PlaneUpper/Material          = "G4_Cu" 
d:Ge/PlaneUpper/HLX               = 7.5 mm
d:Ge/PlaneUpper/HLY               = 7.5 mm
d:Ge/PlaneUpper/HLZ               = 0.25 mm
d:Ge/PlaneUpper/TransX            = 0 mm
d:Ge/PlaneUpper/TransY            = 0 mm
d:Ge/PlaneUpper/TransZ            = -7.25 mm
s:Ge/PlaneUpper/DrawingStyle = "solid"

s:Ge/PlaneLower/Parent            = "EnclosingBox"
s:Ge/PlaneLower/Type              = "TsBox"
s:Ge/PlaneLower/Material          = "G4_Cu" 
d:Ge/PlaneLower/HLX               = 7.5 mm
d:Ge/PlaneLower/HLY               = 7.5 mm
d:Ge/PlaneLower/HLZ               = 0.25 mm
d:Ge/PlaneLower/TransX            = 0 mm
d:Ge/PlaneLower/TransY            = 0 mm
d:Ge/PlaneLower/TransZ            = 7.25 mm
s:Ge/PlaneLower/DrawingStyle = "solid"
            
s:Ge/MiddleSphere/Parent            = "EnclosingBox"
s:Ge/MiddleSphere/Type              = "TsSphere"
s:Ge/MiddleSphere/Material          = "G4_Cu" 
d:Ge/MiddleSphere/RMax              = 6 mm
d:Ge/MiddleSphere/TransX            = 0 mm
d:Ge/MiddleSphere/TransY            = 0 mm
d:Ge/MiddleSphere/TransZ            = 0 mm
s:Ge/MiddleSphere/DrawingStyle = "solid"

''')

    return

if __name__ == '__main__':

    ######################################################
    ###               GENERAL PARAMETERS               ###
    ######################################################
    
    # Possible values: 'uniform', 'linear', 'non-uniform'
    MODE = 'non-uniform'       
    output_file_name = r'sphere_plates_testing'



    ######################################################
    ###           NON-UNIFORM MODE PARAMETERS          ###
    ######################################################

    ansys_file_nonuniform = r'ansys_files/sphere_plates_10-3-23.fld'
    ansys_position_scale_nonuniform = 1000  # factor that will give positions in mm when multiplied with ansys position data
    repetition_axis = ''   # '' if no repetition necessary, otherwise 'x', 'y', or 'z'
    repetition_depth = 0   # 0 if no repetition necessary, otherwise units in mm
    






    ######################################################
    ###             UNIFORM TEST PARAMETERS            ###
    ######################################################
    
    box_dimensions_uniform = np.array([5, 10, 3])             # mm
    box_placement_uniform = np.array([4, -4, 10])             # mm
    box_rot_uniform = np.array([0, 0, 0])                     # deg

    efield_vector_uniform = np.array([1, 0, 0])               # kV/mm
    efield_table_points_uniform = np.array([3, 3, 3])         # num of points/axis



    ######################################################
    ###              LINEAR TEST PARAMETERS            ###
    ######################################################

    box_dimensions_linear = np.array([10, 5, 7])              # mm
    box_placement_linear = np.array([4, -4, 5])               # mm
    box_rot_linear = np.array([0, 0, 0])                      # deg

    efield_dir_linear = np.array([1, 1, 0])                   # kV/mm
    efield_grad_linear = 2                                    # kV/mm per mm
    efield_table_points_linear = np.array([11, 11, 11])       # num of points/axis



    ######################################################
    ###           NON-UNIFORM TEST PARAMETERS          ###
    ######################################################

    construct_wire_dynodes_topas_file_nonuniform = False
    construct_sphere_plates_topas_file_nonuniform = False






    ######################################################
    ###                       CODE                     ###
    ######################################################
    if MODE == 'uniform':
        create_topas_files_uniform(output_file_name, box_placement_uniform, box_rot_uniform, box_dimensions_uniform, efield_vector_uniform, efield_table_points_uniform)
    elif MODE == 'linear':
        create_topas_file_linear(output_file_name, box_placement_linear, box_rot_linear, box_dimensions_linear, efield_dir_linear, efield_grad_linear, efield_table_points_linear)
    elif MODE == 'non-uniform':
        create_opera_3d_file_nonuniform(output_file_name, ansys_file_nonuniform, repetition_axis, repetition_depth, ansys_position_scale_nonuniform)
        if construct_wire_dynodes_topas_file_nonuniform:
            depth_wire_dynodes = 64      # mm
            create_topas_file_wire_dynodes_nonuniform(output_file_name, depth_wire_dynodes)
        elif construct_sphere_plates_topas_file_nonuniform:
            create_topas_file_sphere_plates_nonuniform(output_file_name)
