import cantera as ct
import numpy as np
import math
import re


ct.suppress_deprecation_warnings()
ct.suppress_thermo_warnings()


# Define functions to be used in IRPA code

def count_atoms(fuel_name):
    # Use regular expressions to find carbon and hydrogen counts
    carbon_pattern = r'C(\d*)'
    hydrogen_pattern = r'H(\d*)'

    # Find all matches
    carbon_match = re.search(carbon_pattern, fuel_name)
    hydrogen_match = re.search(hydrogen_pattern, fuel_name)

    # Count carbon atoms
    if carbon_match:
        carbon_count = int(carbon_match.group(1)) if carbon_match.group(1) else 1
    else:
        carbon_count = 0

    # Count hydrogen atoms
    if hydrogen_match:
        hydrogen_count = int(hydrogen_match.group(1)) if hydrogen_match.group(1) else 1
    else:
        hydrogen_count = 0

    return carbon_count, hydrogen_count


def calculate_oxygen_moles(fuel):
    n_carbon, n_hydrogen = count_atoms(fuel)
    oxygen_moles = n_carbon + n_hydrogen/4
    return oxygen_moles


def freeflame_simulation(fuel, model, phi, T, P, o2_frac, domain, renormalize=False):
    # model_name = input('Give the model name with _fne if PD model (ex. FFCM, ANL, USCII or FFCM_fne, etc).')
    # Module calculates free flame simulation and returns flame object
    width = domain
    loglevel = 1

    # Use O2 fraction from inputs to calculate N2 moles
    o2_moles = calculate_oxygen_moles(fuel)
    n2_moles = (o2_moles * (1 - o2_frac)) / o2_frac

    # Create gas solution object
    gas = ct.Solution(model)
    # gas.TPX = T, P, {'CH4': 1 * phi, 'O2': 2.0, 'N2': n2_moles}
    gas.TPX = T, P, {fuel: 1 * phi, 'O2': o2_moles, 'N2': n2_moles}

    # Set up flame object
    flame = ct.FreeFlame(gas, width=width)
    flame.set_refine_criteria(ratio=3, slope=0.02, curve=0.02)
    # f.show_solution()

    # Solve with mixture-averaged transport model
    flame.transport_model = 'Mix'
    if renormalize:
        flame.solve(loglevel=loglevel, refine_grid=True, auto=True)
    else:
        flame.solve(loglevel=loglevel, refine_grid=True, auto=False)

    return flame


def molecule_list(gas_obj, element):
    # Create a list of possible species from the gas object
    model_species = gas_obj.species_names
    # Create an empty 2D tuple to populate with species names and corresponding indices
    molecules_list_e = [[] for g in range(2)]

    # Iterate through the list of possible species
    for h in range(len(model_species)):
        temp_str = model_species[h]

        # Look for species that contain tracing element and add them to the list, along with their corresponding index
        if element in temp_str:
            molecules_list_e[0].append(temp_str)
            molecules_list_e[1].append(h)

    return molecules_list_e


def max_flux_calc(flx_max, int_rpd_data_list):
    for i in range(len(int_rpd_data_list)):        # Iterate through the list of fluxes
        temp_flux = int_rpd_data_list[i]           # Create a temporary variable for the flux
        if temp_flux < 0.0:
            temp_flux = -1*temp_flux               # Make negative fluxes temporarily positive for comparison

        if temp_flux > flx_max:                    # Compare temporary flux with max_flux value
            flx_max = temp_flux                    # If temporary flux > max_flux, reassign max_flux value
    # For very small max fluxes, set a baseline
    if flx_max < 1E-10:
        flx_max = 1E-10

    return flx_max


def total_flux_calc(primary_sp, secondary_sp, net_flux, fuel):
    total_flux = 0
    for i in range(len(primary_sp)):
        if primary_sp[i] == fuel:
            # print(net_flux[i])
            total_flux = total_flux + net_flux[i]
        if secondary_sp[i] == fuel:
            total_flux = total_flux - net_flux[i]
    return total_flux


def path_weight_calc(flx_ratio, thresh):
    # Taken from Cantera source code to scale path weights by flux
    path_width = 1 - 4 * math.log10(flx_ratio / thresh) / math.log10(thresh) + 1.0

    return path_width


def arrow_size_calc(path_width):
    # Taken from Cantera source code to scale arrow size by flux in diagram
    arrow_size = 0.5 * path_width

    return arrow_size


def output_file_create(filename, name):
    # Create the output file to give to DOT to generate diagram
    output_file = open(filename, "a")

    # Write first two lines
    output_file.write('digraph ' + name + ' {\nnode[fontsize=28, shape="box"]\n')

    return output_file


def get_diagram_data(flame_obj, index, model, element):
    gas_obj = ct.Solution(model)
    # Set gas object to properties of point in the flame
    gas_obj.TPX = flame_obj.T[index], flame_obj.P, flame_obj.X[:, index]
    diagram = ct.ReactionPathDiagram(gas_obj, element)  # Create Reaction Path Diagram object
    rpd_data = diagram.get_data()

    return rpd_data


def irpa(modelname, flame, tracing_element):
    domain = flame.grid

    # Establish a gas object
    gas = ct.Solution(modelname)

    # Returns a tuple of molecules names and indices containing the desired tracing element
    carbon_molecules = molecule_list(gas, tracing_element)

    # Create empty matrices to fill with data from the reaction path diagram function
    # flux_matrix_forward = []
    # flux_matrix_backward = []
    flux_matrix_net = []

    # Loop through every grid point in the flame domain
    for i in range(len(domain)):
        # if i > 75:
        #     break
        print('Grid Point: ', int(i))
        # flux_matrix_forward.append(np.zeros((len(carbon_molecules[0]), len(carbon_molecules[0]))))
        # flux_matrix_backward.append(np.zeros((len(carbon_molecules[0]), len(carbon_molecules[0]))))
        flux_matrix_net.append(np.zeros((len(carbon_molecules[0]), len(carbon_molecules[0]))))

        # Get the flux data from the reaction path diagram object at a given point 'i' in the flame
        # temp_i = flame.T[i]
        # species_i = flame.X[:,i]
        data = get_diagram_data(flame, i, modelname, tracing_element)

        temp_data = data.split('\n')
        del temp_data[:2]               # Delete the first two lines of data object since they don't contain usable data
        for j in range(len(temp_data)):
            temp_var = temp_data[j].split(' ')
            for k in range(len(carbon_molecules[0])):
                for l in range(len(carbon_molecules[0])):
                    if temp_var[0] == carbon_molecules[0][k]:
                        if temp_var[1] == carbon_molecules[0][l]:
                            # flux_matrix_forward[-1][k, l] = float(temp_var[2])
                            # flux_matrix_backward[-1][k, l] = float(temp_var[3])
                            flux_matrix_net[-1][k, l] = float(temp_var[2]) + float(temp_var[3])

        del data

    print('Data Parsed')
    integration_matrix = [[] for i in range(len(flux_matrix_net))]
    for i in range(len(integration_matrix)):
        integration_matrix[i] = np.zeros((len(carbon_molecules[0]), len(carbon_molecules[0])))

    for i in range(len(carbon_molecules[0])):
        for j in range(len(carbon_molecules[0])):
            temp_var = []
            for k in range(len(integration_matrix)):
                temp_var.append(flux_matrix_net[k][i, j])
            for k in range(len(temp_var)):
                if k > 0:
                    # Perform a trapz integration over all fluxes for each species combination
                    integration_matrix[k][i,j] = np.trapz(temp_var[0:k], domain[0:k])

    # The last matrix in the integration tuple will contain the integration from the entire flame domain
    total_fluxes = integration_matrix[-1]

    # Create a 3D tuple of lists to put species pairs and fluxes into
    integral_rpd_data = [[] for i in range(3)]

    for i in range(len(carbon_molecules[0])):
        sp1 = carbon_molecules[0][i]                   # Create variable for first species of interest
        for j in range(len(carbon_molecules[0])):
            sp2 = carbon_molecules[0][j]               # Create variable for second species of interest
            #
            # Filter out instances where sp1 is the same as sp2
            if sp1 != sp2:
                # Create a net flow variable that corresponds to a given species pair i,j
                net_flow_var = total_fluxes[i,j]

                if net_flow_var != 0.0:

                    integral_rpd_data[0].append(sp1)
                    integral_rpd_data[1].append(sp2)
                    integral_rpd_data[2].append(net_flow_var)

    sp1_list = integral_rpd_data[0]
    sp2_list = integral_rpd_data[1]
    flux_data = integral_rpd_data[2]

    return sp1_list, sp2_list, flux_data


def irpa_data_process_write(sp1_list, sp2_list, flux_data, modelname, threshold, outfile_name, fuelname, max_flux=False):
    ##################################################################################################################
    # Write to a .txt file readable by DOT
    # Start by outlining default variables - taken from Cantera source code - ReactionPath.cpp
    name = 'reaction_paths'
    m_font = "Helvetica"              # Reaction path font
    fontsize = 24
    hue = 0.7
    bright = 0.9

    gas2 = ct.Solution(modelname)

    # Scaling/normalization parameter can either be max_flux or total_flux --> this is user-specified
    if max_flux:
        # Find the max flux
        scaling_parameter = max_flux_calc(flx_max=0.0, int_rpd_data_list=flux_data)
        print('Max flux is:', scaling_parameter)
    else:
        scaling_parameter = total_flux_calc(sp1_list, sp2_list, flux_data, fuelname)    #, fuel=fuel_name
        print('Total carbon flux', scaling_parameter)


    begin_node = [[] for i in range(2)]
    end_node = [[] for i in range(2)]

    flux_percent_list = []
    path_weight = []

    # Set beginning and end of path based on sign of the flux
    # def node_sorting(node_begin_list, node_end_list):
    for i in range(len(flux_data)):
        sp1 = sp1_list[i]
        sp2 = sp2_list[i]
        flx = flux_data[i]
        if flx > 0.0:          # If flux is positive
            flux_ratio = (flx / scaling_parameter)
            if flux_ratio >= threshold:         # If flux ratio exceeds threshold
                begin_node[0].append(sp1)
                begin_node[1].append(gas2.species_index(sp1))
                end_node[0].append(sp2)
                end_node[1].append(gas2.species_index(sp2))

                flux_percent_list.append(flux_ratio*100)
                path_weight.append(path_weight_calc(flux_ratio, threshold))

        else:                  # If flux is negative, need to switch the start and end node
            flux_ratio = -flx / scaling_parameter
            if flux_ratio >= threshold:         # If flux ratio exceeds threshold
                begin_node[0].append(sp2)
                begin_node[1].append(gas2.species_index(sp2))
                end_node[0].append(sp1)
                end_node[1].append(gas2.species_index(sp1))

                flux_percent_list.append(flux_ratio*100)
                path_weight.append(path_weight_calc(flux_ratio, threshold))

    output_file = output_file_create(outfile_name, name)

    # def write_paths_to_file(output_file, l_width):
    # Start by writing the paths to the file
    for r in range(len(path_weight)):
        start_node = 's' + str(begin_node[1][r])
        last_node = 's' + str(end_node[1][r])

        arrow_width = arrow_size_calc(path_weight[r])        # Calculate arrow size from the path weight
        #
        output_file.write('{0} -> {1}'.format(start_node, last_node))
        output_file.write('[fontname="{0}", penwidth={1}, arrowsize={2}'.format(m_font, round(path_weight[r], 2),
                                                                                round(arrow_width, 2)))
        output_file.write(', color="{0}, {1}, {2}", label=" {3}", fontsize="{4}"];\n'.format(hue, round(flux_percent_list[r] + 1, 2), bright,
                                                                             (round(flux_percent_list[r])), fontsize))

    used_nodes = []
    # Then write the list of nodes to the file
    # Check for species in begin and end nodes lists
    for m in range(len(begin_node[0])):
        sp_name = begin_node[0][m]
        sp_node = 's' + str(begin_node[1][m])
        if sp_name in used_nodes:
            None
            # print('Species {0} already written to DOT file'.format(sp_name))

        # If species not yet added to file, create a node for species in the outfile
        else:
            output_file.write('{0} [ fontname="{1}", label="{2}"];\n'.format(sp_node, m_font, sp_name))
            used_nodes.append(sp_name)
            # print('Species {0} added to DOT file'.format(sp_name))

    for m in range(len(end_node[0])):
        sp_name = end_node[0][m]
        sp_node = 's' + str(end_node[1][m])
        if sp_name in used_nodes:
            None
            # print('Species {0} already written to DOT file'.format(sp_name))

        # If species not yet added to file, create a node for species in the outfile
        else:
            output_file.write('{0} [ fontname="{1}", label="{2}"];\n'.format(sp_node, m_font, sp_name))
            used_nodes.append(sp_name)
            # print('Species {0} added to DOT file'.format(sp_name))

    output_file.write('}')
    output_file.close()

    print('Integrated reaction path analysis complete.')

    return None


def irpa_main_function(model, flame, tracing_element, thresh, outfile, fuel):
    # Calls two main functions to run IRPA
    # This function is called separately for fluidity of use
    sp1, sp2, flux = irpa(model, flame, tracing_element)
    irpa_data_process_write(sp1, sp2, flux, modelname=model, threshold=thresh, outfile_name=outfile, fuelname=fuel,
                            max_flux=False)
    return None
