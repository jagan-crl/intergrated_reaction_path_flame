from irpa_functions import irpa_main_function, freeflame_simulation

########################################################################################################################
# Main script
# Change necessary parameters below
model = 'ffcm1.yaml'
outfile = 'output.txt'      # Change to desired name, saved to same folder where Python scripts are saved
tracing_element = 'C'
fuel_name = 'C2H4'

# Create a flame object from Cantera's FreeFlame simulation - can use functionality provided or reload another solution
flame_obj = freeflame_simulation(fuel_name, model, 1.0, 298, 101325*1.0, 0.21, 0.03)

# Run main reaction path analysis code
irpa_main_function(model, flame_obj, tracing_element, thresh=0.03, outfile=outfile, fuel=fuel_name)
