
import numpy as np
import random
from scipy.optimize import minimize
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import csv



############################################################################
####################################### Constants ######################################
############################################################################
p0 = 101325 # Pa
T = 288 # K
k = 1.3806e-23 # Boltzmann constant
m = 5.6e-26 # Mass of molecule of air
alt = 6091 # m = 20,000ft
velocity = 27.78 #m/s = 100 km/h
g = 9.81 # kg/m^2
density_air = 1.225 # kg/m^3
#w = 25.0
#l = 100.0
#t = 1
He = 1.114 # kg/m^3 where Helium density = 0.178 kg/m^3 and air density = 1.292 kg/m^3
Re = "1.118e9"
Mach = "0.082"
fiber_glass_strength = 800.0e6
fibre_glass_density = 1500.0
cuboct_scaling_factor = 3.0/2.0
k_const = fiber_glass_strength / fibre_glass_density**cuboct_scaling_factor
scaling_factor_mass = 1.2
l_scaling_factor = 1000.0
density_scaling_factor = 10.0

# proportion:[CD,CDf,CDp]
drag = {10: [1.13444, 0.11138, 1.02306],
15: [1.18496, 0.11262, 1.07233],
20: [1.19983, 0.11319, 1.08664],
25: [1.67896, 0.13436, 1.54460],
30: [1.77221, 0.14106, 1.63116],
35: [1.99681, 0.14574, 1.85107],
40: [2.41140, 0.14942, 2.26198],
45: [2.76347, 0.15309, 2.61038],
50: [3.07368, 0.15632, 2.91737],
55: [3.62485, 0.15818, 3.46667],
60: [3.79110, 0.16237, 3.62873],
65: [4.21302, 0.16458, 4.04844],
70: [4.43665, 0.16834, 4.26831],
75: [5.18289, 0.16682, 5.01607],
80: [5.25104, 0.17005, 5.08099],
85: [4.95898, 0.17758, 4.78140],
90: [4.92524, 0.18233, 4.74292],
95: [4.81300, 0.18779, 4.62521],
100:[6.46638, 0.16600, 6.30039]}

scale_objectives = [10.0e6,10.0e8,10.0e4]
#scale_objectives = [10.0e8,10.0e-2,10.0e0]


def of_1_drag(w_to_l,l):
	CD = 0
	CDf = 1
	CDp = 2
	w = w_to_l * l
	proportion = int(round(w_to_l * 100/5.0)*5.0)
	# check the w:l falls within the range of the drag table
	if proportion < 10 or proportion > 100:
		print "out of proportion " + str(w_to_l)
		return "Error Grace, sorry"
	drag_coefficient = drag[proportion][CDf] 	# lookup drag value from the table
	surface_area = (4 * np.pi * (w/2)**2) + (2 * np.pi * (w/2) * (l-w))	
	drag_force = drag_coefficient * density_air * surface_area * (velocity**2/2)
	return drag_force


def of_2_bending(w_to_l,l,t_to_w):
	w = w_to_l * l
	t = t_to_w * w
	# assuming that the thickness of digital material is filled with air
	volume = np.pi * ((w**2/4)*(l-w)) + ((4/3) * (w/2)**3)
	M = 0.005 * density_air * velocity**2 * volume**2/3 * l
	I = np.pi/4 * ((w/2)**4 - (w/2 - t)**4)
	stress = (M*w)/I
	return stress


def of_3_buoyancy(w_to_l,l,t_to_w,density):
	w = w_to_l * l
	t = t_to_w * w
	mass_cylinder = np.pi * density * (l-w) * ((w/2)**2 - ((w-2*t)/2)**2)
	mass_end_caps = np.pi * density * ((w**3/6) - ((4/3) * (w/2 - t)**3))
	mass_total = (mass_cylinder + mass_end_caps) * scaling_factor_mass
	lifting_volume = np.pi * ((((w-2*t)/2)**2)*(l-w) + ((4/3) * (w/2 - t)**3))
	mass_lifted = lifting_volume * He
	buoyancy = mass_lifted - mass_total
	return -buoyancy


def generate_pareto_data(data_points):
    pareto_variables = []
    for i in range(data_points):
        l_random = 0.02 + (random.random()*0.28)
        w_to_l = 0.1 + (random.random()*0.9)
        t_to_w = 0.1 + (random.random() * 0.5)
        density = 0.0431 + (random.random() * 0.9569)
        variables = [w_to_l,l_random,t_to_w,density]
        pareto_variables.append(variables)

    pareto_objective_functions = []
    for j in range(data_points):
        of_1 = of_1_drag(pareto_variables[j][0], pareto_variables[j][1])/scale_objectives[0]
        of_2 = of_2_bending(pareto_variables[j][0], pareto_variables[j][1], pareto_variables[j][2])/scale_objectives[1]
        of_3 = of_3_buoyancy(pareto_variables[j][0], pareto_variables[j][1], pareto_variables[j][2], pareto_variables[j][3])/scale_objectives[2]	
        if of_3 <= 0:
        	function_evaluations = [of_1,of_2,of_3]
        	pareto_objective_functions.append(function_evaluations)
    print len(pareto_objective_functions)
    	
    return pareto_objective_functions


def plot_biobjective_with_pareto_fronts(J1,J2,J3):
    fig = plt.figure(figsize=(10,10), facecolor="white")   

    plt.subplot(311)
    plt.plot(J1,J2,'c.',  color='magenta')
    #plt.title('Drag vs bending')
    plt.xlabel('Drag')
    plt.ylabel('Bending')
   
    plt.subplot(312)
    plt.plot(J2,J3,'c.',  color='magenta')
    #plt.title('Bending vs Buoyancy')
    plt.xlabel('Bending')
    plt.ylabel('Buoyancy')
    
    plt.subplot(313)
    plt.plot(J1,J3,'c.',  color='magenta')
    #plt.title('Drag vs buoyancy')
    plt.xlabel('Drag')
    plt.ylabel('Buoyancy')
    
    plt.show()

def plot_all_objectives(J1,J2,J3):
    fig = plt.figure(figsize=(10,6), facecolor="white")
    ax = fig.add_subplot(111, projection='3d')
    t1 = "Objective space"
    ax.set_title(t1, fontsize=15)
    ax.scatter(J1,J2,J3)
    ax.set_xlabel('Drag', fontsize = 16)
    ax.set_ylabel('Bending stress',  fontsize = 16)
    ax.set_zlabel('Buoyancy', fontsize = 16)

    plt.show()


#pareto_data = np.array(generate_pareto_data(50000))
#of1 = pareto_data[:,0]
#of2 = pareto_data[:,1]
#of3 = pareto_data[:,2]
#plot_biobjective_with_pareto_fronts(of1,of2,of3)
#plot_all_objectives(of1,of2,of3)


def objective_functions_summed(variables, weighting):
	w_to_l = variables[1]
	t_to_w = variables[2]
	l = variables[0] * l_scaling_factor
	w = w_to_l*l
	t = t_to_w*w
	density = variables[3] * density_scaling_factor

	# Objective function 1
	of1 = of_1_drag(w_to_l,l)/scale_objectives[0]

	# Objective function 2
	of2 = of_2_bending(w_to_l,l,t_to_w)/scale_objectives[1]

	# Objective function 3
	val3 =of_3_buoyancy(w_to_l,l,t_to_w,density)
	if (val3 >-20):
		of3 = 10000
	else:
		of3 = val3/scale_objectives[2]

	#print str(of1) + ", " + str(of2) + ", " + str(of3)
	#print str(of1*weighting[0]) + ", " + str(of2*weighting[1]) + ", " + str(of3*weighting[2])
	sum_ofs = (of1*weighting[0]) + (of2*weighting[1]) + (of3*weighting[2])
	#print sum_ofs
	#print ""
	return (of1*weighting[0]) + (of2*weighting[1]) + (of3*weighting[2])


def objective_functions_separate(variables, weighting):
	w_to_l = variables[1]
	t_to_w = variables[2]
	l = variables[0] * l_scaling_factor
	w = w_to_l*l
	t = t_to_w*w
	density = variables[3] * density_scaling_factor

	# Objective function 1
	of1 = of_1_drag(w_to_l,l)/scale_objectives[0]

	# Objective function 2
	of2 = of_2_bending(w_to_l,l,t_to_w)/scale_objectives[1]

	# Objective function 3
	val3 =of_3_buoyancy(w_to_l,l,t_to_w,density)
	if (val3 >-20):
		of3 = 1000000
	else:
		of3 = val3/scale_objectives[2]
	return of1, of2, of3


def minimise_with_weighting(weighting):
    starting_values = [0.1,0.4,0.01,0.1]
    outputs = []
    for g in range(len(weighting)):
		#print "weighting = " + str(weighting[g])
		result = minimize(objective_functions_summed, starting_values, 
			args=(weighting[g],),
            method='SLSQP',
            bounds=((0.02,0.3),(0.1,1.0),(0.01,0.5),(0.0431,1.0)),
            options={'maxiter':100, 'disp':False, 'eps':0.00001})
		outputs.append(result)
    return outputs
        

weightings = [[1,0.01,0],
           	  [0,1,0],
              [0,0,1],
              [1,0.5,0.5],
              [0.5,1,0.5],
              [0.5,0.5,1.0],
              [0.5,0.5,0.5],
              [0.1,0.6,0.3],
              [0.3,0.1,0.6],
              [0.6,0.3,0.1],
              [0,0.9,0.25],
              [0.25,0.9,0.5]]    

# 0.01 l=119.98420136  w_to_l=0.1  t_to_w=0.01  density=5.0
# 0.1 l=119.983262528  w_to_l=0.1  t_to_w=0.5  density=5.0
# 0 l=119.98420136  w_to_l=0.1  t_to_w=0.01  density=5.0




def generate_weightings(quantity):
	random_weights = []
	for i in range(quantity):
		random_weights.append([random.random(), random.random(), random.random()])
	return random_weights


def minimisation_and_plot(quantity):
	weightings = generate_weightings(quantity) # comment this line out in order to just use the 12 weightings above
	outputs = minimise_with_weighting(weightings)
	OF1 = []
	OF2 = []
	OF3 = []
	base_coordinates = []
	arc1_coordinates = []
	arc2_coordinates = []
	polar_plane = []
	L_coordinates = []
	W_coordinates = []
	values = []
	T = []
	DENSITY = []
	spacing = 0
	of1_counter = 0
	for i in range(len(outputs)):
		of1,of2,of3 = objective_functions_separate(outputs[i].x,weightings[i])
		w_to_l = outputs[i].x[1]
		t_to_w = outputs[i].x[2]
		l = outputs[i].x[0] * l_scaling_factor
		w = w_to_l * l
		t = t_to_w * w
		density = outputs[i].x[3] * density_scaling_factor

		base_coordinates.append("0, " + str(spacing) + ", 0")
		L_coordinates.append(str(l-w) + ", " + str(spacing) + ", 0")
		W_coordinates.append("0, " + str(w+spacing) + ", 0")
		arc1_coordinates.append("0, " + str((w/2)+spacing) + ", 0")
		arc2_coordinates.append(str(l-w) + ", " + str((w/2)+spacing) + ", 0")
		T.append(str(t) + ", 0, 0")
		DENSITY.append(str(density) + ", 0, 0")
		of1,of2,of3 = objective_functions_separate(outputs[i].x,weightings[i])
		OF1.append(of1)
		OF2.append(of2)
		OF3.append(of3)
		values.append([l-w,w,t,density])
		spacing += w+10
		#print ""
		#print "Weighting = " + str(weightings[i])
		#print "Drag = " + str(of1) + " Bending = " + str(of2) + " Buoyancy = " + str(of3)
		#print "l=" + str(l) + "  " + "w_to_l=" + str(w_to_l) + "  " + "t_to_w=" + str(t_to_w) + "  " + "density=" + str(density) 
	
		if of1 < 0.021:
			print "LOW variables = " + str(outputs[i].x)
			of1_counter += 1
		else:
			print "HIG variables = " + str(outputs[i].x)
	print "number of drag samples = " + str(of1_counter	)

	base_coordinates_file = open('airship_generator_data/base_coordinates', 'w')
	for item in base_coordinates:
 		base_coordinates_file.write("%s\n" % item)

	L_coordinates_file = open('airship_generator_data/L_coordinates', 'w')
	for item in L_coordinates:
 		L_coordinates_file.write("%s\n" % item)

	W_coordinates_file = open('airship_generator_data/W_coordinates', 'w')
	for item in W_coordinates:
 		W_coordinates_file.write("%s\n" % item)

	arc1_coordinates_file = open('airship_generator_data/arc1_coordinates', 'w')
	for item in arc1_coordinates:
 		arc1_coordinates_file.write("%s\n" % item)

	arc2_coordinates_file = open('airship_generator_data/arc2_coordinates', 'w')
	for item in arc2_coordinates:
 		arc2_coordinates_file.write("%s\n" % item)

	T_file = open('airship_generator_data/thickneses', 'w')
	for item in T:
 		T_file.write("%s\n" % item)

	density_file = open('airship_generator_data/density', 'w')
	for item in DENSITY:
  		density_file.write("%s\n" % item)


  	#for h in range(len(OF3)):
  	#	print "density:" + str(DENSITY[h]) + " strength:" + str(OF2[h]) 
  	#print len(OF3)
  	plot_biobjective_with_pareto_fronts(OF1,OF2,OF3)
	plot_all_objectives(OF1,OF2,OF3)
	#print DENSITY


minimisation_and_plot(100)

