import openmc
import numpy as np
import re

# Parameters
cube_length = 60 #cm
bp_thickness = 10 #cm
lead_shield_thickness = 10.16 #cm
air_gap = 5 #distance between shield and beam port

# Constants
#neutrons
source_energy = [0.025, 1e6]  # Energy of neutrons in eV \
source_distribution = [1,1]
source_strength_thermal = 2.0e10  # Source strength for thermal neutrons in neutrons/sec
source_strength_fast = 6.0e1  # Source strength for fast neutrons in neutrons/sec

# gammas prompt fission
bins1 = np.linspace(0.1, 0.6, 10, endpoint=False)
stren1 = bins1*6.6*bins1 #photons per fission for those energy bins
bins2 = np.linspace(0.6, 1.5, 10, endpoint=False)
stren2 = bins2*(20.2*np.exp(-1.78*bins2))*bins2 #photons per fission for those energy bins
bins3 = np.linspace(1.5, 10.5, 10, endpoint=False)
stren3 = bins3*(7.2*np.exp(-1.09*bins3))*bins3 #photons per fission for those energy bins
# Combine into one array
beam_energy = np.concatenate((bins1, bins2, bins3))*1e6 #energy bins
beam_strength = np.concatenate((stren1, stren2, stren3)) #photons per fission for those energy bins
max_strength = np.max(beam_strength)
source_strength = 1e3  # Source strength in photons/sec

# Define Materials
concrete = openmc.Material(name="Concrete")
concrete.add_element("O", 0.52)
concrete.add_element("Si", 0.325)
concrete.add_element("Ca", 0.06)
concrete.add_element("Al", 0.035)
concrete.add_element("Fe", 0.06)
concrete.set_density("g/cm3", 2.3)

polyethylene = openmc.Material(name="Polyethylene")
polyethylene.add_element("C", 1.0)
polyethylene.add_element("H", 2.0)
polyethylene.set_density("g/cm3", 0.95)

lead = openmc.Material(name="Lead")
lead.add_element("Pb", 1.0)
lead.set_density("g/cm3", 11.34)

air = openmc.Material(name="Air")
air.add_element("N", 0.784)
air.add_element("O", 0.21)
air.add_element("Ar", 0.006)
air.set_density("g/cm3", 0.001225)

borated_polyethylene = openmc.Material(name="Borated Polyethylene")
borated_polyethylene.add_element("C", 1.0)
borated_polyethylene.add_element("H", 2.0)
borated_polyethylene.add_element("B", 0.1)
borated_polyethylene.set_density("g/cm3", 1.2)

steel = openmc.Material(name="Steel") #not used
steel.add_element("Fe", 0.98)
steel.add_element("C", 0.02)
steel.set_density("g/cm3", 7.87)

# Materials
person = openmc.Material(name="person") #using water to approximate a person
person.add_nuclide("H1", 2.0)
person.add_nuclide("O16", 1.0)
person.set_density("g/cm3", 1.0)

materials = openmc.Materials([concrete, polyethylene, lead, borated_polyethylene, steel, air, person])
materials.export_to_xml()


# Define Geometry
# Surfaces
concrete_left = openmc.XPlane(x0=-100.0, boundary_type="vacuum")
concrete_right = openmc.XPlane(x0=100.0, boundary_type="vacuum")
concrete_front = openmc.YPlane(y0=-200.0, boundary_type="vacuum")
concrete_back = openmc.YPlane(y0=0.0)
concrete_bottom = openmc.ZPlane(z0=0.0)
concrete_top = openmc.ZPlane(z0=200.0, boundary_type="vacuum")
concrete_floor_bottom = openmc.ZPlane(z0=-200.0, boundary_type="vacuum")

# Hole along the y-axis, centered at z=1 m)
radius = 3*2.54 #6 inch diameter
hole_outer_radius = openmc.YCylinder(r=radius, x0=0.0, z0=100)
hole_front = concrete_front
hole_back = concrete_back

poly_front = openmc.YPlane(y0=-170.0)
poly_back = openmc.YPlane(y0=-100.0)
# Polyethylene Plug inside the hole
poly_region = -hole_outer_radius & +poly_front & -poly_back

# lead shutter approximation
lead_front = openmc.YPlane(y0=-20.0)
lead_back = openmc.YPlane(y0=-10.0)
lead_region = -hole_outer_radius & +lead_front & -lead_back


# Concrete Cube Blocking the Hole
cube_left = openmc.XPlane(x0=-cube_length/2)
cube_right = openmc.XPlane(x0=cube_length/2)
cube_front = openmc.YPlane(y0=air_gap)
cube_back = openmc.YPlane(y0=air_gap+cube_length)
cube_bottom = openmc.ZPlane(z0=100-cube_length/2)
cube_top = openmc.ZPlane(z0=100+cube_length/2)
cube_region = +cube_left & -cube_right & +cube_front & -cube_back & +cube_bottom & -cube_top

# Borated Polyethylene around the Concrete Cube
bp_left = openmc.XPlane(x0=-(cube_length/2+bp_thickness/2))
bp_right = openmc.XPlane(x0=(cube_length/2+bp_thickness/2))
#bp_front = openmc.YPlane(y0=5.0)
bp_back = openmc.YPlane(y0=air_gap+(cube_length+bp_thickness))
bp_bottom = openmc.ZPlane(z0=100-(cube_length/2+bp_thickness/2))
bp_top = openmc.ZPlane(z0=100+(cube_length/2+bp_thickness/2))
bp_region = +bp_left & -bp_right & +cube_front & -bp_back & +bp_bottom & -bp_top


# lead_shield Layer outside the Borated Polyethylene
lead_shield_left = openmc.XPlane(x0=-(cube_length/2+bp_thickness/2+lead_shield_thickness/2))
lead_shield_right = openmc.XPlane(x0=(cube_length/2+bp_thickness/2+lead_shield_thickness/2))
#lead_shield_front = openmc.YPlane(y0=5.0)
lead_shield_back = openmc.YPlane(y0=air_gap+(cube_length+bp_thickness+lead_shield_thickness))
lead_shield_bottom = openmc.ZPlane(z0=100.0 - (cube_length/2+bp_thickness/2+lead_shield_thickness/2))
lead_shield_top = openmc.ZPlane(z0=100 + (cube_length/2+bp_thickness/2+lead_shield_thickness/2))
lead_shield_region = +lead_shield_left & -lead_shield_right & +cube_front & -lead_shield_back & +lead_shield_bottom & -lead_shield_top & ~bp_region & ~cube_region


person_height = 5.7*12*2.54 #average height of 5'7"
person_radius = 12*2.54 # radius of 1 ft
#phantom_cyl_1 = openmc.ZCylinder(r=3, x0=0, y0=2)
phantom_cyl_1 = openmc.ZCylinder(r=person_radius, x0=(cube_length/2+bp_thickness+lead_shield_thickness+person_radius+1), y0=person_radius+1)
person_top = openmc.ZPlane(z0=person_height)
phantom_cyl_2 = openmc.ZCylinder(r=person_radius, x0=0.0, y0=(air_gap+(cube_length+bp_thickness+lead_shield_thickness)+person_radius+1))
phantom_sphere_1 = openmc.Sphere(x0=0,y0=-5,z0=100,r=5)
room_back = openmc.YPlane(y0=air_gap+(cube_length+bp_thickness+lead_shield_thickness)+ 2*person_radius + 20, boundary_type="vacuum")



# Define regions
concrete_region = +concrete_left & -concrete_right & +concrete_front & -concrete_back & +concrete_bottom & -concrete_top & +hole_outer_radius
floor_region = +concrete_floor_bottom & -concrete_bottom & +concrete_left & -concrete_right & +concrete_front & -room_back
phantom_1 = -phantom_cyl_1 & +concrete_bottom & -person_top
#phantom_1 = -phantom_sphere_1
phantom_2 = -phantom_cyl_2 & +concrete_bottom & -person_top
air_region = ~phantom_1 & ~phantom_2 & ~lead_region & ~floor_region & ~concrete_region & ~cube_region & ~lead_shield_region & ~bp_region 
#air_region = (-hole_outer_radius & ((-lead_front & +poly_back) | (-poly_front & +concrete_front) | (-concrete_back & +lead_back))) | (+concrete_left & -concrete_right & +concrete_back & -room_back & +concrete_bottom & -concrete_top & ~lead_shield_region & ~bp_region & ~cube_region & ~phantom_1 & ~phantom_2)




# Cells
concrete_cell = openmc.Cell(name="Concrete", region=concrete_region, fill=concrete)
poly_cell = openmc.Cell(name="Polyethylene Plug", region=poly_region, fill=polyethylene)
lead_cell = openmc.Cell(name="Lead Block", region=lead_region, fill=lead)
floor_cell = openmc.Cell(name="Concrete Floor", region=floor_region, fill=concrete)
air_cell = openmc.Cell(name="Air Filling Rest of universe", region=air_region, fill=air)
cube_cell = openmc.Cell(name="Concrete Cube", region=cube_region, fill=air)
bp_cell = openmc.Cell(name="Borated Polyethylene", region=bp_region, fill=borated_polyethylene)
lead_shield_cell = openmc.Cell(name="lead_shield Layer", region=lead_shield_region, fill=lead)
phantom_1_cell = openmc.Cell(name="Person standing to the side of shield", region=phantom_1, fill=person)
phantom_2_cell = openmc.Cell(name="Person standing behind shield", region=phantom_2, fill=person)

# Universe and Geometry
universe = openmc.Universe(cells=[concrete_cell, poly_cell, lead_cell, floor_cell, air_cell, cube_cell, bp_cell, lead_shield_cell, phantom_1_cell, phantom_2_cell])
geometry = openmc.Geometry(universe)
geometry.export_to_xml()


plot = openmc.Plot()
plot.filename = "geometry_plot"
plot.basis = "xy"  # Choose 'xy', 'xz', or 'yz' for the slicing plane
plot.origin = (0, 0, 100)  # Center of the slice (x, y, z in cm)
plot.width = (400, 400)  # Width and height of the plot in cm
plot.pixels = (800, 800)  # Resolution of the plot in pixels
plot.color_by = "material"  # Options: 'material' or 'cell'

# Add the plot to a collection and export
plots = openmc.Plots([plot])
plots.export_to_xml()

# Generate the plot (requires OpenMC to be installed)
openmc.plot_geometry()

# thermal neutron source
source_thermal = openmc.Source()
source_thermal.space = openmc.stats.Point((0, -10, 100))  # cone source of neutrons in fornt of lead plug
source_thermal.energy = openmc.stats.Discrete(source_energy[0], source_distribution[0])  #source
source_thermal.angle = openmc.stats.PolarAzimuthal(mu=openmc.stats.Uniform(0, 0.1), phi=openmc.stats.Uniform(0.0, 2 * np.pi), reference_uvw=(1.0, 0.0, 0.0))
source_thermal.particle = 'neutron'
source_thermal.strength = source_strength_thermal  # Source strength

# fast neutron source
source_fast = openmc.Source()
source_fast.space = openmc.stats.Point((0, -10, 100))  # cone source of neutrons in fornt of lead plug
source_fast.energy = openmc.stats.Discrete(source_energy[1], source_distribution[1])  #source
source_fast.angle = openmc.stats.PolarAzimuthal(mu=openmc.stats.Uniform(0, 0.1), phi=openmc.stats.Uniform(0.0, 2 * np.pi), reference_uvw=(1.0, 0.0, 0.0))
source_fast.particle = 'neutron'
source_fast.strength = source_strength_fast  # Source strength

#gamma source
source_gamma = openmc.Source()
source_gamma.space = openmc.stats.Point((0, -10, 100))  # point source of gammas behind lead plug
source_gamma.energy = openmc.stats.Discrete(beam_energy, beam_strength/max_strength) # prompt fission source
source_gamma.angle = openmc.stats.PolarAzimuthal(mu=openmc.stats.Uniform(0, 0.1), phi=openmc.stats.Uniform(0.0, 2 * np.pi), reference_uvw=(1.0, 0.0, 0.0))
source_gamma.particle = 'photon'
source_gamma.strength = source_strength*max_strength
sources = [source_thermal, source_fast, source_gamma]

# Settings
settings = openmc.Settings()
settings.run_mode = "fixed source"
settings.batches = 10  # Number of batches
settings.particles = 1000000  # Number of particles per batch
settings.max_lost_particles = 100000  # Increase this to a high value if necessary



# Define energy groups in MeV
neutron_energy_groups = [
    0,       # Group 1: Minumum energy in MeV
    2.5e-8,  # Group 1: Maximum energy in MeV
    1.0e-7,  # Group 2
    1.0e-6,  # Group 3
    1.0e-5,  # Group 4
    1.0e-4,  # Group 5
    1.0e-3,  # Group 6
    1.0e-2,  # Group 7
    2.0e-2,  # Group 8
    5.0e-2,  # Group 9
    1.0e-1,  # Group 10
    2.0e-1,  # Group 11
    5.0e-1,  # Group 12
    1.0,     # Group 13
]

# Convert MeV to eV
neutron_energy_groups = [group * 1e6 for group in neutron_energy_groups]
# Create an EnergyFilter with these bins
neutron_energy_filter = openmc.EnergyFilter(neutron_energy_groups)
particle_filter = openmc.ParticleFilter(['neutron'])  # Tally only gamma particles

# Flux tally for neutrons in the detector region with energy filter
neutron_point_filter_phantom_1 = openmc.CellFilter(phantom_1_cell)  # Phantom 1 (next to shield)
neutron_flux_tally_phantom_1 = openmc.Tally(name="Point Flux by Energy Phantom 1")
neutron_flux_tally_phantom_1.filters = [neutron_point_filter_phantom_1, particle_filter, neutron_energy_filter]  # Add filters
neutron_flux_tally_phantom_1.scores = ["flux"]


neutron_point_filter_phantom_2 = openmc.CellFilter(phantom_2_cell)  # phantom 2 (behind shield)
neutron_flux_tally_phantom_2 = openmc.Tally(name="Point Flux by Energy Phantom 2")
neutron_flux_tally_phantom_2.filters = [neutron_point_filter_phantom_2, particle_filter, neutron_energy_filter]  # Add filters
neutron_flux_tally_phantom_2.scores = ["flux"]


# Define gamma energy groups in MeV
gamma_energy_groups = [
    0,
    1.0e-2,
    1.50e-2,
    2.0e-2,
    3.0e-2,
    4.0e-2,
    5.0e-2,
    6.00e-2,
    8.00e-2,
    1e-1,
    1.5e-1,
    2.0e-1,
    0.3,
    0.4,
    0.5,
    0.6,
    0.8,
    1.0,
    2.0,
    4.0,
    6.0,
    8.0,
    10.0,
]

# Convert MeV to eV
gamma_energy_groups = [group * 1e6 for group in gamma_energy_groups]
# Flux tally for photons in the detector region with energy filter
particle_filter = openmc.ParticleFilter(['photon'])  # Tally only gamma particles
gamma_energy_filter = openmc.EnergyFilter(gamma_energy_groups)

gamma_point_filter_phantom_1 = openmc.CellFilter(phantom_1_cell)  # Phantom 1
gamma_flux_tally_phantom_1 = openmc.Tally(name="Point Flux by Energy Phantom 1")
gamma_flux_tally_phantom_1.filters = [gamma_point_filter_phantom_1, particle_filter, gamma_energy_filter]  # Add filters
gamma_flux_tally_phantom_1.scores = ["flux"]


gamma_point_filter_phantom_2 = openmc.CellFilter(phantom_2_cell)  # Phantom 2
gamma_flux_tally_phantom_2 = openmc.Tally(name="Point Flux by Energy Phantom 2")
gamma_flux_tally_phantom_2.filters = [gamma_point_filter_phantom_2, particle_filter, gamma_energy_filter]  # Add filters
gamma_flux_tally_phantom_2.scores = ["flux"]


tally_types = [[neutron_flux_tally_phantom_1, neutron_flux_tally_phantom_2], [neutron_flux_tally_phantom_1, neutron_flux_tally_phantom_2], [gamma_flux_tally_phantom_1, gamma_flux_tally_phantom_2]]





# Neutron Dose Response function, Sv-cm^2
neutron_data = [
    (2.5e-8,4.0),  # Group 1: Maximum energy in MeV
    (1.0e-7,4.40),  # Group 2
    (1.0e-6,4.82),  # Group 3
    (1.0e-5,4.46),  # Group 4
    (1.0e-4,4.14),  # Group 5
    (1.0e-3,3.83),  # Group 6
    (1.0e-2,4.53),  # Group 7
    (2.0e-2,5.87),  # Group 8
    (5.0e-2,10.9),  # Group 9
    (1.0e-1,19.8),  # Group 10
    (2.0e-1,38.6),  # Group 11
    (5.0e-1,87.0),  # Group 12
    (1.0,143),     # Group 13
]

gamma_data = [
    (1.0e-2, 0.0496),  # Energy (MeV), Dose Response Function (Sv·cm²)
    (1.50e-2, 0.129),
    (2.0e-2, 0.211),
    (3.0e-2, 0.307),
    (4.0e-2, 0.345),
    (5.0e-2, 0.363),
    (6.00e-2, 0.382),
    (8.00e-2, 0.441),
    (1e-1, 0.519),
    (1.5e-1, 0.754),
    (2.0e-1, 1.00),
    (0.3,1.51),
    (0.4,2.00),
    (0.5,2.47),
    (0.6,2.91),
    (0.8,3.73),
    (1.0,4.48),
    (2.0,7.45),
    (4.0,11.9),
    (6.0,15.7),
    (8.0,19.3),
    (10.0,23.0),
]

datas = [neutron_data, neutron_data, gamma_data]
total_dose_1 = 0
total_uncertainty_1 = 0
total_dose_2 = 0
total_uncertainty_2 = 0

# Run OpenMC
for x in range(len(datas)):
    settings.source = sources[x]
    settings.export_to_xml()
    # Tallies
    tallies = openmc.Tallies(tally_types[x])
    tallies.export_to_xml()
    openmc.run(threads=110, output=False)

    # Post-Processing

    flux_value_regex = r"Flux\s+([\d\.e\+\-]+)\s\+\/\-\s([\d\.e\+\-]+)"

    # Extract tally results
    flux_mean = []  # Array of mean current values for each energy group
    flux_std_dev = []  # Array of uncertainties for each group

    with open("tallies.out", 'r') as file:
        for line in file:
            flux_match = re.search(flux_value_regex, line)
            if flux_match:
                flux, uncert = map(float, flux_match.groups())
                flux_mean.append(flux)
                flux_std_dev.append(uncert)

    

    data = datas[x]
    dose_eq_1 = np.zeros(len(data))
    uncertainty_1 = np.zeros(len(data))
    dose_eq_2 = np.zeros(len(data))
    uncertainty_2 = np.zeros(len(data))
    for i, (energy, response) in enumerate(data):
        dose_1 = flux_mean[i] * energy * (response * 2)* 1e-12 * 1e2 * 1e3 * 3600 #flux_energy (MeV) * response * 2 (to follow 1985 ICRP) * 1e5 (comvert to mrem/s) * 3600 (mrem/hr)
        uncert_1 = flux_std_dev[i] * energy * (response * 2) * 1e-12 * 1e2 * 1e3 * 3600
        dose_2 = flux_mean[i+len(data)] * energy * (response * 2)* 1e-12 * 1e2 * 1e3 * 3600 #flux_energy (MeV) * response * 2 (to follow 1985 ICRP) * 1e5 (comvert to mrem/s) * 3600 (mrem/hr)
        uncert_2 = flux_std_dev[i+len(data)] * energy * (response * 2) * 1e-12 * 1e2 * 1e3 * 3600
        if np.isclose(dose_1, uncert_1, atol=1e-12):
            dose_eq_1[i] += 0
            uncertainty_1[i] += 0
        else:
            dose_eq_1[i] += dose_1
            uncertainty_1[i] += uncert_1
        
        if np.isclose(dose_2, uncert_2, atol=1e-12):
            dose_eq_2[i] += 0
            uncertainty_2[i] += 0
        else:
            dose_eq_2[i] += dose_2
            uncertainty_2[i] += uncert_2
    
    
    rad_type = 'neutrons'
    if x==0:
        neutron_energy = 'thermal'
    elif x==1:
        neutron_energy = 'fast'
    if x == 2:
        neutron_energy = ''
        dose_eq_1 /= 2
        uncertainty_1 /= 2
        dose_eq_2 /= 2
        uncertainty_2 /= 2
        rad_type = 'gammas'
    
    print(f'Ambient Dose Equivalent for {neutron_energy} {rad_type} Transmitted by Energy Group Phantom 1:')
    for i, (energy, response) in enumerate(data):
        print(f"Group {i+1}: {energy:.3f} MeV: {abs(dose_eq_1[i]):.3e} ± {uncertainty_1[i]:.3e}")
    print(f"Total dose: {abs(np.sum(dose_eq_1)):.3e} mrem/hr ± {np.sqrt(np.sum(uncertainty_1**2)):.3e} \n")
    
    print(f'Ambient Dose Equivalent for {neutron_energy} {rad_type} Transmitted by Energy Group Phantom 2:')
    for i, (energy, response) in enumerate(data):
        print(f"Group {i+1}: {energy:.3f} MeV: {abs(dose_eq_2[i]):.3e} ± {uncertainty_2[i]:.3e}")
    print(f"Total dose: {abs(np.sum(dose_eq_2)):.3e} mrem/hr ± {np.sqrt(np.sum(uncertainty_2**2)):.3e}")
    
    total_dose_1 += np.sum(dose_eq_1)
    total_uncertainty_1 += np.sum(uncertainty_1**2)
    
    total_dose_2 += np.sum(dose_eq_2)
    total_uncertainty_2 += np.sum(uncertainty_2**2)

print(f"Total dose from all sources for phantom 1: {(total_dose_1):.3e} mrem/hr ± {np.sqrt(total_uncertainty_1):.3e} \n")
print(f"Total dose from all sources for phantom 2: {(total_dose_2):.3e} mrem/hr ± {np.sqrt(total_uncertainty_2):.3e}")


