#Code to ensure that the source used generates 100 mrem worth of gamma radiation at the beam port outlet

import openmc
import numpy as np
import re

# Constants
#beam_energy = 200e3  # Photon energy in eV (200 keV)
# Generate the bins
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

materials = openmc.Materials([concrete, polyethylene, lead, air])
materials.export_to_xml()

# Define Geometry
# Surfaces
concrete_left = openmc.XPlane(x0=-100.0, boundary_type="vacuum")
concrete_right = openmc.XPlane(x0=100.0, boundary_type="vacuum")
concrete_front = openmc.YPlane(y0=-200.0, boundary_type="vacuum")
concrete_back = openmc.YPlane(y0=0.0, boundary_type="vacuum")
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



# Define regions
concrete_region = +concrete_left & -concrete_right & +concrete_front & -concrete_back & +concrete_bottom & -concrete_top & +hole_outer_radius
floor_region = +concrete_floor_bottom & -concrete_bottom & +concrete_left & -concrete_right & +concrete_front & -concrete_back
air_region = -hole_outer_radius & ((-lead_front & +poly_back) | (-poly_front & +concrete_front) )
air_flux_region = -hole_outer_radius & -concrete_back & +lead_back

# Cells
concrete_cell = openmc.Cell(name="Concrete", region=concrete_region, fill=concrete)
poly_cell = openmc.Cell(name="Polyethylene Plug", region=poly_region, fill=polyethylene)
lead_cell = openmc.Cell(name="Lead Block", region=lead_region, fill=lead)
floor_cell = openmc.Cell(name="Concrete Floor", region=floor_region, fill=concrete)
air_cell = openmc.Cell(name="Air Filling Rest of Hole", region=air_region, fill=air)
air_flux_cell = openmc.Cell(name="air at exit of beam port", region=air_flux_region, fill=air)

# Universe and Geometry
universe = openmc.Universe(cells=[concrete_cell, poly_cell, lead_cell, floor_cell, air_cell, air_flux_cell])
geometry = openmc.Geometry(universe)
geometry.export_to_xml()



# Source Definition (Pencil Beam)
source = openmc.Source()
source.space = openmc.stats.Point((0, -10, 100))  # point source of neutrons behind lead plug
source.energy = openmc.stats.Discrete(beam_energy, beam_strength/max_strength) # prompt fission source
source.angle = openmc.stats.PolarAzimuthal(mu=openmc.stats.Uniform(0, 0.1), phi=openmc.stats.Uniform(0.0, 2 * np.pi), reference_uvw=(1.0, 0.0, 0.0))
source.particle = 'photon'
source.strength = source_strength*max_strength



# Settings
settings = openmc.Settings()
settings.run_mode = "fixed source"
settings.batches = 100
settings.particles = 100000
settings.max_lost_particles = 100000  
settings.source = source
settings.export_to_xml()

# Tallies
tallies = openmc.Tallies()

# Define gamma energy groups in MeV
gamma_energy_groups = [
    0,       # Group 1: Minumum energy in MeV
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
energy_bins = [group * 1e6 for group in gamma_energy_groups]
# Create an EnergyFilter with these bins
energy_filter = openmc.EnergyFilter(energy_bins)

# Flux tally for photons in the detector region with energy filter
point_filter = openmc.CellFilter(air_flux_cell)  # Use the small pseudo-point cell
particle_filter = openmc.ParticleFilter(['photon'])  # Tally only gamma particles
energy_filter = openmc.EnergyFilter(energy_bins)

flux_tally = openmc.Tally(name="Point Flux by Energy")
flux_tally.filters = [point_filter, particle_filter, energy_filter]  # Add filters
flux_tally.scores = ["flux"]
tallies.append(flux_tally)


# Export the tally to the XML file
tallies.export_to_xml()

# Run the simulation
openmc.run(threads=110)

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

print(flux_mean)

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

print("Ambient Dose Equivalent for neutrons Transmitted by Energy Group:")
dose_eq = np.zeros(len(gamma_data))
uncertainty = np.zeros(len(gamma_data))
for i, (energy, response) in enumerate(gamma_data):
    dose_eq[i] = flux_mean[i] * energy * response * 1e-12 * 1e2 * 1e3 * 3600
    uncertainty[i] = flux_std_dev[i] * energy * response * 1e-12 * 1e2 * 1e3 * 3600
    print(f"Group {i+1}: {energy:.3f} MeV: {abs(dose_eq[i]):.3e} ± {uncertainty[i]:.3e}")

print(f"Total dose: {abs(np.sum(dose_eq)):.3e} mrem/hr ± {np.sqrt(np.sum(uncertainty**2)):.3e}")
