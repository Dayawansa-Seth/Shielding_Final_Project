#Code to ensure that the source used generates 100 mrem worth of neutron radiation at the beam port outlet (50 thermal, 50 fast)

import openmc
import numpy as np
import re
import os

# Constants
#neutrons
source_energy = [0.025, 1e6]  # Energy of neutrons in eV \
source_distribution = [1,1]
source_strength_thermal = 4.3e10  # Source strength for thermal neutrons in neutrons/sec
source_strength_fast = 1.25e2  # Source strength for fast neutrons in neutrons/sec



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
source_thermal = openmc.Source()
source_thermal.space = openmc.stats.Point((0, -10, 100))  # point source of neutrons behind lead plug
source_thermal.energy = openmc.stats.Discrete(source_energy[0], source_distribution[0])  #source
source_thermal.angle = openmc.stats.PolarAzimuthal(mu=openmc.stats.Uniform(0, 0.1), phi=openmc.stats.Uniform(0.0, 2 * np.pi), reference_uvw=(1.0, 0.0, 0.0))
source_thermal.particle = 'neutron'
source_thermal.strength = source_strength_thermal  # Source strength

# Source Definition (Pencil Beam)
source_fast = openmc.Source()
source_fast.space = openmc.stats.Point((0, -10, 100))  # point source of neutrons behind lead plug
source_fast.energy = openmc.stats.Discrete(source_energy[1], source_distribution[1])  #source
source_fast.angle = openmc.stats.PolarAzimuthal(mu=openmc.stats.Uniform(0, 0.1), phi=openmc.stats.Uniform(0.0, 2 * np.pi), reference_uvw=(1.0, 0.0, 0.0))
source_fast.particle = 'neutron'
source_fast.strength = source_strength_fast  # Source strength
# Settings
settings = openmc.Settings()
settings.run_mode = "fixed source"
settings.batches = 10  # Number of batches
settings.particles = 1000000  # Number of particles per batch
settings.source = source_thermal
settings.export_to_xml()

# Tallies
tallies = openmc.Tallies()

# Define energy groups in eV (convert from MeV to eV)
energy_groups = [
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
energy_groups = [group * 1e6 for group in energy_groups]
# Create an EnergyFilter with these bins
energy_filter = openmc.EnergyFilter(energy_groups)

# Flux tally for photons in the detector region with energy filter
point_filter = openmc.CellFilter(air_flux_cell)  # Use the small pseudo-point cell
particle_filter = openmc.ParticleFilter(['neutron'])  # Tally only gamma particles
flux_tally = openmc.Tally(name="Point Flux by Energy")
flux_tally.filters = [point_filter, particle_filter, energy_filter]  # Add filters
flux_tally.scores = ["flux"]
tallies.append(flux_tally)

tallies.export_to_xml()

# Run OpenMC
for x in range(len(source_energy)):
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

    # Dose Response function, Sv-cm^2
    data = [
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


    if x == 0:
        dose_eq = np.zeros(len(data))
        uncertainty = np.zeros(len(data))
    for i, (energy, response) in enumerate(data):
        dose = flux_mean[i] * energy * response * 1e-12 * 1e2 * 1e3 * 3600
        uncert = flux_std_dev[i] * energy * response * 1e-12 * 1e2 * 1e3 * 3600
        if np.isclose(dose, uncert, atol=1e-12):
            dose_eq[i] += 0
            uncertainty[i] += 0
        else:
            dose_eq[i] += dose
            uncertainty[i] += uncert
        
    
    if x != (len(source_energy)-1):
        settings.source = source_fast
        settings.export_to_xml()

print("Ambient Dose Equivalent for neutrons Transmitted by Energy Group:")
for i, (energy, response) in enumerate(data):
    print(f"Group {i+1}: {energy:.3f} MeV: {abs(dose_eq[i]):.3e} ± {uncertainty[i]:.3e}")
print(f"Total dose: {abs(np.sum(dose_eq)):.3e} mrem/hr ± {np.sqrt(np.sum(uncertainty**2)):.3e}")


