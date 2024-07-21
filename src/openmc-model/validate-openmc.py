import openmc
import numpy as np
import os
import matplotlib.pyplot as plt
from numpy.polynomial.legendre import Legendre

class ControlRod:
    def __init__(self, height, radius, num_sections, coefficients):
        self.height = height
        self.radius = radius
        self.num_sections = num_sections
        self.section_height = height / (num_sections + 1)
        self.coefficients = coefficients
        self.radii = self._generate_radii()

    def _generate_radii(self):
        x = np.linspace(0, self.height, self.num_sections + 1)
        xi = 2 * x / self.height - 1  # map to Legendre polynomial interval [-1, 1]
        scaling_factor = Legendre(self.coefficients)(xi) * 0.1 + 1
        radii = self.radius * scaling_factor

        # Scale radii to maintain the volume
        V_target = np.pi * self.height * self.radius ** 2
        V_radii = np.sum(np.pi * self.section_height * radii ** 2)
        radii *= (V_target / V_radii) ** (1 / 2)
        
        return radii

    def get_sections(self, fuel_height, translation_height):
        sections = []
        num = 1
        for i, radius_section in enumerate(self.radii):
            if self.section_height * i <= translation_height:
                z_top = fuel_height / 2 - translation_height + self.section_height * (num)
                z_bottom = fuel_height / 2 - translation_height + self.section_height * (num-1)
                section_surface = openmc.ZCylinder(r=self.radii[-num], boundary_type='transmission')
                sections.append((section_surface, z_top, z_bottom))
                num += 1
            else:
                break
        return sections

def create_openmc_model(translation_height, coefficients):
    # Path to cross-section data
    notebook_path = os.getcwd()
    path = notebook_path + '/cross/cross_sections.xml'
    openmc.config['cross_sections'] = path

    # Define materials
    control_rod_material = openmc.Material(name = 'control_rod_material')
    control_rod_material.add_element('B', 1.0)
    control_rod_material.set_density('g/cm3', 20.52)

    moderator = openmc.Material(name="moderator")
    moderator.set_density('g/cm3', 1.00)
    moderator.add_nuclide('H1', 2.0)
    moderator.add_nuclide('O16', 1.0)
    moderator.add_s_alpha_beta('c_H_in_H2O')

    fuel = openmc.Material(name="fuel")
    fuel.set_density('g/cm3', 10.0)
    fuel.add_nuclide('U235', 0.02)
    fuel.add_nuclide('U238', 0.2)
    fuel.add_nuclide('O16', 0.78)

    cladding = openmc.Material(name="clad")
    cladding.set_density('g/cm3', 6.56)
    cladding.add_element('Zr', 1.0)

    # Create Materials collection
    mats = openmc.Materials([control_rod_material, moderator, fuel, cladding])
    mats.export_to_xml()

    # Define the cylindrical geometry
    control_rod_radius = 5  # in cm
    fuel_inner_radius = control_rod_radius + 3 # in cm
    fuel_outer_radius = fuel_inner_radius + 20  # in cm
    cladding_outer_radius = fuel_outer_radius +2  # in cm
    moderator_outer_radius = cladding_outer_radius + 5
    fuel_height = 50  # in cm

    # Surfaces
    control_rod_cylinder = openmc.ZCylinder(r=control_rod_radius, boundary_type = 'transmission')
    fuel_inner_cylinder = openmc.ZCylinder(r=fuel_inner_radius, boundary_type = 'transmission')
    fuel_outer_cylinder = openmc.ZCylinder(r=fuel_outer_radius, boundary_type = 'transmission')
    cladding_outer_cylinder = openmc.ZCylinder(r=cladding_outer_radius, boundary_type = 'transmission')
    moderator_outer_cylinder = openmc.ZCylinder(r=moderator_outer_radius, boundary_type='reflective')
    top = openmc.ZPlane(z0=fuel_height / 2, boundary_type='reflective')
    bottom = openmc.ZPlane(z0=-fuel_height / 2, boundary_type='reflective')

    # Parameters for the control rod
    control_rod_height = 50  # in cm
    num_sections = 399

    control_rod = ControlRod(control_rod_height, control_rod_radius, num_sections,
                            coefficients)

    control_rod_sections = control_rod.get_sections(control_rod_height, translation_height)

    # Create cells for each section of the control rod
    control_rod_cells = []

    moderator_region = -fuel_inner_cylinder & -top & +bottom
    for section_surface, z_top, z_bottom in control_rod_sections:
        top_plane = openmc.ZPlane(z0=z_top, boundary_type = 'transmission') 
        bottom_plane = openmc.ZPlane(z0=z_bottom, boundary_type = 'transmission')
        control_rod_cell = openmc.Cell()
        control_rod_cell.region = -section_surface & +bottom_plane & -top_plane
        control_rod_cell.fill = control_rod_material
        control_rod_cells.append(control_rod_cell)
        moderator_region & +section_surface
        moderator_region & -top_plane
        moderator_region & +bottom_plane

    # Moderator cell (central region)
    moderator_cell = openmc.Cell(fill=moderator, region=moderator_region)

    # Fuel cell (hollow cylinder region)
    fuel_region = +fuel_inner_cylinder & -fuel_outer_cylinder & -top & +bottom
    fuel_cell = openmc.Cell(fill=fuel, region=fuel_region)

    # Cladding cell
    cladding_region = +fuel_outer_cylinder & -cladding_outer_cylinder & -top & +bottom
    cladding_cell = openmc.Cell(fill=cladding, region=cladding_region)

    # Outer moderator cell (outside cladding)
    outer_moderator_region = +cladding_outer_cylinder & -top & +bottom & -moderator_outer_cylinder
    outer_moderator_cell = openmc.Cell(fill=moderator, region=outer_moderator_region)

    # Create a root universe
    root = openmc.Universe(cells=[moderator_cell, fuel_cell, cladding_cell, outer_moderator_cell] + control_rod_cells)

    # Create geometry and export
    geometry = openmc.Geometry(root)
    geometry.export_to_xml()

    # Create mesh tallies for axial power
    # Define cylindrical mesh parameters
    r_grid = np.linspace(0.0, moderator_outer_radius, 20)  # Radial grid from 0 to cladding_outer_radius
    z_grid = np.linspace(-fuel_height / 2, fuel_height / 2, 30)  # Axial grid from -fuel_height/2 to fuel_height/2
    phi_grid = np.linspace(0.0, 2 * np.pi, 60)  # Azimuthal grid covering 0 to 2*pi
    origin = (0.0, 0.0, 0.0)  # Origin of the cylindrical mesh

    # Create cylindrical mesh
    cylindrical_mesh = openmc.CylindricalMesh(r_grid=r_grid, z_grid=z_grid, phi_grid=phi_grid, origin=origin)

    # Define mesh tally filter using the cylindrical mesh
    mesh_filter = openmc.MeshFilter(mesh=cylindrical_mesh)

    # Define the tally
    t = openmc.Tally(name='axial peaking')
    t.filters = [mesh_filter]
    t.scores = ['flux']
    tallies = openmc.Tallies([t])
    tallies.export_to_xml()

    # Create source
    source = openmc.IndependentSource()
    source.space = openmc.stats.CylindricalIndependent(r=openmc.stats.Uniform(np.max(control_rod.radii), fuel_outer_radius), 
                                                       phi=openmc.stats.Uniform(0, 2 * np.pi), 
                                                       z=openmc.stats.Uniform(-fuel_height / 2, fuel_height / 2))
    source.angle = openmc.stats.Isotropic()
    source.energy = openmc.stats.Discrete([1.0e6], [1.0])  # 1 MeV neutrons

    # Set up OpenMC settings
    settings = openmc.Settings()
    settings.batches = 100
    settings.inactive = 50
    settings.particles = 100000
    settings.source = source
    settings.export_to_xml()
    
    # Create voxel plot
    vox_plot = openmc.Plot()
    vox_plot.type = 'voxel'
    vox_plot.filename = 'vox_plot_{}.h5'.format(translation_height)
    
    vox_plot.colors = {
    moderator: 'blue',
    fuel: 'magenta',
    cladding: 'pink',
    control_rod_material: 'yellow'}
    vox_plot.color_by = 'material'
    vox_plot.width = (2*moderator_outer_radius, 2*moderator_outer_radius, fuel_height)
    vox_plot.pixels = (100, 100, 100)

    plots = openmc.Plots([vox_plot])
    plots.export_to_xml()
    
    # Plot geometry
    openmc.plot_geometry()
    
    return cylindrical_mesh.dimension, openmc.model.Model(geometry, mats, settings, tallies)

# Example coeffs for Legendre expansion
coeffs = [-0.01059336,  0.81348262, -0.13965943, -0.39435563,  0.62933526]

# Example of how to use the function
apf = []
h = []
for i in range(20,50,10):
    os.system('rm summary.h5')
    os.system('rm statepoint.100.h5')
    cylindrical_mesh_dims, mod = create_openmc_model(translation_height=i, coefficients = coeffs)
    h.append(i)
    mod.run()
    with openmc.StatePoint('statepoint.100.h5') as sp:
        tally_v = sp.get_tally(name='axial peaking')
        fission_rates = tally_v.get_values(scores=['flux'])
        fission_rates.shape = cylindrical_mesh_dims
        # Calculate the average fission rate along the height (z-axis)
        average_axial_power = np.zeros(cylindrical_mesh_dims[2])
        for i in range(0,cylindrical_mesh_dims[2]):
            average_axial_power[i] = np.mean(fission_rates[:,:,i])
        apf.append(np.max(average_axial_power/np.mean(average_axial_power)))
        
        print(apf)
    
plt.plot(h, apf)
