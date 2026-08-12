using Revise
using Radiant
using DelimitedFiles
using Dates
using QuadGK ############################

using PyCall
using PyPlot
pygui(true)

function main()

    # User parameters ---------------
    beam_energy = 0.1 # MeV
    minimum_energy = 0.001 # MeV
    number_of_groups = 20 ######### was 80. changed to 20 for faster runtime
    group_structure = "linear"
    target_thickness = 10.0 # cm
    number_of_voxels = 100
    number_of_spectrum_positions = 10

    # Define materials ---------------
    water = Material("water")
    water.set_density(1.0)
    water.set_state_of_matter("liquid")
    water.add_element("H", 0.1111)
    water.add_element("O", 0.8889)

    aluminum = Material("aluminum")
    aluminum.add_element("Al")

    # Define particles ---------------
    electron = Electron()
    positron = Positron()
    photon = Photon()

    # Initalize Compton set ria model --------------- #############################
    compton_interaction = Compton()
    compton_interaction.model = "impulse_approximation"
    compton_interaction.is_subshells_dependant = true # ria needs subshell depedency

    # Define cross-sections ---------------
    cs = Cross_Sections()
    cs.set_source("physics-models")
    cs.set_materials([aluminum])
    cs.set_particles([electron, photon, positron])
    cs.set_group_structure(group_structure, number_of_groups, beam_energy, minimum_energy)
    cs.set_interactions([
        Elastic_Collision(),
        Inelastic_Collision(),
        Bremsstrahlung(),
        Annihilation(),
        compton_interaction, #############################
        Rayleigh(),
        Photoelectric(),
        Pair_Production(),
        Relaxation(),
    ])
    cs.set_legendre_order(7)
    cs.build()

    # Define geometry ---------------
    geo = Geometry()
    geo.set_type("Cartesian")
    geo.set_dimension(1)
    geo.set_material_per_region([aluminum])
    geo.set_boundary_conditions("X-", "void")
    geo.set_boundary_conditions("X+", "void")
    geo.set_number_of_regions("X", 1)
    geo.set_voxels_per_region("X", [number_of_voxels])
    geo.set_region_boundaries("X", [0.0, target_thickness])
    geo.build(cs)

    # Define methods ----------------
    electron_solver = Discrete_Ordinates()
    electron_solver.set_particle(electron)
    electron_solver.set_solver_type("BFP")
    electron_solver.set_acceleration("livolant")
    electron_solver.set_quadrature("gauss-lobatto", 8)
    electron_solver.set_angular_boltzmann("galerkin")
    electron_solver.set_angular_fokker_planck("finite-difference")
    electron_solver.set_convergence_criterion(1e-7)
    electron_solver.set_maximum_iteration(300)
    electron_solver.set_scheme("E", "DG", 2)
    electron_solver.set_scheme("x", "DG", 2)
    electron_solver.set_is_full_coupling(false)

    positron_solver = Discrete_Ordinates()
    positron_solver.set_particle(positron)
    positron_solver.set_solver_type("BFP")
    positron_solver.set_acceleration("livolant")
    positron_solver.set_quadrature("gauss-lobatto", 8)
    positron_solver.set_angular_boltzmann("galerkin")
    positron_solver.set_angular_fokker_planck("finite-difference")
    positron_solver.set_convergence_criterion(1e-7)
    positron_solver.set_maximum_iteration(300)
    positron_solver.set_scheme("E", "DG", 2)
    positron_solver.set_scheme("x", "DG", 2)
    positron_solver.set_is_full_coupling(false)

    photon_solver = Discrete_Ordinates()
    photon_solver.set_particle(photon)
    photon_solver.set_solver_type("BTE")
    photon_solver.set_acceleration("livolant")
    photon_solver.set_quadrature("gauss-lobatto", 8)
    photon_solver.set_angular_boltzmann("galerkin")
    photon_solver.set_convergence_criterion(1e-7)
    photon_solver.set_maximum_iteration(300)
    photon_solver.set_scheme("E", "DG", 2)
    photon_solver.set_scheme("x", "DG", 2)

    solvers = Solvers()
    solvers.add_solver(electron_solver)
    solvers.add_solver(positron_solver)
    solvers.add_solver(photon_solver)
    solvers.set_maximum_number_of_generations(20)
    solvers.set_convergence_criterion(1e-7)

    # Define primary photon beam ----------------
    source = Surface_Source()
    source.set_particle(photon)
    source.set_intensity(1.0)
    source.set_energy_group(1)
    source.set_direction([1.0, 0.0, 0.0])
    source.set_location("X-")

    fixed_sources = Fixed_Sources(cs, geo, solvers)
    fixed_sources.add_source(source)

    # Define computation settings --------------
    computation = Computation_Unit()
    computation.set_cross_sections(cs)
    computation.set_geometry(geo)
    computation.set_solvers(solvers)
    computation.set_sources(fixed_sources)
    @time computation.run()

    # Output dose and photon spectrum --------------
    x = computation.get_voxels_position("x")
    dose = computation.get_energy_deposition()
    energies = computation.get_energies(photon)
    photon_flux = computation.get_flux(photon)

    spectrum_indices = unique(round.(Int, range(1, size(photon_flux, 2), length=number_of_spectrum_positions)))
    spectrum_positions = x[spectrum_indices]
    sampled_photon_flux = photon_flux[:, spectrum_indices]

    sanitize(value) = replace(string(value), "." => "p")
    timestamp = replace(string(now()), r"[:\.T]" => "_")
    stem = "general_1D_photon_beam_E$(sanitize(beam_energy))MeV_$(sanitize(compton_interaction.model))_$(timestamp)"

    dose_path = "$(stem)_dose.csv"
    spectrum_path = "$(stem)_photon_spectrum.csv"

    open(dose_path, "w") do io
        write(io, "x_cm,dose_MeV_per_g_cm2\n")
        writedlm(io, [x dose], ',')
    end

    #open(spectrum_path, "w") do io
    #    spectrum_header = join(["energy_MeV"; ["flux_x$(sanitize(position))cm" for position in spectrum_positions]], ",")
    #    write(io, spectrum_header * "\n")
    #    writedlm(io, [energies sampled_photon_flux], ',')
    #end

    #println("Saved dose to: ", dose_path)
    #println("Saved photon spectrum to: ", spectrum_path)

    figure(1)
    clf()
    subplot(211)
    title("Dose")
    xlabel("Depth (cm)")
    ylabel("Energy deposition (MeV/g * cm^2)")
    PyPlot.plot(x, dose)

    subplot(212)
    title("Photon spectrum")
    xlabel("Energy (MeV)")
    ylabel("Photon flux")
    for i in 1:length(spectrum_indices)
        PyPlot.plot(energies, sampled_photon_flux[:, i], label="x=$(round(spectrum_positions[i], digits=2)) cm")
    end
    legend(fontsize=7)
    tight_layout()
    PyPlot.show()

    return nothing
end

#main()
Base.invokelatest(main)
