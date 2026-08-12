using Revise
using Radiant
using DelimitedFiles
using Dates
using QuadGK

using PyCall
using PyPlot
pygui(true)

function run_case(model_name::String, is_subshells_dependant::Bool;
                   beam_energy::Float64=0.1,
                   minimum_energy::Float64=0.001,
                   number_of_groups::Int=20,
                   group_structure::String="linear",
                   target_thickness::Float64=10.0,
                   number_of_voxels::Int=100,
                   number_of_spectrum_positions::Int=10)

    println("=== Running Compton model: $model_name (subshells=$is_subshells_dependant) ===")

    # Define materials ---------------
    aluminum = Material("aluminum")
    aluminum.add_element("Al")

    # Define particles ---------------
    electron = Electron()
    positron = Positron()
    photon = Photon()

    # Compton interaction for this case ---------------
    compton_interaction = Compton()
    compton_interaction.model = model_name
    compton_interaction.is_subshells_dependant = is_subshells_dependant

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
        compton_interaction,
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

    # Extract outputs --------------
    x = computation.get_voxels_position("x")
    dose = computation.get_energy_deposition()
    energies = computation.get_energies(photon)
    photon_flux = computation.get_flux(photon)

    spectrum_indices = unique(round.(Int, range(1, size(photon_flux, 2), length=number_of_spectrum_positions)))
    spectrum_positions = x[spectrum_indices]
    sampled_photon_flux = photon_flux[:, spectrum_indices]

    return (
        model = model_name,
        x = x,
        dose = dose,
        energies = energies,
        spectrum_positions = spectrum_positions,
        spectrum_indices = spectrum_indices,
        photon_flux = sampled_photon_flux,
    )
end

sanitize(value) = replace(string(value), "." => "p")

function main()

    # User parameters ---
    beam_energy = 0.1 # MeV
    minimum_energy = 0.001 # MeV
    number_of_groups = 80
    group_structure = "linear"
    target_thickness = 10.0 # cm
    number_of_voxels = 100
    number_of_spectrum_positions = 10

    common_kwargs = (
        beam_energy=beam_energy, minimum_energy=minimum_energy,
        number_of_groups=number_of_groups, group_structure=group_structure,
        target_thickness=target_thickness, number_of_voxels=number_of_voxels,
        number_of_spectrum_positions=number_of_spectrum_positions,
    )

    # RIA run ---------------
    ria = run_case("impulse_approximation", true; common_kwargs...)

    # Klein-Nishina run ---------------
    reference = run_case("klein-nishina", false; common_kwargs...)

    # Dose relative difference in percent ---------------
    relative_diff = @. (ria.dose - reference.dose) / reference.dose * 100.0
    max_abs_idx = argmax(abs.(relative_diff))
    println("Max |relative dose difference| = $(round(relative_diff[max_abs_idx], digits=2))% at x = $(round(ria.x[max_abs_idx], digits=3)) cm")

    # Save comparison table --------------
    timestamp = replace(string(now()), r"[:\.T]" => "_")
    stem = "compton_compare_$(sanitize(ria.model))_vs_$(sanitize(reference.model))_$(timestamp)"
    dose_path = "$(stem)_dose.csv"

    open(dose_path, "w") do io
        write(io, "x_cm,dose_$(ria.model)_MeV_per_g_cm2,dose_$(reference.model)_MeV_per_g_cm2,relative_diff_percent\n")
        writedlm(io, [ria.x ria.dose reference.dose relative_diff], ',')
    end
    println("Saved dose comparison to: ", dose_path)

    # Plot comparison --------------
    figure(1)
    clf()

    subplot(311)
    title("Dose comparison")
    xlabel("Depth (cm)")
    ylabel("Energy deposition (MeV/g * cm^2)")
    PyPlot.plot(ria.x, ria.dose, label=ria.model)
    PyPlot.plot(reference.x, reference.dose, label=reference.model, linestyle="--")
    legend(fontsize=7)

    subplot(312)
    title("Relative dose difference (RIA vs reference)")
    xlabel("Depth (cm)")
    ylabel("Difference (%)")
    PyPlot.plot(ria.x, relative_diff)
    PyPlot.axhline(0.0, color="gray", linewidth=0.8, linestyle=":")

    subplot(313)
    pos_idx = ceil(Int, number_of_spectrum_positions / 2)
    title("Photon spectrum at x = $(round(ria.spectrum_positions[pos_idx], digits=2)) cm")
    xlabel("Energy (MeV)")
    ylabel("Photon flux")
    PyPlot.plot(ria.energies, ria.photon_flux[:, pos_idx], label=ria.model)
    PyPlot.plot(reference.energies, reference.photon_flux[:, pos_idx], label=reference.model, linestyle="--")
    legend(fontsize=7)

    tight_layout()
    PyPlot.show()

    return nothing
end

#main()
Base.invokelatest(main)