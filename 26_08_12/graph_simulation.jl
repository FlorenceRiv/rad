using Revise
using Radiant
using DelimitedFiles
using Dates
using JLD2

using PyCall
using PyPlot


no_ria = ["no ria", "cs_noRIA_Au_0to1MeV_10KeVbinWidth_2026-07-21_23_14_06_341.jld2"]
ria    = ["ria",    "cs_RIA_Au_0to1MeV_10KeVbinWidth_2026-07-22_01_34_30_636.jld2"]


function run_case(version)

    label, filename = version[1], version[2]

    # User parameters ---------------
    target_thickness = 0.5 # cm
    number_of_voxels = 100
    number_of_spectrum_positions = 10

    # Define particles ---------------
    electron = Electron()
    positron = Positron()
    photon = Photon()

    # Load pre-built cross-sections ---------------
    data = JLD2.load(filename)
    cs = data["xs_full"]
    cs.set_particles([electron, photon, positron])

    # Define material for geometry ---------------
    gold = Material("gold")
    gold.add_element("Au")

    # Define geometry ---------------
    geo = Geometry()
    geo.set_type("Cartesian")
    geo.set_dimension(1)
    geo.set_material_per_region([gold])
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
    println("=== Running case: $label ($filename) ===")
    @time computation.run()

    # Extract outputs --------------
    x = computation.get_voxels_position("x")
    dose = computation.get_energy_deposition()
    energies = computation.get_energies(photon)
    photon_flux = computation.get_flux(photon)

    spectrum_indices = unique(round.(Int, range(1, size(photon_flux, 2), length=number_of_spectrum_positions)))
    spectrum_positions = x[spectrum_indices]
    sampled_photon_flux = photon_flux[:, spectrum_indices]

    # Save this run's own CSVs, same as before --------------
    timestamp = replace(string(now()), r"[:\.T]" => "_")
    top_energy = round(maximum(energies), digits=3)
    stem = "Au_$(replace(filename, ".jld2" => ""))_E$(sanitize(top_energy))MeV_$(timestamp)"

    dose_path = "$(stem)_dose.csv"
    spectrum_path = "$(stem)_photon_spectrum.csv"

    open(dose_path, "w") do io
        write(io, "x_cm,dose_MeV_per_g_cm2\n")
        writedlm(io, [x dose], ',')
    end

    open(spectrum_path, "w") do io
        spectrum_header = join(["energy_MeV"; ["flux_x$(sanitize(position))cm" for position in spectrum_positions]], ",")
        write(io, spectrum_header * "\n")
        writedlm(io, [energies sampled_photon_flux], ',')
    end

    println("Saved dose to: ", dose_path)
    println("Saved photon spectrum to: ", spectrum_path)

    return (
        label = label,
        filename = filename,
        x = x,
        dose = dose,
        energies = energies,
        spectrum_positions = spectrum_positions,
        spectrum_indices = spectrum_indices,
        photon_flux = sampled_photon_flux,
    )
end

function main()

    no_ria_result = run_case(no_ria)
    ria_result = run_case(ria)

    # Dose relative difference (RIA vs no-RIA), in percent ---------------
    relative_diff = @. (ria_result.dose - no_ria_result.dose) / no_ria_result.dose * 100.0
    max_abs_idx = argmax(abs.(relative_diff))
    println("Max |relative dose difference| = $(round(relative_diff[max_abs_idx], digits=2))% at x = $(round(ria_result.x[max_abs_idx], digits=3)) cm")

    # Save comparison table --------------
    timestamp = replace(string(now()), r"[:\.T]" => "_")
    stem = "Au_ria_vs_noria_compare_$(timestamp)"
    dose_compare_path = "$(stem)_dose.csv"

    open(dose_compare_path, "w") do io
        write(io, "x_cm,dose_ria_MeV_per_g_cm2,dose_no_ria_MeV_per_g_cm2,relative_diff_percent\n")
        writedlm(io, [ria_result.x ria_result.dose no_ria_result.dose relative_diff], ',')
    end
    println("Saved dose comparison to: ", dose_compare_path)

    # Plot comparison --------------
    figure(1)
    clf()

    subplot(311)
    title("Dose in gold: RIA vs no-RIA")
    xlabel("Depth (cm)")
    ylabel("Energy deposition (MeV/g * cm^2)")
    plot(ria_result.x, ria_result.dose, label=ria_result.label)
    plot(no_ria_result.x, no_ria_result.dose, label=no_ria_result.label, linestyle="--")
    legend(fontsize=7)

    subplot(312)
    title("Relative dose difference (RIA vs no-RIA)")
    xlabel("Depth (cm)")
    ylabel("Difference (%)")
    plot(ria_result.x, relative_diff)
    axhline(0.0, color="gray", linewidth=0.8, linestyle=":")

    subplot(313)
    pos_idx = ceil(Int, length(ria_result.spectrum_positions) / 2)
    title("Photon spectrum at x = $(round(ria_result.spectrum_positions[pos_idx], digits=3)) cm")
    xlabel("Energy (MeV)")
    ylabel("Photon flux")
    plot(ria_result.energies, ria_result.photon_flux[:, pos_idx], label=ria_result.label)
    plot(no_ria_result.energies, no_ria_result.photon_flux[:, pos_idx], label=no_ria_result.label, linestyle="--")
    legend(fontsize=7)

    tight_layout()

    figure_path = "$(stem)_comparison.png"
    savefig(figure_path, dpi=150)
    println("Saved figure to: ", figure_path)


    PyPlot.show()
    

    return nothing
end

Base.invokelatest(main)