using Revise
println("revise loaded")
using Radiant
println("radiant loaded")
using Dates
println("dates loaded")
using QuadGK ############################
println("quadgk loaded")
using JLD2
println("jld2 loaded")

using PyCall
println("pycall loaded")
using PyPlot
pygui(true)
println("pyplot loaded")



case1 = ["no ria 10kevBW", "cs_noRIA_Au_0to1MeV_10KeVbinWidth_2026-07-21_23_14_06_341.jld2"]
case2 = ["ria",    "cs_RIA_Au_0to1MeV_10KeVbinWidth_2026-08-12_14_06_58_603.jld2"]
#case1 = ["no ria 1kevBW",    "cs_noRIA_Au_0to1MeV_1KeVbinWidth_2026-08-12_20_27_37_766.jld2"]
#case2 = ["no ria 5kevBW",    "cs_noRIA_Au_0to1MeV_5KeVbinWidth_2026-08-13_12_46_43_377.jld2"]

println("=== Graph simulation ===")

function run_case(version)

    println("=== Running case: $(version[1]) ===")

    label, filename = version[1], version[2]

    # User parameters ---------------
    target_thickness = 0.2 # cm
    number_of_voxels = 100

    println("Target thickness: $target_thickness cm")

    # Define particles ---------------
    electron = Electron()
    positron = Positron()
    photon = Photon()

    println("Defined particles: electron, positron, photon")

    # Load pre-built cross-sections ---------------
    data = JLD2.load(filename)
    cs = data["xs_full"]
    cs.set_particles([electron, photon, positron])

    println("Loaded cross-sections from: $filename")

    # Define material for geometry ---------------
    gold = Material("gold")
    gold.add_element("Au")

    println("Defined material: gold")

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

    println("Defined geometry with $number_of_voxels voxels in gold region of thickness $target_thickness cm")

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

    println("Defined solvers for electron, positron, and photon")

    # Define primary photon beam ----------------
    source = Surface_Source()
    source.set_particle(photon)
    source.set_intensity(1.0)
    source.set_energy_group(1) # number of energy group. = beam energy. try lower energies ================================================
    source.set_direction([1.0, 0.0, 0.0])
    source.set_location("X-")

    fixed_sources = Fixed_Sources(cs, geo, solvers)
    fixed_sources.add_source(source)

    println("Defined primary photon beam source at X- boundary")

    # Define computation settings --------------
    computation = Computation_Unit()
    computation.set_cross_sections(cs)
    computation.set_geometry(geo)
    computation.set_solvers(solvers)
    computation.set_sources(fixed_sources)
    println("=== Running case: $label ($filename) ===")
    @time computation.run()

    println("Computation completed for case: $label")

    # Extract outputs --------------
    x = computation.get_voxels_position("x")
    dose = computation.get_energy_deposition()
    
    # Photons
    energies = computation.get_energies(photon)
    photon_flux = computation.get_flux(photon)

    # Electrons
    electron_energies = computation.get_energies(electron)
    electron_flux = computation.get_flux(electron)

    # Pick the start / middle / end voxels for the spectrum plots ---------
    start_idx = 1
    mid_idx = cld(number_of_voxels, 2)
    end_idx = number_of_voxels
    spectrum_indices = [start_idx, mid_idx, end_idx]
    spectrum_positions = x[spectrum_indices]
    
    sampled_photon_flux = photon_flux[:, spectrum_indices]
    sampled_electron_flux = electron_flux[:, spectrum_indices]

    println("Extracted dose, photon flux, and electron flux data")

    return (
        label = label,
        filename = filename,
        x = x,
        dose = dose,
        energies = energies,
        electron_energies = electron_energies,
        spectrum_positions = spectrum_positions,
        spectrum_indices = spectrum_indices,
        photon_flux = sampled_photon_flux,
        electron_flux = sampled_electron_flux,
    )
end

function main()

    println("running main")
    println("case 1")
    case1_result = run_case(case1)
    println("case 2")
    case2_result = run_case(case2)

    # Dose relative difference in percent ---------------
    relative_diff = @. (case2_result.dose - case1_result.dose) / case1_result.dose * 100.0
    max_abs_idx = argmax(abs.(relative_diff))
    println("Max |relative dose difference| = $(round(relative_diff[max_abs_idx], digits=2))% at x = $(round(case2_result.x[max_abs_idx], digits=3)) cm")

    println("=== Dose comparison completed ===")

    timestamp = replace(string(now()), r"[:\.T]" => "_")
    stem = "Au_ria_vs_noria_compare_$(timestamp)"

    # Plot dose comparison --------------

    println("plotting dose comparison figure")

    figure(1)
    clf()

    subplot(211)
    title("Dose in gold")
    xlabel("Depth (cm)")
    ylabel("Energy deposition (MeV/g * cm^2)")
    plot(case2_result.x, case2_result.dose, label=case2_result.label)
    plot(case1_result.x, case1_result.dose, label=case1_result.label, linestyle="--")
    legend(fontsize=7)

    subplot(212)
    title("Relative dose difference")
    xlabel("Depth (cm)")
    ylabel("Difference (%)")
    plot(case2_result.x, relative_diff)
    axhline(0.0, color="gray", linewidth=0.8, linestyle=":")

    tight_layout()

    dose_figure_path = "$(stem)_dose_comparison.png"
    savefig(dose_figure_path, dpi=150)
    println("Saved dose figure to: ", dose_figure_path)

    # Plot photon spectra (log scale) at start / middle / end --------------

    println("plotting photon spectrum figures (log scale)")

    # Compton backscatter (mu = -1) for the source energy ---------------
    mₑc² = 0.510999 # MeV
    source_energy = maximum(case2_result.energies) # top of the energy grid = beam energy
    backscatter_energy = source_energy / (1 + 2 * (source_energy / mₑc²))
    println("Compton backscatter (mu=-1) for E0=$(round(source_energy, digits=3)) MeV: $(round(backscatter_energy, digits=4)) MeV")

    position_labels = ["Start of material", "Middle of material", "End of material"]

    figure(2)
    clf()

     for i in 1:3
        subplot(3, 1, i)
        title("$(position_labels[i]) (x ~ $(round(case2_result.spectrum_positions[i], digits=3)) cm)")
        xlabel("Energy (MeV)")
        ylabel("Photon flux (log scale)")
        semilogy(case2_result.energies, case2_result.photon_flux[:, i], label=case2_result.label)
        semilogy(case1_result.energies, case1_result.photon_flux[:, i], label=case1_result.label, linestyle="--")
        axvline(backscatter_energy, color="gray", linewidth=1.2, linestyle=":", label="Backscatter ($(round(backscatter_energy, digits=3)) MeV)")
        legend(fontsize=7)
    end

    tight_layout()

    spectrum_figure_path = "$(stem)_spectrum_log.png"
    savefig(spectrum_figure_path, dpi=150)
    println("Saved spectrum figure to: ", spectrum_figure_path)

    # Plot electron spectra (log scale) at start / middle / end --------------

    println("plotting electron spectrum figures (log scale)")

    figure(3)
    clf()

    for i in 1:3
        subplot(3, 1, i)
        title("Electron Flux: $(position_labels[i]) (x ~ $(round(case2_result.spectrum_positions[i], digits=3)) cm)")
        xlabel("Energy (MeV)")
        ylabel("Electron flux (log scale)")
        semilogy(case2_result.electron_energies, case2_result.electron_flux[:, i], label=case2_result.label)
        semilogy(case1_result.electron_energies, case1_result.electron_flux[:, i], label=case1_result.label, linestyle="--")
        legend(fontsize=7)
    end

    tight_layout()

    electron_figure_path = "$(stem)_electron_spectrum_log.png"
    savefig(electron_figure_path, dpi=150)
    println("Saved electron spectrum figure to: ", electron_figure_path)
  
    PyPlot.show()

    return nothing
end

Base.invokelatest(main)
