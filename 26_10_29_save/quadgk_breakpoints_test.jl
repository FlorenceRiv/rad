using Revise
using Radiant
using DelimitedFiles
using Dates
using QuadGK
using JLD2

# Define materials ---------------
water = Radiant.Water()
air = Air_Dry_Near_Sea_Level()
aluminum = Aluminum()
tungsten = Tungsten()
gold = Gold()

# Define particles ---------------
electron = Electron()
positron = Positron()
photon = Photon()

# Set group structure
number_of_groups = 10
upper_Ei_midpoint = 0.01 # MeV
cutoff_E = 0.001
bin_width = (upper_Ei_midpoint - cutoff_E) / (number_of_groups - 0.5)

function make_cs_object()
    cs = Cross_Sections()
    cs.set_source("physics-models")
    cs.set_materials([gold])
    cs.set_particles([electron, photon, positron])
    cs.set_group_structure("linear", number_of_groups, upper_Ei_midpoint, cutoff_E)

    compton_interaction = Compton()
    compton_interaction.model = "impulse_approximation"
    compton_interaction.is_subshells_dependant = true

    cs.set_interactions([compton_interaction])
    cs.set_legendre_order(7)
    return cs
end

# ------------------------------------------------------------
# 1. Warmup Pass: Force JIT compilation before taking measurements
# ------------------------------------------------------------
println("JIT compilation warmup...")
Radiant.TEST_BREAKPOINT[] = false
make_cs_object().build()
Radiant.TEST_BREAKPOINT[] = true
make_cs_object().build()
Radiant.TEST_BREAKPOINT[] = false   # reset before the real runs
GC.gc()
println("=============================================================================================")

# ------------------------------------------------------------
# 2. Run 1: No Breakpoint
# ------------------------------------------------------------
Radiant.TEST_BREAKPOINT[] = false  # Set module-level variable explicitly

cs1 = make_cs_object()
time_no_breakpoint = @elapsed cs1.build() 

GC.gc() # Clear heap memory before next run

# ------------------------------------------------------------
# 3. Run 2: With Breakpoint
# ------------------------------------------------------------
Radiant.TEST_BREAKPOINT[] = true

cs2 = make_cs_object()
time_with_breakpoint = @elapsed cs2.build()

# ------------------------------------------------------------
# Output Results
# ------------------------------------------------------------
println("no breakpoint:   $(round(time_no_breakpoint, digits=4)) seconds")
println("with breakpoint: $(round(time_with_breakpoint, digits=4)) seconds")