using Revise
using Radiant
using DelimitedFiles
using Dates
using QuadGK
using JLD2

# Define materials ---------------
gold = Gold()

# Define particles ---------------
electron = Electron()
positron = Positron()
photon = Photon()

# set group structure
number_of_groups = 1000
upper_Ei_midpoint = 1 # MeV
cutoff_E = 0.001
bin_width = (upper_Ei_midpoint - cutoff_E) / (number_of_groups - 0.5)
println("bin width = ", bin_width, "MeV")

cs = Cross_Sections()
cs.set_source("physics-models")
cs.set_materials([gold])
cs.set_particles([electron,photon,positron])
cs.set_group_structure("linear",number_of_groups,upper_Ei_midpoint,cutoff_E) 

# Initalize Compton set ria model ---------------
compton_interaction = Compton() 
compton_interaction.model = "impulse_approximation"
compton_interaction.is_subshells_dependant = true # ria needs subshell depedency

cs.set_interactions([compton_interaction])
cs.set_legendre_order(7)
cs.build()
JLD2.save("test_only_compton.jld2", Dict("xs_full" => cs))

println("bin width = ", bin_width, "MeV")



