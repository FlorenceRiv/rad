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
titanium = Titanium()
cobalt = Cobalt()
copper = Copper()
nickel = Nickel()
niobium = Niobium()
thulium = Thulium()
yttrium = Yttrium()
manganese = Manganese()
uranium = Uranium()

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
cs.set_materials([#water,air,aluminum,tungsten,
                    gold
                    #,titanium,cobalt,copper,nickel,niobium,thulium,yttrium,manganese,uranium
                ])
cs.set_particles([electron,photon,positron])
cs.set_group_structure("linear",number_of_groups,upper_Ei_midpoint,cutoff_E) 
# linear, number of groups, midpoint energy of the highest-energy group (MeV), cutoff energy (lower bound of the lowest-energy group, MeV)

# Initalize Compton set ria model ---------------
compton_interaction = Compton() 
compton_interaction.model = "impulse_approximation"
compton_interaction.is_subshells_dependant = true # ria needs subshell depedency

cs.set_interactions([Elastic_Collision(),Inelastic_Collision(),Bremsstrahlung(),Annihilation(),compton_interaction,Rayleigh(),Photoelectric(),Pair_Production(),Relaxation()])
cs.set_legendre_order(7)
cs.build()
JLD2.save("cross_sections_RIA_Au_0to1MeV_1KeVbinWidth.jld2", Dict("xs_full" => cs))

println("bin width = ", bin_width, "MeV")
