"""
trying multigroup discretization integration

integral of final energy in in the middle test

definitions:

    - E : initial photon energy in mₑc²
    - Ef : outgoing photon energy in mₑc²
    - Ui : binding energy of subshell i in mₑc²
    - Jio : part of orbital compton profile of subshell i in mₑc²
    - l : Legendre order
    - g : incoming photon energy group index
    - g_f : outgoing photon energy group index
    - Z : atomic number (Al)
    - E_max : (keV) midpoint energy of the highest energy group.
    - E_min : (keV) cutoff energy for multigroup structure

"""

using Plots
using QuadGK
using Radiant

include("Radiant.jl/src/cross_sections/electron_subshells.jl")
include("Radiant.jl/src/cross_sections/orbital_compton_profiles.jl")
include("Radiant.jl/src/cross_sections/energy_group_structure.jl")
include("recurrent_functions.jl") # my functions


Z = 13 # atomic number (Al)
E_max = 50. # keV. midpoint energy of the highest energy group.
E_min = 0.01 # keV. cutoff energy for multigroup structure
Ef_max = 50.
Ef_min = 0.01

g = 2 # initial photon energy group index
g_f = 2 # outgoing photon energy group index
l = 0 # Legendre moment order

E_max = keV_to_mec2(E_max) # keV → mₑc²
E_min = keV_to_mec2(E_min) # keV → mₑc²
Ef_max = keV_to_mec2(Ef_max) # keV → mₑc²
Ef_min = keV_to_mec2(Ef_min) # keV → mₑc²

_, Zi, Ui, _, _, _ = electron_subshells(Z,false)
Jio = inverse_mₑc²_units_orbital_compton_profiles(Z)



""" Radiant code:

    linear_energy_group_structure(Ng::Int64,E::Number,Ec::Number)

Generate a linear group structure for multigroup calculations.

# Input Argument(s)
- `Ng::Int64`: number of groups.
- `E::Number`: midpoint energy of the highest energy group.
- `Ec::Number`: cutoff energy.

# Output Argument(s)
- `Eᵇ::Vector{Float64}`: vector containing the (Ng+1)-boundaries of the group structure.

"""

E_grid = linear_energy_group_structure(10, E_max, E_min)
Ef_grid = linear_energy_group_structure(10, Ef_max, Ef_min)

function print_E_structure(E_structure)
    println("Energy group structure (keV):")
    for i in 1:length(E_structure)-1
        println("Group $i: [$(mec2_to_keV(E_structure[i])), $(mec2_to_keV(E_structure[i+1]))] keV")
    end
end

print_E_structure(E_grid)
print_E_structure(Ef_grid)



function Ef_integral(l, g_f, Ef_grid, E, Ui, Jio)
    
    group_index_upper_bound_f = g_f
    group_index_lower_bound_f = g_f + 1

    lower_bound = Ef_grid[group_index_lower_bound_f]

    upper_bound = min(Ef_grid[group_index_upper_bound_f], E)

    int, err = quadgk(Ef -> begin
        legendre_moment(l, E, Ef, Ui, Jio)
    end, lower_bound, upper_bound)
    
    return int

end


function E_integral(l, g, g_f, E_grid, Ef_grid, Ui, Jio)
    
    group_index_upper_bound_f = g
    group_index_lower_bound_f = g + 1

    lower_bound = E_grid[group_index_lower_bound_f]

    upper_bound = E_grid[group_index_upper_bound_f]

    int, err = quadgk(E -> begin
        Ef_integral(l, g_f, Ef_grid, E, Ui, Jio)
    end, lower_bound, upper_bound)

    return int

end


function deltaEf(g_f, Ef_grid)

    return Ef_grid[g_f] - Ef_grid[g_f+1]

end

group_legendre_moment = E_integral(l, g, g_f, E_grid, Ef_grid, Ui, Jio) / deltaEf(g_f, Ef_grid)

println("Group Legendre moment for l=$l, g=$g and g_f=$g_f: $group_legendre_moment")
