"""
trying multigroup discretization integration

integral of final energy in in the middle test

definitions:

    - E : initial photon energy in mₑc²
    - Ef : outgoing photon energy in mₑc²
    - Ui : binding energy of subshell i in mₑc²
    - Jio : part of orbital compton profile of subshell i in mₑc²
    - l : Legendre moment order
    - g : initial photon energy group index
    - g_f : outgoing photon energy group indexa
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
include("recurrent_functions.jl") 

Z = 13 # (Al)
E_max_keV = 50.
E_min_keV = 0.01
Ef_max_keV = 50. 
Ef_min_keV = 0.01

g = 2
g_f = 2
l = 100

E_max = keV_to_mec2(E_max_keV)
E_min = keV_to_mec2(E_min_keV)
Ef_max = keV_to_mec2(Ef_max_keV)
Ef_min = keV_to_mec2(Ef_min_keV)

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

#print_E_structure(E_grid)
print_E_structure(Ef_grid)



function Ef_integral(l, g_f, Ef_grid, E, Zi, Ui, Jio)
    
    if E < Ef_grid[g_f + 1]
        println("Initial photon energy E=$(mec2_to_keV(E)) keV is less than the lower bound Ef=$(mec2_to_keV(Ef_grid[g_f + 1])) keV of outgoing photon energy group g'=$g_f .")
        return 0.0
    end
        
    group_index_upper_bound_f = g_f
    group_index_lower_bound_f = g_f + 1

    Ef_lower_bound = Ef_grid[group_index_lower_bound_f]

    # final photon energy cannot exceed initial photon energy
    Ef_upper_bound = min(Ef_grid[group_index_upper_bound_f], E)

    int, err = quadgk(Ef -> begin
        legendre_moment(l, E, Ef, Zi, Ui, Jio)
    end, Ef_lower_bound, Ef_upper_bound)
    
    return int

end


function E_integral(l, g, g_f, E_grid, Ef_grid, Zi, Ui, Jio)
    
    group_index_upper_bound = g
    group_index_lower_bound = g + 1

    E_lower_bound = E_grid[group_index_lower_bound]

    E_upper_bound = E_grid[group_index_upper_bound]

    int, err = quadgk(E -> begin
        Ef_integral(l, g_f, Ef_grid, E, Zi, Ui, Jio)
    end, E_lower_bound, E_upper_bound)

    return int

end


function deltaEf(g, E_grid)

    return E_grid[g] - E_grid[g+1]

end

#group_legendre_moment = E_integral(l, g, g_f, E_grid, Ef_grid, Zi, Ui, Jio) / deltaEf(g, E_grid)
#println("Group Legendre moment for l=$l, g=$g and g_f=$g_f: $group_legendre_moment")

Specific_energy_legendre_moment = Ef_integral(l, g_f, Ef_grid, keV_to_mec2(50), Zi, Ui, Jio)
println("Legendre moment for l=$l, g_f=$g_f and E=50 keV: $Specific_energy_legendre_moment")
