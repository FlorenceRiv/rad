

using QuadGK
using Radiant

include("Radiant.jl/src/cross_sections/electron_subshells.jl")
include("Radiant.jl/src/cross_sections/orbital_compton_profiles.jl")
include("Radiant.jl/src/cross_sections/energy_group_structure.jl")
include("recurrent_functions.jl") 

"""
For compton scattering with impulse approximation

      only doing photons so far. no results for electron (type='P')

Calculate the feed function 𝓕 (normalized probability of scattering from Ei into each
group gf) for each Legendre moment up to order L. Also calculate the energy weighted
feed function 𝓕ₑ for energy-deposition cross section.

# Input Argument(s)
- `Z::Vector{Int64}` : atomic number of the element(s) composing the material.
- `ωz::Vector{Float64}` : weight fraction of the element(s) composing the material.
- `ρ::Float64` : density of the material [in g/cm³].
- `L::Int64` : Legendre truncation order.
- `Ei::Float64` : energy of the incoming particle [in mₑc²].
- `Eout::Vector{Float64}` : energy group boundaries [in mₑc²].
- `Ng::Int64` : number of groups.
- `particles::Vector{Particle}` : list of the particles imply in the interaction.
- `type::String` : type of interaction (scattering 'S' or production 'P').
- `incoming_particle::Particle` : incoming particle.
- `scattered_particle::Particle` : scattered particle.


# Output Argument(s)
- `𝓕::Array{Float64}` : feed function.
- `𝓕ₑ::Vector{Float64}` : energy weighted feed function.

"""

function feed_ria(
    Z::Vector{Int64},
    ωz::Vector{Float64},
    ρ::Float64,
    L::Int64,
    Ei::Float64,
    Eout::Vector{Float64},
    Ng::Int64#,
#    particles::Vector{Particle},
#    type::String,
#    incoming_particle::Particle,
#    scattered_particle::Particle,
)

interaction = Compton() # hardcoding instead of input argument. this struct has default parameters
interaction.model = "impulse_approximation" # default is "waller-hartree"
interaction.is_subshells_dependant = true # ria needs subshell depedency
type = "S" # only doing photon for now

𝓕 = zeros(Ng+1,L+1)
𝓕ₑ = zeros(Ng+1)


Nz = length(Z)
for i in 1:Nz # for every element

    Nshells,Zi,Ui,Ti,ri,Jio = electron_subshells(Z[i],false) # electrons not free

    macro_factor = nuclei_density(Z[i], ρ) * ωz[i]

    for gf in 1:Ng # Iterate through each final energy group
 
        # Final energy group
        Ef⁻ = Eout[gf]; Ef⁺ = Eout[gf+1]
        Ef⁻,Ef⁺,isSkip = bounds(interaction,Ef⁻,Ef⁺,Ei,type)
        if isSkip continue end
        

        # integrate over final energy group
        group_moments, _ = quadgk(Ef⁺, Ef⁻) do Ef
        
            # integrate legendre moment (ddcs over mu) for each l :
            # ddcs_ria goes over all subshells
            dcs_moment_array, _ = quadgk(-1.0, 1.0) do μ
                return Float64[ddcs_ria(Ei, Ef, Zi, Ui, Jio, μ) * Radiant.legendre_polynomials(l, μ) for l in 0:L]
            end
            dcs_moment_array .*= 2π 

            # energy-weighted L=0 moment for 𝓕ₑ
            energy_weighted = dcs_moment_array[1] * Ef
            
            # push the energy-weighted value to the end of the vector
            # allows quadgk to integrate 𝓕 and 𝓕ₑ in a single pass
            return push!(dcs_moment_array, energy_weighted)

        end

        # Legendre moments (indices 1 to L+1) and scale
        𝓕[gf, :] .+= group_moments[1:L+1] .* macro_factor

        # energy-weighted value (index L+2) and scale
        𝓕ₑ[gf] += group_moments[L+2] * macro_factor

    end
end

return 𝓕, 𝓕ₑ

end
