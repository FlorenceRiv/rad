
"""
    feed(Z::Vector{Int64},ωz::Vector{Float64},ρ::Float64,L::Int64,Ei::Float64,
    Eout::Vector{Float64},Ng::Int64,interaction::Interaction,gi::Int64,Ngi::Int64,
    particles::Vector{Particle},type::String,incoming_particle::Particle,
    scattered_particle::Particle,Ein::Vector{Float64},Ec::Float64)

Calculate the feed function 𝓕 (normalized probability of scattering from Ei into each
group gf) for each Legendremoments up to order L. Also calculate the energy weighted
feed function 𝓕ₑ for energy-deposition cross section.

# Input Argument(s)
- `Z::Vector{Int64}` : atomic number of the element(s) composing the material.
- `ωz::Vector{Float64}` : weight fraction of the element(s) composing the material.
- `ρ::Float64` : density of the material [in g/cm³].
- `L::Int64` : Legendre truncation order.
- `Ei::Float64` : energy of the incoming particle [in mₑc²].
- `Eout::Vector{Float64}` : energy group boundaries [in mₑc²].
- `Ng::Int64` : number of groups.
- `interaction::Interaction` : interaction informations.
- `gi::Int64` : incoming particle group index.
- `Ngi::Int64` :  number of groups for the incoming particle.
- `particles::Vector{Particle}` : list of the particles imply in the interaction.
- `type::String` : type of interaction (scattering or production).
- `incoming_particle::Particle` : incoming particle.
- `scattered_particle::Particle` : scattered particle.
- `Ein::Vector{Float64}` : energy group boundaries corresponding to the incoming
  particle [in mₑc²].
- `Ec::Float64` : cutoff energy between soft and catastrophic interaction.
- `is_elastic_scattering::Bool` : boolean indicating if the scattering is elastic.
- `is_subshells::Bool` : boolean indicating if the cross-sections are subshells dependant.

# Output Argument(s)
- `𝓕::Array{Float64}` : feed function.
- `𝓕ₑ::Vector{Float64}` : energy weighted feed function.

# Reference(s)
- MacFarlane et al. (2021) : The NJOY Nuclear Data Processing System, Version 2012.

"""
function feed(Z::Vector{Int64},ωz::Vector{Float64},ρ::Float64,L::Int64,Ei::Float64,Eout::Vector{Float64},Ng::Int64,interaction::Interaction,gi::Int64,Ngi::Int64,particles::Vector{Particle},type::String,incoming_particle::Particle,scattered_particle::Particle,Ein::Vector{Float64},Ec::Float64,is_elastic::Bool,is_subshells::Bool)

#----
# Initialization
#----
𝓕 = zeros(Ng+1,L+1)
𝓕ₑ = zeros(Ng+1)

if ! ((interaction isa Compton) && (interaction.model == "impulse_approximation")) # check for special case of compton ria, otherwise do the regular feed calculation

    ΔQ = get_mass_energy_variation(interaction,type,true)

    # Outgoing particle energy spectrum
    is_dirac, Np, q_type = out_distribution_dispatch(interaction,type)
    if is_dirac Np = 1; u = [0]; w = [2] else u,w = quadrature(Np,q_type) end

    #----
    # Feed function over all groups and under the cutoff energy
    #----

    # Loop over the coumpound elements
    Nz = length(Z)
    for i in range(1,Nz)

        # Loop over subshells
        Nshells,Zi,Ui,Ti,ri,_ = electron_subshells(Z[i],~is_subshells)
        for gf in range(1,Ng), δi in range(1,Nshells)
            
            # Final energy group
            Ef⁻ = Eout[gf]; Ef⁺ = Eout[gf+1]
            Ef⁻,Ef⁺,isSkip = bounds_dispatch(interaction,Ef⁻,Ef⁺,Ei,gi,gf,type,Ui[δi],Ec,incoming_particle)
            if isSkip continue end
            ΔEf = Ef⁻ - Ef⁺
            
            # Integration over the energy group
            𝓕i = zeros(L+1)
            𝓕iₑ = 0
            for n in range(1,Np)

                # Outgoing particle energy group
                if (is_elastic) Ef = Ei else Ef = (u[n]*ΔEf + (Ef⁻+Ef⁺))/2 end

                # Compute Legendre angular flux moments
                Σsᵢ = ΔEf .* w[n]/2 .* dcs_dispatch(interaction,L,Ei,Ef,Z[i],scattered_particle,type,i,particles,Ein,Ef⁻,Ef⁺,δi,Ui[δi],Zi[δi],Ti[δi],ri[δi],Ec,incoming_particle) * nuclei_density(Z[i],ρ) * ωz[i]
                if is_dirac Σsᵢ /= ΔEf  end
                𝓕i .+= Σsᵢ
                𝓕iₑ += Σsᵢ[1] * (Ef+ΔQ)

            end
            𝓕[gf,:] .+= 𝓕i
            𝓕ₑ[gf] += 𝓕iₑ
        end
    end

elseif ((interaction isa Compton) && (interaction.model == "impulse_approximation")) # special case: random impulse approximation for Compton scattering
        
    interaction.is_subshells_dependant = true # ria needs subshell depedency
    type = "S" # only doing photon for now

    Nz = length(Z)
    for i in 1:Nz # for every element

        _,Zi,Ui,_,ri,_ = electron_subshells(Z[i],false) # electrons not free
        # Jio is necessary for compton profile. Scale from atomic units to inverse electron rest mass units 1/mₑc² by multiplying 1/α
        Jio = orbital_compton_profiles(Z[i]) .* 137.035999177

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
                dcs_moment_array .*= 2π # integrate over phi

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
end

return 𝓕, 𝓕ₑ
end
