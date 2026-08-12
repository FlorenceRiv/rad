

"""
    feed(Z::Vector{Int64},atz::Vector{Float64},L::Int64,Ei::Float64,
    Eout::Vector{Float64},Ng::Int64,interaction::Interaction,gi::Int64,Ngi::Int64,
    particles::Vector{Particle},type::String,incoming_particle::Particle,
    scattered_particle::Particle,Ein::Vector{Float64},Ec::Float64,
    is_elastic::Bool,is_subshells::Bool)

Calculate the feed function 𝓕 (normalized probability of scattering from Ei into each
group gf) for each Legendre moment up to order L. Also calculate the energy weighted
feed function 𝓕ₑ for energy-deposition cross section.

# Input Argument(s)
- `Z::Vector{Int64}` : atomic number of the element(s) composing the material.
- `atz::Vector{Float64}` : atomic percent of the element(s) composing the material.
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
- `is_elastic::Bool` : boolean indicating if the scattering is elastic.
- `is_subshells::Bool` : boolean indicating if the cross-sections are subshells dependant.

# Output Argument(s)
- `𝓕::Array{Float64}` : feed function (per unit nuclei density).
- `𝓕ₑ::Vector{Float64}` : energy weighted feed function (per unit nuclei density).

# Reference(s)
- MacFarlane et al. (2021) : The NJOY Nuclear Data Processing System, Version 2012.

"""
function feed(Z::Vector{Int64},atz::Vector{Float64},L::Int64,Ei::Float64,Eout::Vector{Float64},Ng::Int64,interaction::Interaction,gi::Int64,Ngi::Int64,particles::Vector{Particle},type::String,incoming_particle::Particle,scattered_particle::Particle,Ein::Vector{Float64},Ec::Float64,is_elastic::Bool,is_subshells::Bool)

#----
# Initialization
#----
𝓕 = zeros(Ng+1,L+1)
𝓕ₑ = zeros(Ng+1)

if ! ((interaction isa Compton) && (interaction.model == "impulse_approximation") && (type == "S")) # check for special case of compton ria, otherwise do the regular feed calculation

#   DOING ELECTRON COMPTON TOO FOR DEBUG. REMOVE THIS LATER: && (type == "S")

    ΔQ = get_mass_energy_variation(interaction,type,true)

    # Outgoing particle energy spectrum
    is_dirac, Np, q_type = out_distribution_dispatch(interaction,type)
    if is_dirac Np = 1; u = [0]; w = [2] else u,w = quadrature(Np,q_type) end

    #----
    # Feed function over all groups and under the cutoff energy
    #----

    # Heavy inelastic S (same scattered particle)
    is_heavy_inelastic_S = (is_proton(incoming_particle) || is_alpha(incoming_particle)) && (incoming_particle == scattered_particle)

    # Loop over the compound elements
    Nz = length(Z)
    
    for i in range(1,Nz)

        # Loop over subshells and outgoing groups
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

        # For heavy particles compute cache once and reuse for all analytic integrals
        analytic_A = 0.0
        M₁ = 0.0
        cache = nothing
        if is_heavy_inelastic_S
            cache = HeavyInelasticCache(Zi[δi], Ei, incoming_particle)
            analytic_A = integrate_A_over_W2_per_subshell(cache, Ef⁻, Ef⁺) * atz[i]
            M₁ = feed_first_moment_heavy_particle(cache, Ef⁻, Ef⁺) * atz[i]
        end

        # Use quadrature integration for all particles
            for n in range(1,Np)
                # Outgoing particle energy group
                if (is_elastic) Ef = Ei else Ef = (u[n]*ΔEf + (Ef⁻+Ef⁺))/2 end

                # Compute Legendre angular flux moments
                Σsᵢ = ΔEf .* w[n]/2 .* dcs_dispatch(interaction,L,Ei,Ef,Z[i],scattered_particle,type,i,particles,Ein,Ef⁻,Ef⁺,δi,Ui[δi],Zi[δi],Ti[δi],ri[δi],Ec,incoming_particle) * atz[i]
                if is_dirac Σsᵢ /= ΔEf end
                𝓕i .+= Σsᵢ
                𝓕iₑ += Σsᵢ[1] * (Ef + ΔQ)
        end

        # Add analytic singular contribution
        if is_heavy_inelastic_S
            𝓕i .+= analytic_A .* ones(L+1)
            for l in range(0,L)
                lead_log = integrate_leading_1overW_per_subshell(cache, Ef⁻, Ef⁺, l) * atz[i]
                𝓕i[l+1] += lead_log
            end
            σ_analytic = feed_analytical_heavy_particle(cache, Ef⁻, Ef⁺) * atz[i]
            𝓕iₑ = Ei * σ_analytic - M₁

            if ~isapprox(𝓕i[1], σ_analytic; rtol=1e-3, atol=1e-12)
                rel = abs(𝓕i[1] - σ_analytic) / max(abs(σ_analytic), 1e-300)
                print("FEED_WARN: Z=$(Z[i]), δi=$(δi), gf=$(gf), σ_analytic=$(σ_analytic), numeric_l0=$(𝓕i[1]), rel_diff=$(rel)\n")
                𝓕i[1] = σ_analytic
            end
            end
            𝓕[gf,:] .+= 𝓕i
            𝓕ₑ[gf] += 𝓕iₑ
        end
    end

elseif ((interaction isa Compton) && (interaction.model == "impulse_approximation") && (type == "S")) 
    # special case: random impulse approximation for Compton scattering
    # only type == 'S' because electron production for RIA not yet implemented
        # 'P' is handled by the regular feed function
        
    # for ria, need subshell depedency: "interaction.is_subshells_dependant = true".
    # but not set here to not modify the interaction object outside of feed()

    Nz = length(Z)
    for i in 1:Nz

        _, Zi, Ui, _, _, _ = electron_subshells(Z[i], false)
        if Ei < minimum(Ui); continue; end

        Jio = orbital_compton_profiles(Z[i]) .* 137.035999177

        for s in range(1,length(Ui)) # go over subshells

            if Ei <= Ui[s]; continue; end # skip if incoming energy smaller than ionization energy

            for gf in 1:Ng

                Ef⁻ = Eout[gf]; Ef⁺ = Eout[gf+1]
                Ef⁻, Ef⁺, isSkip = bounds_dispatch(interaction, Ef⁻, Ef⁺, Ei, gi, gf, type, Ui[s], Ec, incoming_particle) # crops Ef- for RIA
                if isSkip; continue; end
                if Ef⁻ <= Ef⁺; continue; end

                group_moments, _ = quadgk(Ef⁺, Ef⁻) do Ef

                    # integrate legendre moment (ddcs over mu) for each l :
                    dcs_moment_array, _ = quadgk(-1.0, 1.0) do μ
                    return [ddcs_ria_subshell(Ei, Ef, Zi[s], Ui[s], Jio[s], μ) * Radiant.legendre_polynomials(l, μ) for l in 0:L]
                    end

                    dcs_moment_array .*= 2π # integrate over phi

                    # energy-weighted L=0 moment for 𝓕ₑ. ΔQ=0
                    energy_weighted = dcs_moment_array[1] * Ef
                
                    # push the energy-weighted value to the end of the vector
                    # allows quadgk to integrate 𝓕 and 𝓕ₑ in a single pass
                    return push!(dcs_moment_array, energy_weighted)
                end

                # Legendre moments (indices 1 to L+1) and scale
                𝓕[gf, :] .+= group_moments[1:L+1] .* atz[i]

                # energy-weighted value (index L+2) and scale
                𝓕ₑ[gf] += group_moments[L+2] * atz[i]

            end
        end
    end

elseif ((interaction isa Compton) && (interaction.model == "impulse_approximation") && (type == "P")) # electron compton RIA

    Nz = length(Z)
    for i in 1:Nz

        _, Zi, Ui, _, _, _ = electron_subshells(Z[i], false)
        if Ei < minimum(Ui); continue; end

        Jio = orbital_compton_profiles(Z[i]) .* 137.035999177

        for s in range(1,length(Ui)) # go over subshells

            if Ei <= Ui[s]; continue; end # skip if incoming energy smaller than ionization energy

            for gf in 1:Ng 

                # final energies of electrons
                Ef⁻ = Eout[gf]; Ef⁺ = Eout[gf+1]

                Ef⁻, Ef⁺, isSkip = bounds_dispatch(interaction, Ef⁻, Ef⁺, Ei, gi, gf, type, Ui[s], Ec, incoming_particle) # crops Ef⁻ at min(Ef⁻, Ei-Ui) for RIA
                if isSkip; continue; end
                if Ef⁻ <= Ef⁺; continue; end

                group_moments, _ = quadgk(Ef⁺, Ef⁻) do Ef

                    Efγ = Ei - Ui[s] -Ef

                    # integrate legendre moment (ddcs over mu) for each l :
                    # electron cosine range: 0 to 1
                    dcs_moment_array, _ = quadgk(0.0, 1.0) do μ
                    return [
                                ddcs_ria_subshell(Ei, Efγ, Zi[s], Ui[s], Jio[s], convert_mu_electron_to_mu_photon(μ, Ei)) *
                                jacobian_for_mu_electron_to_mu_photon_conversion(μ, Ei)* 
                                Radiant.legendre_polynomials(l, μ) 

                            for l in 0:L]
                    end

                    dcs_moment_array .*= 2π # integrate over phi

                    # energy-weighted L=0 moment for 𝓕ₑ. ΔQ=0
                    energy_weighted = dcs_moment_array[1] * Ef ################ 
                
                    # push the energy-weighted value to the end of the vector
                    # allows quadgk to integrate 𝓕 and 𝓕ₑ in a single pass
                    return push!(dcs_moment_array, energy_weighted)
                end

                # Legendre moments (indices 1 to L+1) and scale
                𝓕[gf, :] .+= group_moments[1:L+1] .* atz[i]

                # energy-weighted value (index L+2) and scale
                𝓕ₑ[gf] += group_moments[L+2] * atz[i]

            end
        end
    end

end
return 𝓕, 𝓕ₑ
end

"""
    feed_elastic_scattering(Z::Vector{Int64},atz::Vector{Float64},L::Int64,Ei::Float64,
    Eout::Vector{Float64},Ng::Int64,interaction::Interaction,gi::Int64,Ngi::Int64,
    particles::Vector{Particle},type::String,incoming_particle::Particle,
    scattered_particle::Particle,Ein::Vector{Float64},Ec::Float64,is_elastic::Bool,
    is_subshells::Bool,A::Vector{Vector{Int64}},
    atpercentA::Vector{Vector{Float64}})

Calculate the elastic-scattering feed function 𝓕 from incident energy `Ei` into each
outgoing energy group and Legendre moment up to order `L`. Also calculate the
energy-weighted feed function 𝓕ₑ for energy-deposition cross sections, including isotope
fractions when isotope-resolved data are provided.

# Input Argument(s)
- `Z::Vector{Int64}` : atomic number of the element(s) composing the material.
- `atz::Vector{Float64}` : atomic percent of the element(s) composing the material.
- `L::Int64` : Legendre truncation order.
- `Ei::Float64` : energy of the incoming particle [in mₑc²].
- `Eout::Vector{Float64}` : energy group boundaries [in mₑc²].
- `Ng::Int64` : number of groups.
- `interaction::Interaction` : interaction information.
- `gi::Int64` : incoming particle group index.
- `Ngi::Int64` : number of groups for the incoming particle.
- `particles::Vector{Particle}` : list of particles involved in the interaction.
- `type::String` : type of interaction (`"S"` for scattering or `"P"` for production).
- `incoming_particle::Particle` : incoming particle.
- `scattered_particle::Particle` : scattered particle.
- `Ein::Vector{Float64}` : energy group boundaries corresponding to the incoming
  particle [in mₑc²].
- `Ec::Float64` : cutoff energy between soft and catastrophic interaction.
- `is_elastic::Bool` : boolean indicating if the outgoing particle energy is equal to `Ei`.
- `is_subshells::Bool` : boolean indicating if subshell-dependent cross sections are used.
- `A::Vector{Vector{Int64}}` : isotope mass numbers per element.
- `atpercentA::Vector{Vector{Float64}}` : isotope atomic fractions per element.

# Output Argument(s)
- `𝓕::Array{Float64}` : feed function.
- `𝓕ₑ::Vector{Float64}` : energy-weighted feed function.

# Reference(s)
- MacFarlane et al. (2021) : The NJOY Nuclear Data Processing System, Version 2012.

"""
function feed_elastic_scattering(Z::Vector{Int64},atz::Vector{Float64},L::Int64,Ei::Float64,Eout::Vector{Float64},Ng::Int64,interaction::Interaction,gi::Int64,Ngi::Int64,particles::Vector{Particle},type::String,incoming_particle::Particle,scattered_particle::Particle,Ein::Vector{Float64},Ec::Float64,is_elastic::Bool,is_subshells::Bool,A::Vector{Vector{Int64}},atpercentA::Vector{Vector{Float64}})

#----
# Initialization
#----
𝓕 = zeros(Ng+1,L+1)
𝓕ₑ = zeros(Ng+1)

#----
# Feed function over all groups and under the cutoff energy
#----
a
# Loop over the compound elements
Nz = length(Z)
for i in range(1,Nz)
    # Loop over isotopes
    for (Ai, atai) in zip(A[i], atpercentA[i])
        if type == "P" && !(Z[i] == 1 && Ai == 1)
            continue
        end
        δi = 0
        Ui = 0.0
        Zi = Z[i]
        Ti = 0.0
        ri = 0.0
        for gf in range(1,Ng)

            # Final energy group
            Ef⁻ = Eout[gf]; Ef⁺ = Eout[gf+1]
            M_target = get_mass(Z[i], Ai)
            Ef⁻,Ef⁺,isSkip = bounds_dispatch(interaction,Ef⁻,Ef⁺,Ei,gi,gf,type,Ui,Ec,incoming_particle,M_target)
            if isSkip continue end
            ΔEf = Ef⁻ - Ef⁺

            # Integration over the energy group
            𝓕i = zeros(L+1)
            𝓕iₑ = 0
            Ef = is_elastic ? Ei : (Ef⁻ + Ef⁺) / 2
            Σsᵢ = dcs_dispatch(interaction,L,Ei,Ef,Z[i],scattered_particle,type,i,particles,Ein,Ef⁻,Ef⁺,δi,Ui,Zi,Ti,ri,Ec,incoming_particle,Ai) * atz[i] * atai
            𝓕i .+= Σsᵢ
            𝓕iₑ += Σsᵢ[1] * Ef
            𝓕[gf,:] .+= 𝓕i
            𝓕ₑ[gf] += 𝓕iₑ
        end
    end
end
return 𝓕, 𝓕ₑ
end
