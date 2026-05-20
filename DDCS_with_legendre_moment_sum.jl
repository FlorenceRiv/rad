using Plots
using QuadGK
using Radiant

include("Radiant.jl/src/cross_sections/electron_subshells.jl")
include("Radiant.jl/src/cross_sections/orbital_compton_profiles.jl")

Z = 13 # Atomic number (Al)
E_i = 10.0 # keV Initial photon energy
Ef_fixed = 0.98*E_i # Fixed outgoing energy 

_, Zi, Ui, _, _, _ = electron_subshells(Z, false)
Jio = orbital_compton_profiles(Z) .* 137.035999177
E0 = 1.0 
rₑ = 2.8179e-13 
mₑc² = 0.510999 # MeV

function keV_to_mec2(E_eV)
    return E_eV / 1e3 / mₑc²
end

E = keV_to_mec2(E_i)
Ef = keV_to_mec2(Ef_fixed)

function ddcs_ria(Ef, μ)
    Eq = sqrt(E^2 + Ef^2 - 2*E*Ef*μ)
    Ez = (E*Ef*(1-μ) - E0*(E - Ef)) / Eq
    Ec = E0*E / (E0 + E*(1-μ))
    F_Eμ = (E0*E*Ef / (Eq * Ec^2)) * inv(sqrt(1 + (Ez/E0)^2))
    summation = 0.
    for i in eachindex(Ui)
        if E - Ef > Ui[i] # Heaviside. otherwise = 0
            Ji = Jio[i] * (1+2*Jio[i]*abs(Ez)) * exp((1-(1+2*Jio[i]*abs(Ez))^2)/2)
            summation += Zi[i]*Ji ; end ; end
    DDCS = (rₑ^2 /2) * (Ec/E)^2 * (Ec/E +E/Ec + μ^2 - 1) * F_Eμ * summation * 1e24 # cm2 to barn
    return DDCS ; end # Omega and E_f derivative of CS

# Equation 15 Legendre Moments
function calculate_legendre_moment(l::Int, Ef::Float64)
    integrand = μ -> ddcs_ria(Ef, μ) * Radiant.legendre_polynomials(l, μ)
    moment, err = quadgk(integrand, -1.0, 1.0)
    return moment
end

# Equation 14: Reconstruct DDCS
function legendre_sum_ddcs(μ, L, legendre_moments)
    summation = 0.0
    for l in 0:L
        Pl = Radiant.legendre_polynomials(l, μ)
        summation += ((2*l + 1) / 2) * legendre_moments[l+1] * Pl
    end
    return summation
end

μ_vals = range(-1.0, 1.0, length=300)

ria_ddcs = [ddcs_ria(Ef, μ) for μ in μ_vals]

# Calc Legendre Moments from 0 to max_L
max_L = 14
moments = [calculate_legendre_moment(l, Ef) for l in 0:max_L]

plot(μ_vals, ria_ddcs, label="RIA DDCS", linewidth=3, color=:black)

for L in [0, 4, 8, 14]
    reconstructed_vals = [legendre_sum_ddcs(μ, L, moments) for μ in μ_vals]
    plot!(μ_vals, reconstructed_vals, label="Lsum DDCS(L=$L)", linewidth=2)
end

xlabel!("μ (cos(θ))")
ylabel!("σ_s (barn/sr)")
title!("DDCS (RIA and legendre moments sum) \n E_i = $E_i keV, E_f = $Ef_fixed keV")
