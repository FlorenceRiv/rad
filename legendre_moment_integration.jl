"""
Legendre moments
integrating DDCS(E→E',μ)*P_l(μ) over μ [-1,1]
"""

using Plots
using QuadGK
using Radiant

include("Radiant.jl/src/cross_sections/electron_subshells.jl")
include("Radiant.jl/src/cross_sections/orbital_compton_profiles.jl")


Z = 13 # atomic number (Al)
E_i = 50. # keV. initial photon energy
Ef_start = 35. # keV (graph)


_, Zi, Ui, _, _, _ = electron_subshells(Z,false)


Jio = orbital_compton_profiles(Z) .* 137.035999177
# divide by fine structure constant to match E_z units


mₑc² = 0.510999 # MeV

function keV_to_mec2(E_eV)
    # keV → MeV , MeV → mₑc²
    return E_eV / 1e3 / mₑc²
end

function mec2_to_keV(E_mec2)
    # mₑc² → MeV, MeV → keV
    return E_mec2 * 1e3 * mₑc²
end

E = keV_to_mec2(E_i) # keV → mₑc²
Ef_start = keV_to_mec2(Ef_start) # eV → mₑc²

E0 = 1. # mₑc² rest mass energy of electron
rₑ = 2.8179e-13 # (cm) electron radius


function ddcs_ria(Ef, μ) # returns σs in barn/sr


    Eq = sqrt(E^2 + Ef^2 - 2*E*Ef*μ)

    Ez = (E*Ef*(1-μ) - E0*(E - Ef)) / Eq

    Ec = E0*E / (E0 + E*(1-μ))

    F_Eμ = (E0*E*Ef / (Eq * Ec^2)) * inv(sqrt(1 + (Ez/E0)^2))

    summation = 0.

    for i in eachindex(Ui)

        if E - Ef > Ui[i] # Heaviside. otherwise = 0

            Ji =
                Jio[i] *
                (1+2*Jio[i]*abs(Ez)) *
                exp((1-(1+2*Jio[i]*abs(Ez))^2)/2)

            summation += Zi[i]*Ji

        end
    end

    DDCS =
            (rₑ^2 /2) *
            (Ec/E)^2 *
            (Ec/E +E/Ec + μ^2 - 1) *
            F_Eμ *
            summation *
            1e24 # cm2 to barn

    return DDCS # Omega and E_f derivative of CS

end


"""
Radiant legendre_polynomials code:
    legendre_polynomials(l::Int64,x::Float64)


    Calculate the Legendre polynomials Pl(x).


    # Input Argument(s)
    - `l::Int64`: Legendre order.
    - `x::Float64`: evaluation points.


    # Output Argument(s)
    - `Pl::Vector{Float64}`: Legendre polynomial of order l evaluated at μ
"""

function calculate_legendre_moment(l::Int, Ef::Float64)
   
    integrand = μ -> ddcs_ria(Ef, μ) * Radiant.legendre_polynomials(l, μ)
   
    moment, err = quadgk(integrand, -1.0, 1.0)
   
    return moment
end


E_f_vals = range(Ef_start, E, length=200) # mₑc²

xvals = mec2_to_keV.(E_f_vals) # eV

# Initialize arrays
moments_l0 = Float64[]
moments_l1 = Float64[]
moments_l2 = Float64[]

for Ef_val in E_f_vals
    push!(moments_l0, calculate_legendre_moment(0, Ef_val))
    push!(moments_l1, calculate_legendre_moment(1, Ef_val))
    push!(moments_l2, calculate_legendre_moment(70, Ef_val))
end

plot(xvals, moments_l0, label="l = 0", linewidth=2)
plot!(xvals, moments_l1, label="l = 1", linewidth=2)
plot!(xvals, moments_l2, label="l = 2", linewidth=2)

xlabel!("Ef (eV)")
ylabel!("σ_s,l (barn/sr)")
title!("Legendre Moments (E → Ef) for Al, E = $E_i eV")

