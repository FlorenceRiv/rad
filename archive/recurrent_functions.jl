"""
Functions I use frequently

    - inverse_mₑc²_units_orbital_compton_profiles(Z)

    - keV_to_mec2(E_keV)
    - eV_to_mec2(E_eV)
    - mec2_to_keV(E_mec2)
    - mec2_to_eV(E_mec2)

    - ddcs_ria(E, Ef, Ui, Zi, Jio, μ)

    - legendre_moment(l, E, Ef, Zi, Ui, Jio)

"""

using QuadGK
using Radiant

include("Radiant.jl/src/cross_sections/orbital_compton_profiles.jl")

function inverse_mₑc²_units_orbital_compton_profiles(Z)

    # Scale from atomic units to inverse electron rest mass units 1/mₑc² by multiplying 1/α
    return orbital_compton_profiles(Z) .* 137.035999177

end

mₑc² = 0.510999 # MeV

function keV_to_mec2(E_keV)
    # keV → MeV , MeV → mₑc²
    return E_keV / 1e3 / mₑc²
end

function eV_to_mec2(E_eV)
    # eV → MeV , MeV → mₑc²
    return E_eV / 1e6 / mₑc²
end

function mec2_to_keV(E_mec2)
    # mₑc² → MeV, MeV → keV
    return E_mec2 * 1e3 * mₑc²
end

function mec2_to_eV(E_mec2)
    # mₑc² → MeV, MeV → eV
    return E_mec2 * 1e6 * mₑc²
end


E0 = 1. # mₑc² rest mass energy of electron
rₑ = 2.8179e-13 # (cm) electron radius


function ddcs_ria(E, Ef, Zi, Ui, Jio, μ) # returns σs in barn/sr

# Ui and Jio vector constants are not calculated in this function to avoid unecessary repetition (and longer runtime)

    # E: initial photon energy in mₑc²
    # Ef: outgoing photon energy in mₑc²
    # Zi: vector of number of electrons in each subshell i
    # Ui: binding energy of subshell i in mₑc²
    # Jio: part of orbital compton profile of subshell i in mₑc²
    # μ: cosine of scattering angle

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
    - Pl: Legendre polynomial of order l evaluated at μ
"""

function legendre_moment(l, E, Ef, Zi, Ui, Jio)

    # l: Legendre order.
    # E: initial photon energy in mₑc²
    # Ef: outgoing photon energy in mₑc²
    # Zi: vector of number of electrons in each subshell i
    # Ui: binding energy of subshell i in mₑc²
    # Jio: part of orbital compton profile of subshell i in mₑc²

    integrand = μ -> ddcs_ria(E, Ef, Zi, Ui, Jio, μ) * Radiant.legendre_polynomials(l, μ)

    moment, err = quadgk(integrand, -1.0, 1.0)
   
    return moment
end


