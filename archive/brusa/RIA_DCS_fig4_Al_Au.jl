"""
integrating DDCS(E,mu) over mu to get DCS(E)
fig 4 Brusa paper
"""

using Plots
using QuadGK

include("Radiant.jl/src/cross_sections/electron_subshells.jl")
include("Radiant.jl/src/cross_sections/orbital_compton_profiles.jl")

Z_Al = 13 # atomic number (Al)
Z_Au = 79 # atomic number (Au)
E_keV = 50. # keV. initial photon energy
x_start = 0.7

function get_subshells_and_profiles(Z)

    _, Zi, Ui, _, _, _ = electron_subshells(Z,false)
    """
    electron_subshells(Z::Int64,is_free::Bool=false)
    - `Z::Int64`: atomic number of the element.
    - `is_free::Bool` : are the electrons free (not bounded) ?
    
    Conversion of Ui in the function
        mₑc² = 0.510999  # MeV
        Ui ./= 1e6       # eV → MeV
        Ui ./= mₑc²      # MeV → mₑc²

    return Nshells,Zi,Ui,Ti,ri,subshells
    """

    Jio = orbital_compton_profiles(Z)
    """
    return J₀ 
    in ħ/(mₑe²) → atomic units
    """
    # divide by fine structure constant to match E_z units
    Jio .*= 137.035999177 

    return Zi, Ui, Jio
end

Zi_Al, Ui_Al, Jio_Al = get_subshells_and_profiles(Z_Al)
Zi_Au, Ui_Au, Jio_Au = get_subshells_and_profiles(Z_Au)

# E unit conversion
E = E_keV
mₑc² = 0.510999 # MeV
E /= 1e3 # keV → MeV
E /= mₑc² # MeV → mₑc²

E_0 = 1. # (mₑc²) rest mass energy of electron 
r_e = 2.8179e-13 # (cm) electron radius


function calc_ddcs_ria(E_f, mu, Zi, Ui, Jio)

    E_q = sqrt(E^2 + E_f^2 - 2*E*E_f*mu)

    E_z = (E*E_f*(1-mu) - E_0*(E - E_f)) / E_q

    E_c = E_0*E / (E_0 + E*(1-mu)) 

    F_Emu = (E_0*E*E_f / (E_q * E_c^2)) * inv(sqrt(1 + (E_z/E_0)^2))

    summation = 0.

    for i in eachindex(Ui)

        if E - E_f > Ui[i] # Heaviside. otherwise = 0

            Ji = 
                Jio[i] * 
                (1+2*Jio[i]*abs(E_z)) * 
                exp((1-(1+2*Jio[i]*abs(E_z))^2)/2)

            summation += Zi[i]*Ji

        end
    end

    DDCS = 
            (r_e^2 /2) * 
            (E_c/E)^2 * 
            (E_c/E +E/E_c + mu^2 - 1) * 
            F_Emu * 
            summation *
            1e24 # cm2 to barn

    return DDCS # Omega and E_f derivative of CS
end

# quadgk integration (time inefficient)
function quadgk_integrate_ddcs_over_mu(E_f, Zi, Ui, Jio)
    
    integral, error = quadgk(mu -> calc_ddcs_ria(E_f, mu, Zi, Ui, Jio), -1, 1)

    return integral *2pi # integration over phi angle
end


# array of E_f values
E_f_vals = range(x_start*E, E, length=200) # paper starts at 0.7E

# x axis
xvals = E_f_vals ./ E
# y axis
Al_quadgk_yvals = (E/Z_Al)*[quadgk_integrate_ddcs_over_mu(E_f, Zi_Al, Ui_Al, Jio_Al) for E_f in E_f_vals]
Au_quadgk_yvals = (E/Z_Au)*[quadgk_integrate_ddcs_over_mu(E_f, Zi_Au, Ui_Au, Jio_Au) for E_f in E_f_vals]

Plots.plot(
    xvals, Al_quadgk_yvals,
    xlabel = "E' / E",
    ylabel = raw"$\frac{E}{Z}\frac{d\sigma}{dE}$ (barn)",
    label = "Al",
    title = "DCS(E'), E = $E_keV keV",
    legend=:topleft
)
Plots.plot!(
    xvals, Au_quadgk_yvals,
    label = "Au",
)
