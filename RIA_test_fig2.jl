"""
making a plot of DDCS(E,mu)
fig 2 Brusa paper
"""

using Plots

include("Radiant.jl/src/cross_sections/electron_subshells.jl")
include("Radiant.jl/src/cross_sections/orbital_compton_profiles.jl")

Z = 13 # atomic number (Al)
E = 10. # keV. initial photon energy

_, Zi, Ui, _, _, _ = electron_subshells(Z,false)
    """  function info:

    electron_subshells(Z::Int64,is_free::Bool=false)
    - `Z::Int64`: atomic number of the element.
    - `is_free::Bool` : are the electrons free (not bounded) ?
    
    Conversions
        mₑc² = 0.510999  # MeV
        Ui ./= 1e6       # eV → MeV
        Ti ./= 1e6       # eV → MeV
        ri .*= 1e-11     # mÅ → cm
        Nshells = length(Zi)
        Ui ./= mₑc² # MeV → mₑc²
        Ti ./= mₑc² # MeV → mₑc² 

    return Nshells,Zi,Ui,Ti,ri,subshells
    """

Jio = orbital_compton_profiles(Z)
    """
    orbital_compton_profiles(Z::Int64)
    - `Z::Int64` : atomic number of the element.

    return J₀ 
    [in ħ/(mₑe²) → atomic units].
    """
"""
alpha = 0.0072973525643 # sommerfeld constant = e²/(ħc)
c = 2.99792458e8 # m/s
Jio .*= (c / alpha ) # ħ/(mₑe²) → 1/mₑc²
"""

# normalisation E
mₑc² = 0.510999 # MeV
E /= 1e3 # keV → MeV
E /= mₑc² # MeV → mₑc²

E_0 = 0.510999 # MeV rest mass energy of electron (mₑc²)
theta = 60 # degrees
mu = cosd(theta)
r_e = 2.8179e-13 # (cm) electron radius


function calc_ddcs_ria(E_f)

    E_q = sqrt(E^2 + E_f^2 - 2*E*E_f*mu)

    E_z = (E*E_f*(1-mu) - E_0*(E - E_f)) / E_q

    E_c = E_0*E / (E_0 + E*(1-mu)) 

    F_Emu = (E_0*E*E_f / E_q * E_c^2) * (1 + (E_z/E_0)^2)^(-1/2)

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
            summation

    return DDCS # Omega and E_f derivative of CS

end


E_f_values = range(0, E, length=200)

xvals = E_f_values ./ E

yvals = [calc_ddcs_ria(E_f) for E_f in E_f_values]


Plots.plot(
    xvals, yvals,

    xlabel = "E_f / E",
    ylabel = "DDCS (cm²)"
)
