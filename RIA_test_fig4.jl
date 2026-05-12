"""
integrating DDCS(E,mu) over mu to get DCS(E)
fig 4 Brusa paper
"""

using Plots
using QuadGK

include("Radiant.jl/src/cross_sections/electron_subshells.jl")
include("Radiant.jl/src/cross_sections/orbital_compton_profiles.jl")

Z = 13 # atomic number (Al)
E = 50. # keV. initial photon energy

_, Zi, Ui, _, _, _ = electron_subshells(Z,false)
    """
    electron_subshells(Z::Int64,is_free::Bool=false)
    - `Z::Int64`: atomic number of the element.
    - `is_free::Bool` : are the electrons free (not bounded) ?
    
    Conversion of Ui in the function
        mₑc² = 0.510999  # MeV
        Ui ./= 1e6       # eV → MeV
        Ui ./= mₑc² # MeV → mₑc²

    return Nshells,Zi,Ui,Ti,ri,subshells
    """

Jio = orbital_compton_profiles(Z)
    """
    return J₀ 
    in ħ/(mₑe²) → atomic units
    """
# attempt at converting units (kills the curve)
"""
alpha = 0.0072973525643 # sommerfeld constant = e²/(ħc)
c = 2.99792458e8 # m/s
Jio .*= (c / alpha ) # ħ/(mₑe²) → 1/mₑc²
"""

# normalisation E
mₑc² = 0.510999 # MeV
E /= 1e3 # keV → MeV
E /= mₑc² # MeV → mₑc²

E_0 = 1. # (mₑc²) rest mass energy of electron 
r_e = 2.8179e-13 # (cm) electron radius


function calc_ddcs_ria(E_f, mu)

    E_q = sqrt(E^2 + E_f^2 - 2*E*E_f*mu)

    E_z = (E*E_f*(1-mu) - E_0*(E - E_f)) / E_q

    E_c = E_0*E / (E_0 + E*(1-mu)) 

    F_Emu = (E_0*E*E_f / E_q * E_c^2) * inv(sqrt(1 + (E_z/E_0)^2))

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

# quadgk integration (time inefficient)
function integrate_ddcs_over_mu(E_f)
    
    integral, error = quadgk(mu -> calc_ddcs_ria(E_f, mu), -1, 1)

    return integral

end

# array of E_f values
ϵ = 1e-6 # avoid division by 0 at E_f=E , mu=1 (see E_q)
E_f_vals = range(0.7*E, E-ϵ, length=200) # paper starts at 0.9?

# axis
xvals = E_f_vals ./ E
quadgk_yvals = [integrate_ddcs_over_mu(E_f) for E_f in E_f_vals]

Plots.plot(
    xvals, quadgk_yvals,
    xlabel = "E' / E",
    ylabel = "DCS (cm²)",
    label = "Al , E = 50keV"
)
