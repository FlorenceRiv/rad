"""
integrating DDCS(E,mu) over mu to get DCS(E)
comparing the different models
fig 4 Brusa paper
"""

using Radiant
using Plots
using QuadGK
using FastGaussQuadrature, LinearAlgebra


Z = 13 # atomic number (Al)
E_keV = 50. # keV. initial photon energy
x_start = 0.7

# number of nodes for Gauss quadrature
nodes = 50
# Gauss-Jacobi alpha, beta > -1
alpha = 1/3 
beta = -1/3

_, Zi, Ui, _, _, _ = Radiant.electron_subshells(Z,false)
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

Jio = Radiant.orbital_compton_profiles(Z)
    """
    return J₀ 
    in ħ/(mₑe²) → atomic units
    """
# divide by fine structure constant to match E_z units
Jio .*= 137.035999177 

# E unit conversion
E = E_keV
mₑc² = 0.510999 # MeV
E /= 1e3 # keV → MeV
E /= mₑc² # MeV → mₑc²

E_0 = 1. # (mₑc²) rest mass energy of electron 
r_e = 2.8179e-13 # (cm) electron radius


function calc_ddcs_ria(E_f, mu)

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
function quadgk_integrate_ddcs_over_mu(E_f)
    
    integral, error = quadgk(mu -> calc_ddcs_ria(E_f, mu), -1, 1)

    return integral
end


# array of E_f values
E_f_vals = range(x_start*E, E, length=200) # paper starts at 0.7E

# x axis
xvals = E_f_vals ./ E
quadgk_yvals = (E/Z)*2*pi*[quadgk_integrate_ddcs_over_mu(E_f) for E_f in E_f_vals]

# Bound Klein-Nishina and Waller-Hartree energies
E_min = E/(1 + 2E)
E_mask = E_f_vals .>= E_min
xvals_masked = E_f_vals[E_mask] ./ E
yvals_KN = [E * Radiant.klein_nishina(E, E_f) * 1e24 for E_f in E_f_vals[E_mask]] # cm2 to barn. K-N is per electron
yvals_WH = (E/Z)*[Radiant.waller_hartree(Z, E, E_f) for E_f in E_f_vals[E_mask]]*1e24

Plots.plot(
    xvals, quadgk_yvals,
    xlabel = "E' / E",
    ylabel = raw"$\frac{E}{Z}\frac{d\sigma}{dE}$ (barn)",
    label = "quadgk",
    title = "DCS(E') for Al, E = $E_keV keV, $nodes nodes",
    legend=:bottom
)
Plots.plot!(xvals_KN, yvals_KN, label = "Klein-Nishina")
Plots.plot!(xvals_masked, yvals_WH, label = "Waller-Hartree")
