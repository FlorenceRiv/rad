"""
making a plot of DDCS(E,mu)
fig 2 Brusa paper
"""

using Plots
include("Radiant.jl/src/cross_sections/electron_subshells.jl")
include("Radiant.jl/src/cross_sections/orbital_compton_profiles.jl")

Z = 13 # atomic number (Al)
E_keV = 10. # keV. initial photon energy

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

Jio .*= 137.035999177 # divide by fine structure constant to match E_z units


# E unit conversion
mₑc² = 0.510999 # MeV
E = E_keV
E /= 1e3 # keV → MeV
E /= mₑc² # MeV → mₑc²

E_0 = 1. # mₑc² rest mass energy of electron
theta1 = 60 # degrees
theta2 = 180
mu1 = cosd(theta1)
mu2 = cosd(theta2)
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

# array of E_f values
E_f_values = range(0.9*E, E, length=200)

# x axis
xvals = E_f_values ./ E
# y axis
yvals1 = (E/Z).*[calc_ddcs_ria(E_f, mu1) for E_f in E_f_values]
yvals2 = (E/Z).*[calc_ddcs_ria(E_f, mu2) for E_f in E_f_values]
Ec1 = E_c = E_0*E / (E_0 + E*(1-mu1))
Ec2 = E_c = E_0*E / (E_0 + E*(1-mu2))

Plots.plot(
    xvals, yvals1,

    xlabel = "E_f / E",
    ylabel = "(E/Z)*DDCS (barn/sr)",
    # Ji units: (mₑc²)⁻¹, E/Z units: mₑc², → (E/Z)*DDCS units: barn/sr

    label = "θ = 60°",
    title = "Al DDCS(E'), E = $E_keV keV"
)

Plots.plot!(xvals, yvals2, label = "θ = 180°")
Plots.vline!([Ec1/E], label = "E_c (θ=60°)", linestyle=:dash)
Plots.vline!([Ec2/E], label = "E_c (θ=180°)", linestyle=:dash)
