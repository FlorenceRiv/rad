"""
    ddcs_ria(E::Float64, Ef::Float64, Zi::Vector{Float64}, Ui::Vector{Float64}, Jio::Vector{Float64}, μ::Float64) 

returns σs (double derivative of the Compton cross-section) in barn/sr
based on the relativistic impulse approximation (RIA) model.

# Input Argument(s)

- `E::Float64`: incoming photon energy [in mₑc²].
- `Ef::Float64`: outgoing photon energy [in mₑc²].
- `Zi::Vector{Float64}`: vector of number of electrons in each subshell i.
- `Ui::Vector{Float64}`: binding energy per subshell [in mₑc²].
- `Jio::Vector{Float64}`: vector of constants necessary for the orbital compton profile of subshell i in mₑc²
- `μ::Float64`: cosine of scattering angle

# Output Argument(s)
- `σs::Float64` : double derivative of the Compton differential cross-section in barn/sr.

# Reference(s)
- Brusa et al. (1996), Fast sampling algorithm for the simulation of photon Compton
  scattering.

"""

E0 = 1. # mₑc² rest mass energy of electron
rₑ = 2.8179e-13 # (cm) electron radius


function ddcs_ria(E::Float64, Ef::Float64, Zi::Vector{Float64}, Ui::Vector{Float64}, Jio::Vector{Float64}, μ::Float64) 

# Ui and Jio vector constants are not calculated in this function to avoid unecessary repetition (and longer runtime)

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

    σs =
            (rₑ^2 /2) *
            (Ec/E)^2 *
            (Ec/E +E/Ec + μ^2 - 1) *
            F_Eμ *
            summation *
            1e24 # cm2 to barn

    return σs # Omega and E_f derivative of CS in barn/sr

end
