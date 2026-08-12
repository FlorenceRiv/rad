
"""
    ddcs_ria_subshell(E::Float64, Ef::Float64, Zi_s::Float64, Ui_s::Float64, Jio_s::Float64, μ::Float64) 

returns single-subshell σs (double derivative of the Compton cross-section) in cm2
based on the relativistic impulse approximation (RIA) model.

The Heaviside condition E - Ef > Ui_s is enforced by the caller via the integration bounds.

# Input Argument(s)

- `E::Float64`: incoming photon energy [in mₑc²].
- `Ef::Float64`: outgoing photon energy [in mₑc²].
- `Zi_s::Float64`: number of electrons in subshell s.
- `Ui_s::Float64`: subshell binding energy [in mₑc²].
- `Jio_s::Float64`: constant necessary for the orbital compton profile of subshell s in mₑc²
- `μ::Float64`: cosine of scattering angle

# Output Argument(s)
- `σs::Float64` : double derivative of the Compton differential cross-section in cm2 for specific subshell.

# Reference(s)
- Brusa et al. (1996), Fast sampling algorithm for the simulation of photon Compton
  scattering.

"""

function ddcs_ria_subshell(E::Float64, Ef::Float64, Zi_s::Float64, Ui_s::Float64, Jio_s::Float64, μ::Float64)

    if E - Ef <= Ui_s return 0 end

    E0 = 1. # mₑc² rest mass energy of electron
    rₑ = 2.8179e-13 # (cm) electron radius


    Eq = sqrt(E^2 + Ef^2 - 2*E*Ef*μ)

    Ez = (E*Ef*(1-μ) - E0*(E - Ef)) / Eq

    Ec = E0*E / (E0 + E*(1-μ))

    F_Eμ = (E0*E*Ef / (Eq * Ec^2)) * inv(sqrt(1 + (Ez/E0)^2))

    Ji =
                Jio_s *
                (1+2*Jio_s*abs(Ez)) *
                exp((1-(1+2*Jio_s*abs(Ez))^2)/2)

    σs =
            (rₑ^2 /2) *
            (Ec/E)^2 *
            (Ec/E +E/Ec + μ^2 - 1) *
            F_Eμ *
            Zi_s*Ji

    return σs # Omega and E_f derivative of CS in cm2
end


"""
    ddcs_ria_heaviside(E::Float64, Ef::Float64, Zi::Vector{Float64}, Ui::Vector{Float64}, Jio::Vector{Float64}, μ::Float64) 

returns σs (double derivative of the Compton cross-section) in barn/sr
based on the relativistic impulse approximation (RIA) model.

# Input Argument(s)

- `E::Float64`: incoming photon energy [in mₑc²].
- `Ef::Float64`: outgoing photon energy [in mₑc²].
- `Zi::Vector{Float64}`: vector of number of electrons in each subshell i.
- `Ui::Vector{Float64}`: vector of binding energy per subshell [in mₑc²].
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


function ddcs_ria_heaviside(E::Float64, Ef::Float64, Zi::Vector{Float64}, Ui::Vector{Float64}, Jio::Vector{Float64}, μ::Float64) 

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
            summation #*
            #1e24 # cm2 to barn # STAY IN CM2   

    return σs # Omega and E_f derivative of CS in barn/sr

end



"""
    convert_mu_electron_to_mu_photon(μ_electron::Float64, Ei_photon::Float64)

Returns the cosine of the photon scattering angle based on the cosine of the electron scattering angle and the incoming photon energy.

# Input Argument(s)
- `μ_electron::Float64`: cosine of the electron scattering angle.
- `Ei_photon::Float64`: incoming photon energy [in mₑc²].

# Output Argument(s)
- `μ_photon::Float64`: cosine of the photon scattering angle.

# Reference(s)
- wikipedia Compton electron scattering angle formula + algebraic manipulations 
"""

function convert_mu_electron_to_mu_photon(μ_electron::Float64, Ei_photon::Float64)

    μ_photon = ( (1-μ_electron^2)*(1+Ei_photon)^2 - μ_electron^2 ) / ( (1-μ_electron^2)*(1+Ei_photon)^2 + μ_electron^2 )
    
    return μ_photon

end    

"""
    jacobian_for_mu_electron_to_mu_photon_conversion(μ_electron::Float64, Ei_photon::Float64)

Returns the jacobian |dμ_photon/dμ_electron| for the electron Legendre moment integral calculation.

# Input Argument(s)
- `μ_electron::Float64`: cosine of the electron scattering angle.
- `Ei_photon::Float64`: incoming photon energy [in mₑc²].

# Output Argument(s)
- `μ_transition_jacobian::Float64`: the absolute value of the derivative of the transformation.

# Reference(s)
- absolute derivative of μ_photon(μ_electron)
"""

function jacobian_for_mu_electron_to_mu_photon_conversion(μ_electron::Float64, Ei_photon::Float64)

    μ_transition_jacobian = 4*μ_electron*(1+Ei_photon)^2 / ( (1-μ_electron^2)*(1+Ei_photon)^2 + μ_electron^2 )^2

    return μ_transition_jacobian

end





"""
    μ_electron_breakpoint(Efγ::Float64, Ei::Float64)

Returns the μ_electron value corresponding to the free-electron Compton locus
(the peak of the Doppler-broadened ridge), by inverting
`convert_mu_electron_to_mu_photon`. Returns `nothing` if the ridge peak in
μ_photon-space falls outside [-1, 1] (no breakpoint needed).

# Input Argument(s)
- `Efγ::Float64`: photon-equivalent outgoing energy [in mₑc²], i.e. Ei - Ui_s - Ef.
- `Ei::Float64`: incoming photon energy [in mₑc²].

# Output Argument(s)
- `μₑ_peak::Union{Float64,Nothing}`: μ_electron breakpoint, or `nothing` if out of range.
"""
function μ_electron_breakpoint(Efγ::Float64, Ei::Float64)

    # same ridge location as the photon-branch feed(), just in Efγ
    μ_peak_photon = 1.0 - 1.0/Efγ + 1.0/Ei

    if !(-1.0 < μ_peak_photon < 1.0)
        return nothing
    end

    A = (1.0 + Ei)^2
    μₑ_peak_sq = A * (1.0 - μ_peak_photon) / ( (A + 1.0) - μ_peak_photon * (A - 1.0) )

    # clamp for floating-point round-off right at the domain edge
    return sqrt(max(μₑ_peak_sq, 0.0))
end


function radiant_test()
    println("this better work")
end