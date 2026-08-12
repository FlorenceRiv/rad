"""
compare heaviside vs subshells runtime
integrating DDCS(E,mu) over mu to get DCS(E)
fig 4 Brusa paper
"""

using Radiant
using Plots
using QuadGK


Z = 13 # atomic number (Al)
E_keV = 50. # keV. initial photon energy
x_start = 0.7

_, Zi, Ui, _, _, _ = Radiant.electron_subshells(Z,false)

Jio = Radiant.orbital_compton_profiles(Z)
# divide by fine structure constant to match E_z units
Jio .*= 137.035999177 

# E unit conversion
E = E_keV
mₑc² = 0.510999 # MeV
E /= 1e3 # keV → MeV
E /= mₑc² # MeV → mₑc²

E_0 = 1. # (mₑc²) rest mass energy of electron 
r_e = 2.8179e-13 # (cm) electron radius


function ddcs_ria(E_f, mu)

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

    return DDCS # Omega and E_f derivative of CS in barn/sr
end


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
            Zi_s*Ji *
            1e24 # cm2 to barn

    return σs # Omega and E_f derivative of CS in barn/sr
end


function quadgk_integrate_ddcs_over_mu(E_f)
    integral, error = quadgk(mu -> ddcs_ria(E_f, mu), -1.0, 1.0)
    return integral
end

function quadgk_integrate_subshells_over_mu(E_f)
    total_integral = 0.0
    for i in eachindex(Ui)
        integral, _ = quadgk(mu -> ddcs_ria_subshell(E, E_f, Zi[i], Ui[i], Jio[i], mu), -1.0, 1.0)
        total_integral += integral
    end
    return total_integral
end

# array of E_f values
E_f_vals = range(x_start*E, E, length=200)

# x axis
xvals = E_f_vals ./ E

quadgk_integrate_ddcs_over_mu(E_f_vals[1])
quadgk_integrate_subshells_over_mu(E_f_vals[1])

y_multiplier = (E/Z) * 2 * pi
quadgk_yvals = y_multiplier .* [quadgk_integrate_ddcs_over_mu(E_f) for E_f in E_f_vals]
subshell_yvals = y_multiplier .* [quadgk_integrate_subshells_over_mu(E_f) for E_f in E_f_vals]

t1 = @elapsed [quadgk_integrate_ddcs_over_mu(E_f) for E_f in E_f_vals]
t2 = @elapsed [quadgk_integrate_subshells_over_mu(E_f) for E_f in E_f_vals]

println("Integration times (seconds):")
println("heaviside: ", t1)
println("subshells: ", t2)
println("ratio: ", t1/t2)

Plots.plot(
    xvals, quadgk_yvals,
    xlabel = "E' / E",
    ylabel = raw"$\frac{E}{Z}\frac{d\sigma}{dE}$ (barn)",
    label = "ddcs_ria_heaviside",
    title = "DCS(E') for Al, E = $E_keV keV", # Removed undefined $nodes variable
    legend=:topleft
)

Plots.plot!(xvals, subshell_yvals, label = "ddcs_ria_subshell", linestyle=:dash)
