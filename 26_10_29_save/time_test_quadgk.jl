using Radiant
using Plots
using QuadGK


Z = 79 # atomic number (Au)
E_keV = 100. # keV. initial photon energy
x_start = 0.7

_, Zi, Ui, _, _, _ = Radiant.electron_subshells(Z,false)

Jio = Radiant.orbital_compton_profiles(Z)
Jio .*= 137.035999177 

E = E_keV
mₑc² = 0.510999 # MeV
E /= 1e3
E /= mₑc²

E_0 = 1.
r_e = 2.8179e-13


function ddcs_ria_subshell(E::Float64, Ef::Float64, Zi_s::Float64, Ui_s::Float64, Jio_s::Float64, μ::Float64)

    if E - Ef <= Ui_s return 0 end

    E0 = 1.
    rₑ = 2.8179e-13

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
            1e24

    return σs
end


# --- variant 1: no breakpoint
function quadgk_integrate_subshells_nobreak(E_f)
    total_integral = 0.0
    for i in eachindex(Ui)
        integral, _ = quadgk(mu -> ddcs_ria_subshell(E, E_f, Zi[i], Ui[i], Jio[i], mu), -1.0, 1.0)
        total_integral += integral
    end
    return total_integral
end

# --- variant 2: with μ_peak breakpoint ---
function quadgk_integrate_subshells_breakpoint(E_f)
    μ_peak = 1.0 - 1.0/E_f + 1.0/E
    μ_pts = -1.0 < μ_peak < 1.0 ? (-1.0, μ_peak, 1.0) : (-1.0, 1.0)

    total_integral = 0.0
    for i in eachindex(Ui)
        integral, _ = quadgk(mu -> ddcs_ria_subshell(E, E_f, Zi[i], Ui[i], Jio[i], mu), μ_pts...)
        total_integral += integral
    end
    return total_integral
end


E_f_vals = range(x_start*E, E, length=200)
xvals = E_f_vals ./ E

# warm up / compile
quadgk_integrate_subshells_nobreak(E_f_vals[1])
quadgk_integrate_subshells_breakpoint(E_f_vals[1])

y_multiplier = (E/Z) * 2 * pi
nobreak_yvals    = y_multiplier .* [quadgk_integrate_subshells_nobreak(E_f)    for E_f in E_f_vals]
breakpoint_yvals = y_multiplier .* [quadgk_integrate_subshells_breakpoint(E_f) for E_f in E_f_vals]

t1 = @elapsed [quadgk_integrate_subshells_nobreak(E_f)    for E_f in E_f_vals]
t2 = @elapsed [quadgk_integrate_subshells_breakpoint(E_f) for E_f in E_f_vals]

println("Integration times (seconds):")
println("no breakpoint:   ", t1)
println("with breakpoint: ", t2)
println("ratio (nobreak/breakpoint): ", t1/t2)

# check whether the breakpoint also changes the *value*, not just the speed
max_reldiff = maximum(abs.(nobreak_yvals .- breakpoint_yvals) ./ abs.(breakpoint_yvals))
println("max relative difference in integral value: ", max_reldiff)

Plots.plot(
    xvals, nobreak_yvals,
    xlabel = "E' / E",
    ylabel = raw"$\frac{E}{Z}\frac{d\sigma}{dE}$ (barn)",
    label = "no breakpoint",
    title = "DCS(E') for Au, E = $E_keV keV",
    legend=:topleft
)


Plots.plot!(xvals, breakpoint_yvals, label = "with μ_peak breakpoint", linestyle=:dash)