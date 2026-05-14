"""
integrating DDCS(E,mu) over mu to get DCS(E)
fig 4 Brusa paper
"""

using Plots
using QuadGK
using FastGaussQuadrature, LinearAlgebra

include("Radiant.jl/src/cross_sections/electron_subshells.jl")
include("Radiant.jl/src/cross_sections/orbital_compton_profiles.jl")

Z = 13 # atomic number (Al)
E_keV = 50. # keV. initial photon energy

# number of nodes for Gauss quadrature
nodes = 5
# Gauss-Jacobi alpha, beta > -1
alpha = 1/3 
beta = -1/3

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
function quadgk_integrate_ddcs_over_mu(E_f)
    
    integral, error = quadgk(mu -> calc_ddcs_ria(E_f, mu), -1, 1)

    return integral
end

# ===== the next functions are from FastGaussQuadrature.jl =====

# Gauss-Legendre quadrature integration (faster but is it precise for this case?)
function GaussLegendre_integrate_ddcs_over_mu(E_f)
    
    mu, weights = gausslegendre(nodes)

    integral = dot(weights, calc_ddcs_ria.(E_f, mu))

    return integral
end

# Gauss-Chebyshev1
function GaussChebyshev1_integrate_ddcs_over_mu(E_f)
    
    mu, weights = gausschebyshevt(nodes)

    integral = dot(weights, calc_ddcs_ria.(E_f, mu).*sqrt.(1 .- mu.^2)) 
    # divide function by weight w1(x) = 1/sqrt(1-x^2)

    return integral
end

# Gauss-Chebyshev2
function GaussChebyshev2_integrate_ddcs_over_mu(E_f)
    
    mu, weights = gausschebyshevu(nodes)

    integral = dot(weights, calc_ddcs_ria.(E_f, mu)./sqrt.(1 .- mu.^2))
    # divide function by weight w2(x) = sqrt(1-x^2)

    return integral
end

# Gauss-Chebyshev3
function GaussChebyshev3_integrate_ddcs_over_mu(E_f)
    
    mu, weights = gausschebyshevv(nodes)

    integral = dot(weights, calc_ddcs_ria.(E_f, mu).*sqrt.((1 .- mu)./(1 .+ mu)))
    # divide function by weight w3(x) = sqrt((1 + mu)/(1 - mu))

    return integral
end

# Gauss-Chebyshev4
function GaussChebyshev4_integrate_ddcs_over_mu(E_f)
    
    mu, weights = gausschebyshevw(nodes)

    integral = dot(weights, calc_ddcs_ria.(E_f, mu).*sqrt.((1 .+ mu)./(1 .- mu)))
    # divide function by weight w4(x) = sqrt((1 - mu)/(1 + mu))

    return integral
end

# Gauss-Jacobi
function GaussJacobi_integrate_ddcs_over_mu(E_f, alpha, beta)
    
    mu, weights = gaussjacobi(nodes, alpha, beta)

    integral = dot(weights, calc_ddcs_ria.(E_f, mu)./((1 .- mu).^alpha.*(1 .+ mu).^beta))
    # divide function by weight w(x) = (1-x)^alpha * (1+x)^beta

    return integral
end

# Gauss-Radau
function GaussRadau_integrate_ddcs_over_mu(E_f)

    mu, weights = gaussradau(nodes)

    integral = dot(weights, calc_ddcs_ria.(E_f, mu))

    return integral
end

# Gauss_Lobatto
function GaussLobatto_integrate_ddcs_over_mu(E_f)

    mu, weights = gausslobatto(nodes)
    
    integral = dot(weights, calc_ddcs_ria.(E_f, mu))

    return integral
end


# array of E_f values
E_f_vals = range(0.7*E, E, length=200) # paper starts at 0.7E

# x axis
xvals = E_f_vals ./ E
# y axis
quadgk_yvals = [quadgk_integrate_ddcs_over_mu(E_f) for E_f in E_f_vals]
GaussLegendre_yvals = [GaussLegendre_integrate_ddcs_over_mu(E_f) for E_f in E_f_vals]
GaussChebyshev1_yvals = [GaussChebyshev1_integrate_ddcs_over_mu(E_f) for E_f in E_f_vals]
GaussChebyshev2_yvals = [GaussChebyshev2_integrate_ddcs_over_mu(E_f) for E_f in E_f_vals]
GaussChebyshev3_yvals = [GaussChebyshev3_integrate_ddcs_over_mu(E_f) for E_f in E_f_vals]
GaussChebyshev4_yvals = [GaussChebyshev4_integrate_ddcs_over_mu(E_f) for E_f in E_f_vals]
GaussJacobi_yvals = [GaussJacobi_integrate_ddcs_over_mu(E_f, alpha, beta) for E_f in E_f_vals]
GaussRadau_yvals = [GaussRadau_integrate_ddcs_over_mu(E_f) for E_f in E_f_vals]
#GaussLobatto_yvals = [GaussLobatto_integrate_ddcs_over_mu(E_f) for E_f in E_f_vals] # diverges when few nodes

Plots.plot(
    xvals, quadgk_yvals,
    xlabel = "E' / E",
    ylabel = "DCS (cm²)",
    label = "quadgk",
    title = "DCS(E') for Al, E = $E_keV keV, $nodes nodes",
    legend=:bottom
)

Plots.plot!(xvals, GaussLegendre_yvals, label = "GaussLegendre")
Plots.plot!(xvals, GaussChebyshev1_yvals, label = "GaussChebyshev1")
Plots.plot!(xvals, GaussChebyshev2_yvals, label = "GaussChebyshev2")
Plots.plot!(xvals, GaussChebyshev3_yvals, label = "GaussChebyshev3")
Plots.plot!(xvals, GaussChebyshev4_yvals, label = "GaussChebyshev4")
a=round(alpha, digits=2); b=round(beta, digits=2)
Plots.plot!(xvals, GaussJacobi_yvals, label = "GaussJacobi, α=$a, β=$b")
Plots.plot!(xvals, GaussRadau_yvals, label = "GaussRadau")
#Plots.plot!(xvals, GaussLobatto_yvals,label = "GaussLobatto")
