"""
Demonstration that integrating over the electron direction cosine μ_e from 0 to 1 
covers the entire photon direction cosine μ_γ from 1 to -1.
"""

using Plots

function convert_mu_electron_to_mu_photon(μ_electron::Float64, Ei_photon::Float64)

    μ_photon = ( (1-μ_electron^2)*(1+Ei_photon)^2 - μ_electron^2 ) / ( (1-μ_electron^2)*(1+Ei_photon)^2 + μ_electron^2 )
    
    return μ_photon

end    


function plot_compton_kinematic_mapping(Ei::Float64)

    # Define integration domain for electron cosine (0.0 to 1.0)
    mu_e_range = range(0.0, 1.0, length=100)
    
    # Initialize empty photon cosine array
    mu_gamma_range = zeros(Float64, length(mu_e_range))
    
    # Calculate mu_gamma for every mu_e in the integration domain
    for (i, mu_e) in enumerate(mu_e_range)
        mu_gamma_range[i] = convert_mu_electron_to_mu_photon(mu_e, Ei)
    end
    
    # Generate plot
    p = plot(mu_e_range, mu_gamma_range, 
        label="E_i = $(Ei) mc^2",
        xlabel="Electron Angle Cosine (μ_e)", 
        ylabel="Photon Angle Cosine (μ_γ)",
        title="Kinematic Mapping: μ_e to μ_γ",
        lw=2, 
        legend=:bottomleft,
        grid=true,
        framestyle=:box
    )
    
    scatter!(p, [0.0, 1.0], [1.0, -1.0], color=:red, label="Kinematic Limits")
    
    display(p)
end

# test Ei = 1 mc^2
plot_compton_kinematic_mapping(1.0)