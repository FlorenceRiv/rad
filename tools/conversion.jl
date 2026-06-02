mₑc² = 0.510999 # MeV

function keV_to_mec2(E_keV)
    # keV → MeV , MeV → mₑc²
    return E_keV / 1e3 / mₑc²
end

function eV_to_mec2(E_eV)
    # eV → MeV , MeV → mₑc²
    return E_eV / 1e6 / mₑc²
end

function mec2_to_keV(E_mec2)
    # mₑc² → MeV, MeV → keV
    return E_mec2 * 1e3 * mₑc²
end

function mec2_to_eV(E_mec2)
    # mₑc² → MeV, MeV → eV
    return E_mec2 * 1e6 * mₑc²
end
