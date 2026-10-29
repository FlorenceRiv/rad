using Plots, CSV, DataFrames

ria = ["Au_cs_RIA_Au_0to1MeV_10KeVbinWidth_2026-07-22_01_34_30_636_E1p0MeV_2026-07-22_12_38_26_747_dose.csv", "Au_cs_RIA_Au_0to1MeV_10KeVbinWidth_2026-07-22_01_34_30_636_E1p0MeV_2026-07-22_12_38_26_747_photon_spectrum.csv"]
no_ria = ["Au_cs_noRIA_Au_0to1MeV_10KeVbinWidth_2026-07-21_23_14_06_341_E1p0MeV_2026-07-22_13_44_08_032_dose.csv", "Au_cs_noRIA_Au_0to1MeV_10KeVbinWidth_2026-07-21_23_14_06_341_E1p0MeV_2026-07-22_13_44_08_032_photon_spectrum.csv"]

function plotthis(version)

    dose = CSV.read(version[1], DataFrame)
    spectrum = CSV.read(version[2], DataFrame)

    p1 = Plots.plot(dose[:, 1], dose[:, 2], 
        title="Dose", xlabel="Depth (cm)", ylabel="Dose", legend=false)

    p2 = Plots.plot(spectrum[:, 1], Matrix(spectrum[:, 2:end]), 
        title="Spectrum", xlabel="Energy (MeV)", ylabel="Flux", 
        labels=reshape(names(spectrum)[2:end], 1, :))

    Plots.plot(p1, p2, layout=(2, 1), size=(800, 800))
end

#plotthis(ria)
plotthis(no_ria)

