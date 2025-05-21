using Pkg;Pkg.activate("./src/parse")
using JSON
using Statistics

path      = "./external_data/smeared/" 
ensembles = [
    "Sp4b6.5nF2nAS3mF-0.71mAS-1.01T48L20",
    "Sp4b6.5nF2nAS3mF-0.71mAS-1.01T64L20",
    "Sp4b6.5nF2nAS3mF-0.71mAS-1.01T96L20",
    "Sp4b6.5nF2nAS3mF-0.7mAS-1.01T64L20",
    "Sp4b6.5nF2nAS3mF-0.72mAS-1.01T64L32",
]

function mass_and_decay_constant_unrenormalized(path,ens)
    json_file_π = joinpath(path,ens,"meson_extraction_f_ps_samples.json")
    json_file_ρ = joinpath(path,ens,"meson_extraction_f_v_samples.json")
    data_π = JSON.parsefile(json_file_π)
    data_ρ = JSON.parsefile(json_file_ρ)
    fπ   = mean(data_π["f_ps_matrix_element_samples"]./sqrt.(data_π["f_ps_mass_samples"]))
    Δfπ  = std( data_π["f_ps_matrix_element_samples"]./sqrt.(data_π["f_ps_mass_samples"]))
    fρ   = mean(data_ρ["f_v_matrix_element_samples"]./sqrt.(data_ρ["f_v_mass_samples"]))
    Δfρ  = std( data_ρ["f_v_matrix_element_samples"]./sqrt.(data_ρ["f_v_mass_samples"]))
    mπ   = data_π["f_ps_mass_value"]
    mρ   = data_ρ["f_v_mass_value"]
    Δmπ  = std(data_π["f_ps_mass_samples"])
    Δmρ  = std(data_ρ["f_v_mass_samples"])
    return mπ, Δmπ, mρ, Δmρ, fπ, Δfπ, fρ, Δfρ
end

using HDF5
h5file = "data_assets/test.hdf5" 
fid    = h5open(h5file)
P      = [ mean(fid["M$(i)FUN"]["plaquette"][]) for i in 1:5 ]
beta   = 6.5
Z(C, β, ΔΣ, Δ, P) = 1 + C * (ΔΣ + Δ) * (8/β) /(16π^2*P)
ZA(C,β,P) = Z(C, β, -12.82, -3, P)
ZV(C,β,P) = Z(C, β, -12.82, -7.75, P)
ZA_FUN = ZA.(5/4,beta,P)
ZV_FUN = ZV.(5/4,beta,P)
ZA_AS  = ZA.(2  ,beta,P)
ZV_AS  = ZV.(2  ,beta,P)
