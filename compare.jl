using Pkg;Pkg.activate("./src/parse")
using JSON
using Statistics

path      = "/home/fabian/Downloads/JSONs" 
ensembles = [
    "Sp4b6.5nF2nAS3mF-0.71mAS-1.01T48L20",
    "Sp4b6.5nF2nAS3mF-0.71mAS-1.01T64L20",
    "Sp4b6.5nF2nAS3mF-0.71mAS-1.01T96L20",
    "Sp4b6.5nF2nAS3mF-0.7mAS-1.01T64L20",
    "Sp4b6.5nF2nAS3mF-0.72mAS-1.01T64L32",
]

for ens in ensembles
    json_file = joinpath(path,ens,"meson_extraction_f_ps_samples.json")
    data = JSON.parsefile(json_file)
    fπ   = data["f_ps_matrix_element_value"]
    Δfπ  = std(data["f_ps_matrix_element_samples"])
    mπ   = data["f_ps_mass_value"]
    Δmπ  = std(data["f_ps_mass_samples"])
    @show ens, fπ, Δfπ, mπ, Δmπ
end

using HDF5
h5file = "data_assets/test.hdf5" 
fid    = h5open(h5file)
Z(C, β, ΔΣ, Δ, P) = 1 + C * (ΔΣ + Δ) * (8/β) /(16π^2*P)
ZA(C,β,P) = Z(C, β, -12.82, -3, P)
ZV(C,β,P) = Z(C, β, -12.82, -7.75, P)

#beta   = 6.5
#P      = [ mean(fid["M$(i)FUN"]["plaquette"][]) for i in 1:5 ]
#ZA_FUN = ZA.(5/4,beta,P)
#ZV_FUN = ZV.(5/4,beta,P)
#ZA_AS  = ZA.(2  ,beta,P)
#ZV_AS  = ZV.(2  ,beta,P)
