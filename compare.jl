using Pkg;Pkg.activate("./src/parse")
using JSON
using Plots
using HDF5
using Statistics
gr(fontfamily="Computer Modern",legend=:topright,frame=:box,titlefontsize=11,legendfontsize=9,labelfontsize=12,left_margin=0Plots.mm)

path      = "./external_data/smeared/" 
ensembles = Dict(
    "M1" => "Sp4b6.5nF2nAS3mF-0.71mAS-1.01T48L20",
    "M2" => "Sp4b6.5nF2nAS3mF-0.71mAS-1.01T64L20",
    "M3" => "Sp4b6.5nF2nAS3mF-0.71mAS-1.01T96L20",
    "M4" => "Sp4b6.5nF2nAS3mF-0.7mAS-1.01T64L20",
    "M5" => "Sp4b6.5nF2nAS3mF-0.72mAS-1.01T64L32",
)

function mass_and_decay_constant_unrenormalized(path,ens,rep,channel)
    json_file = joinpath(path,ens,"meson_extraction_$(rep)_$(channel)_samples.json")
    data = JSON.parsefile(json_file)
    f    = mean(data["$(rep)_$(channel)_matrix_element_samples"]./sqrt.(data["$(rep)_$(channel)_mass_samples"]))
    Δf   = std( data["$(rep)_$(channel)_matrix_element_samples"]./sqrt.(data["$(rep)_$(channel)_mass_samples"]))
    m    = data["$(rep)_$(channel)_mass_value"]
    Δm   = std(data["$(rep)_$(channel)_mass_samples"])
    return f, Δf, m, Δm
end

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

plt = plot(xticks=(1:5,["M$i" for i in 1:5]))
for i in 1:5
    ens = "M$i"
    fPS, ΔfPS = mass_and_decay_constant_unrenormalized(path,ensembles[ens],"f" ,"ps")[1:2]
    fps, Δfps = mass_and_decay_constant_unrenormalized(path,ensembles[ens],"as","ps")[1:2]
    fV,  ΔfV  = mass_and_decay_constant_unrenormalized(path,ensembles[ens],"f" ,"v" )[1:2]
    fv,  Δfv  = mass_and_decay_constant_unrenormalized(path,ensembles[ens],"as","v" )[1:2]
    scatter!(plt,[i],[fPS], yerr=[ΔfPS], label="", alpha=0.8, shape=:rect,      color=:green)
    scatter!(plt,[i],[fps], yerr=[Δfps], label="", alpha=0.8, shape=:circ,      color=:blue)
    scatter!(plt,[i],[fV],  yerr=[ΔfV],  label="", alpha=0.8, shape=:hexagon,   color=:pink)
    scatter!(plt,[i],[fv],  yerr=[Δfv],  label="", alpha=0.8, shape=:utriangle, color=:orange)
end
plt