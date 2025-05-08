using LatticeUtils
using HDF5
using Plots
using LaTeXStrings
gr(fontfamily="Computer Modern",legend=:topright,frame=:box,titlefontsize=11,legendfontsize=9,labelfontsize=12,left_margin=0Plots.mm)

h5file = "data_assets/test.hdf5" 
fid    = h5open(h5file)
m, Δm  = log_meff(corr')

ensembles = keys(fid)
for ens in ensembles
    eid  = fid[ens]
    T, L = fid[ensembles[1]]["lattice"][1:2]
    
    corr_fπ = read(eid["TRIPLET"],"g5_g0g5_re")
    corr_g5 = read(eid["TRIPLET"],"g5")
    corr_gi = (read(eid["TRIPLET"],"g1") .+ read(eid["TRIPLET"],"g2") .+ read(eid["TRIPLET"],"g3"))./3
    mg5, Δmg5 = implicit_meff(corr_g5') 
    mgi, Δmgi = implicit_meff(corr_gi')
    mfπ, Δmfπ = implicit_meff(corr_fπ')
    
    plt = plot(title="$ens: T=$T, L=$L, m_as=")
    scatter!(plt,mg5,yerr=Δmg5,label="pseudoscalar")
    scatter!(plt,mgi,yerr=Δmgi,label="vector")
    scatter!(plt,mfπ,yerr=Δmfπ,label="g5_g0g5")
    plot!(plt,ylims=(0,1))
    display(plt)
end
