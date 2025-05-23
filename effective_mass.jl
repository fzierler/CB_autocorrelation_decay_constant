using Pkg;Pkg.activate("./src/parse")
using LatticeUtils
using HDF5
using Plots
using LaTeXStrings
gr(fontfamily="Computer Modern",legend=:topright,frame=:box,titlefontsize=11,legendfontsize=9,labelfontsize=12,left_margin=0Plots.mm)

function effective_mass_plots(h5file,savedir)
    isdir(savedir) || mkpath(savedir)
    fid    = h5open(h5file)
    ensembles = keys(fid)
    for ens in ensembles
        eid  = fid[ens]
        T, L = fid[ens]["lattice"][1:2]
        
        corr_fπ = read(eid["TRIPLET"],"g5_g0g5_re")
        corr_g5 = read(eid["TRIPLET"],"g5")
        corr_gi = (read(eid["TRIPLET"],"g1") .+ read(eid["TRIPLET"],"g2") .+ read(eid["TRIPLET"],"g3"))./3

        corr_g5 = correlator_folding(corr_g5';t_dim=1)
        corr_gi = correlator_folding(corr_gi';t_dim=1)
        corr_fπ = correlator_folding(corr_fπ';t_dim=1,sign=-1)
        mg5, Δmg5 = implicit_meff(corr_g5) 
        mgi, Δmgi = implicit_meff(corr_gi)
        mfπ, Δmfπ = implicit_meff(corr_fπ)
        
        plt = plot(title="$ens: T=$T, L=$L, m_as=")
        scatter!(plt,mg5,yerr=Δmg5,label="pseudoscalar")
        scatter!(plt,mgi,yerr=Δmgi,label="vector")
        scatter!(plt,mfπ,yerr=Δmfπ,label="g5_g0g5")
        if endswith(ens,"AS")
            plot!(plt,ylims=(0.55,0.75))
        elseif endswith(ens,"FUN")
            plot!(plt,ylims=(0.3,0.50))
        end
        savefig(plt,joinpath(savedir,"$ens.pdf"))
    end
end
h5file  = "data_assets/test.hdf5"
savedir = "data_assets/effective_mass_plots/"
effective_mass_plots(h5file,savedir)
