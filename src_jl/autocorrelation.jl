using Pkg;Pkg.activate("./src_jl")
using MadrasSokal
using Plots
using LaTeXStrings
using DelimitedFiles
using HDF5
using Statistics
default(fontfamily="Computer Modern", frame=:box, top_margin=4Plots.mm, left_margin=6Plots.mm, plot_titlefontsize=12)

function full_observable_from_hdf5(file, outfile; label, name, plot_label, outdir, index=nothing, group="", therm = 1)
    dir1 = joinpath(outdir,"$(label)_history")
    dir2 = joinpath(outdir,"$(label)_full")
    for dir in [dir1,dir2]
        ispath(dir) || mkpath(dir)
    end

    f = h5open(file)
    for ensemble in keys(f)
        ens = joinpath(ensemble,group)
        obs = read(f[ens],name)
        if !isnothing(index)
            obs = obs[index...]
        end
        cfgn = 1:length(obs)
        
        plt1,τ, Δτ,τexp = MadrasSokal.publication_plot(cfgn,obs,plot_label,therm)
        plt2 = autocorrelation_overview(cfgn,obs,plot_label,therm;with_exponential=true)
 
        β    = read(f[ens],"beta")
        T, L = read(f[ens],"lattice")[1:2]
        mas  = read(f[ens],"quarkmassesmas")[1]
        mf   = read(f[ens],"quarkmassesmf")[1]
        title = latexstring(L"\beta\!=\!%$(β),~ T\!\times\!L^3\!=\!\!%$(T)\!\times\!%$(L)^3,-\!(am_0^{\rm f},am_0^{\rm as})\!=\!(%$(mf),%$(mas))")
        
        h5write(outfile,joinpath(ensemble,label,"tau_exp"),τexp)
        h5write(outfile,joinpath(ensemble,label,"tau_int"),τ)
        h5write(outfile,joinpath(ensemble,label,"Delta_tau_int"),Δτ)

        plot!(plt1,plot_title=title,size=(800,300))  
        plot!(plt2,plot_title=title)  
        savefig(plt1,joinpath(dir1,ensemble*".pdf"))
        savefig(plt2,joinpath(dir2,ensemble*".pdf"))
    end
end
# Write results into csv file
fmt(x,Δx) = errorstring(x,Δx)
fmt(x)    = string(round(x,sigdigits=2))
function write_tau_csv(h5file,flowdata,out_tex)
    io  = open(out_tex,"w")
    f = h5open(h5file)
    g = h5open(flowdata)

    header = L"""
            \begin{tabular}{|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|}
            \hline \hline
            Label & $\beta$ & $am_0^{\rm as}$ & $am_0^{\rm f}$ & $N_t$ & $N_s$ & $N_{\rm therm}$ & $n_{\rm skip}$ & $N_{\rm conf}$ & $\langle P \rangle$ & $w_0 / a$ & $\tau_{\rm int}^{\langle P \rangle}$ & $\tau_{\rm int}^{w_0}$ & $\tau_{\rm int}^{\rm ps}$ & $\tau_{\rm int}^Q$ & $\bar{Q}$ \\ 
            \hline
            """
    footer = """
            \\hline \\hline
            \\end{tabular}
            """
    write(io,header)
    for e in keys(f)
        β    = read(g[e],"beta")
        T, L = read(g[e],"lattice")[1:2]
        mas  = read(g[e],"quarkmassesmas")[1]
        mf   = read(g[e],"quarkmassesmf")[1]
        plaq = read(g[e],"plaquette")
        topo = read(g[e],"Q")
        w0   = fmt(read(g[e],"w0_val"),read(g[e],"w0_std"))
        P    = fmt(mean(plaq),std(plaq)/sqrt(length(plaq)))
        Q    = fmt(mean(topo),std(topo)/sqrt(length(topo)))
        N1   = last(split(read(g[e],"configurations")[1],"n"))
        N2   = last(split(read(g[e],"configurations")[2],"n"))
        n_skip  = parse(Int,N2) - parse(Int,N1)
        n_therm = parse(Int,N1)
        n_conf  = length(read(g[e],"configurations"))
        τQ, ΔτQ = read(f[e]["topology"],"tau_int"),       read(f[e]["topology"],"Delta_tau_int")
        τE, ΔτE = read(f[e]["energy_density"],"tau_int"), read(f[e]["energy_density"],"Delta_tau_int")
        τP, ΔτP = read(f[e]["plaquette"],"tau_int"),      read(f[e]["plaquette"],"Delta_tau_int")
        τπ, Δτπ = read(f[e]["PS_correlator"],"tau_int"),  read(f[e]["PS_correlator"],"Delta_tau_int")
        println(io,"$e & $β & -$mf & -$mas & $T & $L & $n_therm & $n_skip & $n_conf & $P & $w0 & $(fmt(τP,ΔτP)) & $(fmt(τQ,ΔτQ)) & $(fmt(τE,ΔτE)) & $(fmt(τπ,Δτπ)) & $Q \\\\")
    end
    write(io,footer)
    close(io)
end
function main()
    file    = "data_assets/topology.hdf5"
    fileCB  = "data_assets/wall_correlators.hdf5" 
    outfile = "data_assets/autocor.hdf5"
    pltdir  = "data_assets/autocorrelation_plots"
    out_tex = "assets/ensembles.tex"

    isfile(outfile) && rm(outfile)
    ispath(dirname(out_tex)) || mkpath(dirname(out_tex))

    full_observable_from_hdf5(file,  outfile; label="topology",      outdir = pltdir, group=""    ,name="Q",                       plot_label="Q",                 therm = 1)
    full_observable_from_hdf5(file,  outfile; label="energy_density",outdir = pltdir, group=""    ,name="energy_density_w0_sym",   plot_label=L"\mathcal{E}(w_0)", therm = 1)
    full_observable_from_hdf5(file,  outfile; label="plaquette",     outdir = pltdir, group=""    ,name="plaquette",               plot_label=L"<\!p\!>",          therm = 1)
    full_observable_from_hdf5(fileCB,outfile; label="PS_correlator", outdir = pltdir, group="AS"  ,name="TRIPLET/g5",index=(:,10), plot_label=L"C_\pi(t=10)",      therm = 1)

    write_tau_csv(outfile,file,out_tex)
end
main()