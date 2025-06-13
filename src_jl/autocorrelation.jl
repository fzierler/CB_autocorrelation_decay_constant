using Plots
using LaTeXStrings
using DelimitedFiles
using HDF5
using Statistics
using ArgParse
using LsqFit
include("autocorrelation_utils.jl")
default(fontfamily="Computer Modern", frame=:box, top_margin=4Plots.mm, left_margin=6Plots.mm, plot_titlefontsize=12)

function full_observable_from_hdf5(file, outfile; label, name, plot_label, outdir, index=nothing, group="", therm = 1)

    f = h5open(file)
    for ensemble in keys(f)
        ens = joinpath(ensemble,group)
        obs = read(f[ens],name)
        if !isnothing(index)
            obs = obs[index...]
        end
        
        τexp    = exponential_autocorrelation_time(obs;minlags=100)        
        τ0, Δτ0 = madras_sokal_windows(obs)
        τ, W    = findmax(τ0)
        Δτ      = Δτ0[W] 

        h5write(outfile,joinpath(ensemble,label,"tau_exp"),τexp)
        h5write(outfile,joinpath(ensemble,label,"tau_int"),τ)
        h5write(outfile,joinpath(ensemble,label,"Delta_tau_int"),Δτ)
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
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--wilson_flow_hdf5"
            help = "HDF5 file containing the gradient flow data"
            required = true
        "--wall_correlators_hdf5"
            help = "HDF5 file containing the correlators from local wall sources"
            required = true
        "--output_hdf5"
            help = "Where to write the resulting HDF5 file"
            required = true
        "--output_tex"
            help = "Where to write the publication ready LaTeX table"
            required = true
        "--plot_dir"
            help = "Where to save additional figures visualising the autocorrelation of all quantities"
            required = true
    end
    return parse_args(s)
end
function main()
    parsed_args = parse_commandline()

    file_flow = parsed_args["wilson_flow_hdf5"]
    file_wall = parsed_args["wall_correlators_hdf5"] 
    outfile = parsed_args["output_hdf5"]
    out_tex = parsed_args["output_tex"]
    pltdir  = parsed_args["plot_dir"]    

    isfile(outfile) && rm(outfile)
    ispath(dirname(out_tex)) || mkpath(dirname(out_tex))

    full_observable_from_hdf5(file_flow, outfile; label="topology",      outdir = pltdir, group=""    ,name="Q",                       plot_label="Q",                 therm = 1)
    full_observable_from_hdf5(file_flow, outfile; label="energy_density",outdir = pltdir, group=""    ,name="energy_density_w0_sym",   plot_label=L"\mathcal{E}(w_0)", therm = 1)
    full_observable_from_hdf5(file_flow, outfile; label="plaquette",     outdir = pltdir, group=""    ,name="plaquette",               plot_label=L"<\!p\!>",          therm = 1)
    full_observable_from_hdf5(file_wall, outfile; label="PS_correlator", outdir = pltdir, group="AS"  ,name="TRIPLET/g5",index=(:,10), plot_label=L"C_\pi(t=10)",      therm = 1)

    write_tau_csv(outfile,file_flow,out_tex)
end
main()