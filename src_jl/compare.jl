using JSON
using Plots
using HDF5
using Statistics
using LaTeXStrings
using DelimitedFiles
using LatticeUtils
using ArgParse
pgfplotsx(size=(500, 300), legend=:topright,frame=:box,titlefontsize=12,legendfontsize=12,labelfontsize=12,left_margin=0Plots.mm)

function mass_and_decay_constant_unrenormalized(path,ens,rep,channel)
    json_file = joinpath(path,ens,"meson_extraction_$(rep)_$(channel)_samples.json")
    data = JSON.parsefile(json_file)
    f    = mean(data["$(rep)_$(channel)_matrix_element_samples"]./sqrt.(data["$(rep)_$(channel)_mass_samples"]))
    Δf   = std( data["$(rep)_$(channel)_matrix_element_samples"]./sqrt.(data["$(rep)_$(channel)_mass_samples"]))
    m    = data["$(rep)_$(channel)_mass_value"]
    Δm   = std(data["$(rep)_$(channel)_mass_samples"])
    return f, Δf, m, Δm
end
function plot_decay(h5file,wall_fit_path,smeared_path)
    fid    = h5open(h5file)
    P      = [ mean(fid["M$(i)/FUN"]["plaquette"][]) for i in 1:5 ]
    beta   = 6.5
    Z(C, β, ΔΣ, Δ, P) = 1 + C * (ΔΣ + Δ) * (8/β) /(16π^2*P)
    ZA(C,β,P) = Z(C, β, -12.82, -3, P)
    ZV(C,β,P) = Z(C, β, -12.82, -7.75, P)
    ZA_FUN = ZA.(5/4,beta,P)
    ZA_AS  = ZA.(2  ,beta,P)
    ZV_FUN = ZV.(5/4,beta,P)
    ZV_AS  = ZV.(2  ,beta,P)

    xticks = (1:5,["M$i" for i in 1:5])
    title  =  "Renormalized decay constants: smeared and local operators"
    xlabel =  "ensemble"
    ylabel = L"renormalized decay constants $[a]$"
    plt  = plot(;xticks, title, xlabel, ylabel)
    for i in 1:5
        ens = "M$i"
        fPS, ΔfPS = mass_and_decay_constant_unrenormalized(smeared_path,ENSEMBLES[ens],"f" ,"ps")[1:2]
        fps, Δfps = mass_and_decay_constant_unrenormalized(smeared_path,ENSEMBLES[ens],"as","ps")[1:2]
        fV,  ΔfV  = mass_and_decay_constant_unrenormalized(smeared_path,ENSEMBLES[ens],"f" ,"v" )[1:2]
        fv,  Δfv  = mass_and_decay_constant_unrenormalized(smeared_path,ENSEMBLES[ens],"as","v" )[1:2]
        scatter!(plt,[i],[ZA_FUN[i]*fPS],yerr=[ZA_FUN[i]*ΔfPS],label="", alpha=0.8, shape=:rect,      color=:green)
        scatter!(plt,[i],[ZA_AS[i]*fps], yerr=[ZA_AS[i]*Δfps], label="", alpha=0.8, shape=:circ,      color=:blue)
        scatter!(plt,[i],[ZV_FUN[i]*fV], yerr=[ZV_FUN[i]*ΔfV], label="", alpha=0.8, shape=:pentagon,  color=:pink)
        scatter!(plt,[i],[ZV_AS[i]*fv],  yerr=[ZV_AS[i]*Δfv],  label="", alpha=0.8, shape=:utriangle, color=:orange)
    end
    # add legend
    xl, yl = xlims(plt), ylims(plt)
    scatter!(plt,[0],[0],label=L"$af_{\rm PS}$ (smeared)", alpha=0.8, shape=:rect,      color=:green)
    scatter!(plt,[0],[0],label=L"$af_{\rm ps}$ (smeared)", alpha=0.8, shape=:circ,      color=:blue)
    scatter!(plt,[0],[0],label=L"$af_{\rm V}$  (smeared)", alpha=0.8, shape=:pentagon,  color=:pink)
    scatter!(plt,[0],[0],label=L"$af_{\rm v}$  (smeared)", alpha=0.8, shape=:utriangle, color=:orange)
    plot!(plt, xlims=xl,ylims=yl,legend=:outerright)

    for i in 1:5
        fPS, ΔfPS = readdlm(joinpath(wall_fit_path,"M$(i)FUN_ps.csv"),',',skipstart=1)[4:5]
        fps, Δfps = readdlm(joinpath(wall_fit_path,"M$(i)AS_ps.csv"),',',skipstart=1)[4:5]
        fV,  ΔfV  = readdlm(joinpath(wall_fit_path,"M$(i)FUN_v.csv"),',',skipstart=1)[4:5]
        fv,  Δfv  = readdlm(joinpath(wall_fit_path,"M$(i)AS_v.csv"),',',skipstart=1)[4:5]
        scatter!(plt,[i],[ZA_FUN[i]*fPS],yerr=[ZA_FUN[i]*ΔfPS],label="", alpha=0.4, shape=:rect,      color=:green)
        scatter!(plt,[i],[ZA_AS[i]*fps], yerr=[ZA_AS[i]*Δfps], label="", alpha=0.4, shape=:circ,      color=:blue)
        scatter!(plt,[i],[ZV_FUN[i]*fV], yerr=[ZV_FUN[i]*ΔfV], label="", alpha=0.4, shape=:pentagon,  color=:pink)
        scatter!(plt,[i],[ZV_AS[i]*fv],  yerr=[ZV_AS[i]*Δfv],  label="", alpha=0.4, shape=:utriangle, color=:orange)
    end
    # add legend
    xl, yl = xlims(plt), ylims(plt)
    scatter!(plt,[0],[0],label=L"$af_{\rm PS}$ (local)", alpha=0.4, shape=:rect,      color=:green)
    scatter!(plt,[0],[0],label=L"$af_{\rm ps}$ (local)", alpha=0.4, shape=:circ,      color=:blue)
    scatter!(plt,[0],[0],label=L"$af_{\rm V}$  (local)", alpha=0.4, shape=:pentagon,  color=:pink)
    scatter!(plt,[0],[0],label=L"$af_{\rm v}$  (local)", alpha=0.4, shape=:utriangle, color=:orange)
    plot!(plt, xlims=xl,ylims=yl,legend=:outerright)
    return plt
end
function table(file,h5file,wall_fit_path,smeared_path)
    io     = open(file,"w+")
    fid    = h5open(h5file)
    
    P      = [ mean(fid["M$(i)/FUN"]["plaquette"][]) for i in 1:5 ]
    beta   = 6.5
    Z(C, β, ΔΣ, Δ, P) = 1 + C * (ΔΣ + Δ) * (8/β) /(16π^2*P)
    ZA(C,β,P) = Z(C, β, -12.82, -3, P)
    ZV(C,β,P) = Z(C, β, -12.82, -7.75, P)
    
    header = "ens,ZA_FUN,ZA_AS,ZV_FUN,ZV_AS,fPS_s,ΔfPS_s,fps_s,Δfps_s,fV_s,ΔfV_s,fv_s,Δfv_s,fPS_l,ΔfPS_l,fps_l,Δfps_l,fV_l,ΔfV_l,fv_l,Δfv_l"
    write(io,header*"\n")
    for i in 1:5
        ens = "M$i"
        ZA_FUN = ZA(5/4,beta,P[i])
        ZA_AS  = ZA(2  ,beta,P[i])
        ZV_FUN = ZV(5/4,beta,P[i])
        ZV_AS  = ZV(2  ,beta,P[i])
        fPS_s, ΔfPS_s = mass_and_decay_constant_unrenormalized(smeared_path,ENSEMBLES[ens],"f" ,"ps")[1:2]
        fps_s, Δfps_s = mass_and_decay_constant_unrenormalized(smeared_path,ENSEMBLES[ens],"as","ps")[1:2]
        fV_s,  ΔfV_s  = mass_and_decay_constant_unrenormalized(smeared_path,ENSEMBLES[ens],"f" ,"v" )[1:2]
        fv_s,  Δfv_s  = mass_and_decay_constant_unrenormalized(smeared_path,ENSEMBLES[ens],"as","v" )[1:2]
        fPS_l, ΔfPS_l = readdlm(joinpath(wall_fit_path,"M$(i)FUN_ps.csv"),',',skipstart=1)[4:5]
        fps_l, Δfps_l = readdlm(joinpath(wall_fit_path,"M$(i)AS_ps.csv"),',',skipstart=1)[4:5]
        fV_l,  ΔfV_l  = readdlm(joinpath(wall_fit_path,"M$(i)FUN_v.csv"),',',skipstart=1)[4:5]
        fv_l,  Δfv_l  = readdlm(joinpath(wall_fit_path,"M$(i)AS_v.csv"),',',skipstart=1)[4:5]
        writedlm(io,[ens ZA_FUN ZA_AS ZV_FUN ZV_AS fPS_s ΔfPS_s fps_s Δfps_s fV_s ΔfV_s fv_s Δfv_s fPS_l ΔfPS_l fps_l Δfps_l fV_l  ΔfV_l fv_l  Δfv_l],',')
    end
    close(io)
end
function tex_table(file,outfile)
    header = L"""
    \begin{tabular}{ |c|c|c||c|c||c|c||c|c||c|c| }
    \hline
    Ensemble &  $af_{\rm PS}^{loc}$ &  $af_{\rm PS}^{smear}$ &  $af_{\rm ps}^{loc}$ &  $af_{\rm ps}^{smear}$ &  $af_{\rm V}^{loc}$ &  $af_{\rm V}^{smear}$ &  $af_{\rm v}^{loc}$ &  $af_{\rm v}^{smear}$ \\
    \hline \hline
    """
    footer = """
    \\hline \\hline
    \\end{tabular}
    """
    io = open(outfile,"w+")
    write(io,header)
    for row in eachrow(readdlm(file,',',skipstart=1))
        ZA_FUN, ZA_AS, ZV_FUN, ZV_AS = row[2:5] 
        ens  = row[1]
        PS_s = errorstring(ZA_FUN*row[6],ZA_FUN*row[7]  ; nsig=1)
        ps_s = errorstring(ZA_AS *row[8],ZA_AS *row[9]  ; nsig=1)
        V_s  = errorstring(ZV_FUN*row[10],ZV_FUN*row[11]; nsig=1)
        v_s  = errorstring(ZV_AS *row[12],ZV_AS *row[13]; nsig=1)
        PS_l = errorstring(ZA_FUN*row[14],ZA_FUN*row[15]  ; nsig=1)
        ps_l = errorstring(ZA_AS *row[16],ZA_AS *row[17]  ; nsig=1)
        V_l  = errorstring(ZV_FUN*row[18],ZV_FUN*row[19]; nsig=1)
        v_l  = errorstring(ZV_AS *row[20],ZV_AS *row[21]; nsig=1)
        println(io,"$ens & $PS_l & $PS_s & $ps_l & $ps_s & $V_l & $V_s & $v_l & $v_s \\\\") 
    end
    write(io,footer)
    close(io)
end
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--wall_correlators_h5"
            help = "HDF5 file containing the wall source correlators"
            required = true
        "--wall_fits"
            help = "Directory containing the results of the wall source correlators fits"
            required = true
        "--smeared_results"
            help = "Directory containing the results of the smeared correlators"
            required = true
        "--decay_output_csv"
            help = "Output location for CSV file for the comparison of the  decay constants"
            required = true
        "--decay_output_tex"
            help = "Output location for LaTex file containing the publication-ready table for the comparison of decay constants"
            required = true
        "--decay_output_pdf"
            help = "Output location for PDF file containing the plot of the comparison of decay constants"
            required = true
    end
    return parse_args(s)
end
function main()
    parsed_args   = parse_commandline()
    h5file        = parsed_args["wall_correlators_h5"]
    wall_fit_path = parsed_args["wall_fits"]
    smeared_path  = parsed_args["smeared_results"]
    file          = parsed_args["decay_output_csv"]
    file_tex      = parsed_args["decay_output_tex"]
    plot_file     = parsed_args["decay_output_pdf"]

    plt = plot_decay(h5file,wall_fit_path,smeared_path)
    ispath("assets") || mkpath("assets")
    savefig(plt,plot_file)
    table(file,h5file,wall_fit_path,smeared_path)
    tex_table(file, file_tex)
end

# Dictionary for looking up data from smeared correlators
const ENSEMBLES = Dict(
    "M1" => "Sp4b6.5nF2nAS3mF-0.71mAS-1.01T48L20",
    "M2" => "Sp4b6.5nF2nAS3mF-0.71mAS-1.01T64L20",
    "M3" => "Sp4b6.5nF2nAS3mF-0.71mAS-1.01T96L20",
    "M4" => "Sp4b6.5nF2nAS3mF-0.7mAS-1.01T64L20",
    "M5" => "Sp4b6.5nF2nAS3mF-0.72mAS-1.01T64L32",
)
main()