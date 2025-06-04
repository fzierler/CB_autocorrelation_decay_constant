using Pkg;Pkg.activate("./src_jl")
using HiRepParsing
using DelimitedFiles

# It creates a single hdf5 file for all log files. Measurements performed on the same ensemble
# are written in distinct hdf5 groups labelled  by the variable `ensemble`
function main(listfile,h5file)
    setup=true
    filter_channels = true
    channels = ["g5", "g1", "g2", "g3", "g5_g0g5_re"]
    
    ispath(dirname(h5file)) || mkpath(dirname(h5file))
    isfile(h5file) && rm(h5file)
    
    println("Parsing wall source correlator data")
    for (file,ensemble) in eachrow(readdlm(listfile,','))
        @show (file,ensemble)
        m     = match(r"(M[1-5])(FUN|AS)",ensemble)
        group = joinpath(m.captures...)
        writehdf5_spectrum(file,h5file,"TRIPLET";mixed_rep=true,h5group=group,setup,filter_channels,channels,sort=true,re_im=false)
    end
end

listfile = "./metadata/ensembles.csv" 
h5file   = "data_assets/wall_correlators.hdf5" 
main(listfile, h5file)