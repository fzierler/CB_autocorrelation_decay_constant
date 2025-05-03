using Pkg;Pkg.activate("./src/parse")
using HiRepParsing
using DelimitedFiles

# It creates a single hdf5 file for all log files. Measurements performed on the same ensemble
# are written in distinct hdf5 groups labelled  by the variable `ensemble`
function main(listfile,h5file;setup=true,filter_channels=false,channels=nothing)
    isfile(h5file) && rm(h5file)
    for (file,ensemble) in eachrow(readdlm(listfile,','))
        @show (file,ensemble)
        writehdf5_spectrum(file,h5file,"TRIPLET";mixed_rep=true,h5group=ensemble,setup,filter_channels,channels,sort=true,re_im=false)
    end
end

listfile = "./metadata/ensembles.csv" 
h5file   = "data_assets/test.hdf5" 
channels = ["g5", "g1", "g2", "g3", "g5_g0g5_re"]
filter_channels = true
isfile(h5file) && rm(h5file)
main(listfile, h5file; filter_channels, channels)