using Pkg;Pkg.activate("./src/parse")
using HiRepParsing


# This script parses the log files in the directory 'dir', and saves them as an hdf5-file 
# in the location provided by 'h5file'.

# It creates a single hdf5 file for all log files. Measurements performed on the same ensemble
# are written in distinct hdf5 groups labelled  by the variable `ensemble`
function main(listfile,h5file;setup=true,filter_channels=false,channels=nothing)
    isfile(h5file) && rm(h5file)
    for (file,ensemble) in eachrow(readdlm(listfile,','))
        @show file,ensemble
        regex = r"TRIPLET"
        writehdf5_spectrum_with_regexp(file,h5file,regex;mixed_rep=true,h5group=ensemble,setup,filter_channels,channels,sort=true)
    end
end
main("./metadata/ensembles.csv","data_assets/test.hdf5")
