using Pkg;Pkg.activate("./src/parse")
using HiRepParsing
using HiRepOutputCleaner
using DelimitedFiles

# This script parses the log files in the directory 'dir', and saves them as an hdf5-file 
# in the location provided by 'h5file'.

# It creates a single hdf5 file for all log files. Measurements performed on the same ensemble
# are written in distinct hdf5 groups labelled  by the variable `ensemble`
function main(listfile,h5file;setup=true,filter_channels=false,channels=nothing)
    isfile(h5file) && rm(h5file)
    for (file,ensemble) in eachrow(readdlm(listfile,','))
        @show file,ensemble
        #regex = r"(DEFAULT SEMWALL TRIPLET|TRIPLET)"
        #writehdf5_spectrum_with_regexp(file,h5file,regex;mixed_rep=true,h5group=ensemble,setup,filter_channels,channels,sort=true)
        tmp_filename = "tmp/tmp.txt"
        clean_hirep_file(file,tmp_filename;checkpoint_pattern="analysed")
        # only try parsing if the filesize is non-vanishing
        if filesize(tmp_filename) > 0
            regex = r"(DEFAULT SEMWALL TRIPLET|TRIPLET)"
            writehdf5_spectrum_with_regexp(tmp_filename,h5file,regex;mixed_rep=true,h5group=ensemble,setup,filter_channels,channels,sort=true)
        end
        rm(tmp_filename)

    end
end

listfile = "./metadata/ensembles.csv" 
h5file   = "data_assets/test.hdf5" 
channels = ["g5, g1, g2, g3", "g5_g0g5_re"]
filter_channels = true

isfile(h5file) && rm(h5file)
main(listfile, h5file; filter_channels, channels)