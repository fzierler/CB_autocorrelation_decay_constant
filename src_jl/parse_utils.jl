function gaugegroup(file)
    if HDF5.ishdf5(file) 
        hdf5 = h5open(file, "r")
        return read(hdf5,"gauge group")
    else 
        return gaugegroup_log(file)
    end
end
function quarkmasses(file)
    if HDF5.ishdf5(file) 
        hdf5 = h5open(file, "r")
        return read(hdf5,"quarkmasses")
    else 
        return quarkmasses_log(file)
    end
end
function quarkmasses_chimera(file)
    if HDF5.ishdf5(file) 
        hdf5 = h5open(file, "r")
        mf = read(hdf5,"quarkmasses_fundamental")
        mas = read(hdf5,"quarkmasses_antisymmetric")
        return mf, mas
    else 
        return quarkmasses_chimera_log(file)
    end
end
function APE_smearing(file)
    if HDF5.ishdf5(file) 
        hdf5 = h5open(file, "r")
        APE_eps = read(hdf5,"APE_eps")
        APE_level = read(hdf5,"APE_levels")
        return APE_eps, APE_level
    else 
        return APE_smearing_logfile(file)
    end
end
function Wuppertal_smearing_mixed(file)
    if HDF5.ishdf5(file) 
        hdf5 = h5open(file, "r")
        antisymmetric_eps = read(hdf5,"Wuppertal_eps_anti")
        fundamental_eps   = read(hdf5,"Wuppertal_eps_fund")
        return antisymmetric_eps, fundamental_eps
    else
        return Wuppertal_smearing_mixed_logfile(file)
    end
end
function latticesize(file)
    if HDF5.ishdf5(file) 
        hdf5 = h5open(file, "r")
        return read(hdf5,"lattice")
    else 
        return latticesize_log(file)
    end
end
function plaquettes(file)
    if HDF5.ishdf5(file) 
        hdf5 = h5open(file, "r")
        return read(hdf5,"plaquette")
    else 
        return plaquettes_log(file)
    end
end
function inverse_coupling(file)
    if HDF5.ishdf5(file) 
        hdf5 = h5open(file, "r")
        return read(hdf5,"beta")
    else 
        return inverse_coupling_log(file)
    end
end
function correlators(file,type,key;kws...)
    if HDF5.ishdf5(file) 
        hdf5 = h5open(file, "r")
        return read(hdf5,key)
    else 
        return correlators_logfile(file,type,key;kws...)
    end
end 
function correlators(file,type,key,nhits;kws...)
    if HDF5.ishdf5(file) 
        hdf5 = h5open(file, "r")
        corr = read(hdf5,key) 
        return corr 
    else 
        d = correlators_logfile(file,type,key;masses,mass,kws...)
        return reshape_disconnected(d,nhits;rescale=1)
    end
end
function reshape_disconnected(d::Matrix{S},nhits;rescale=1) where S
    # the rounding automatically removes incomplete calcluations
    T, N = size(d)
    n = div(N,nhits)
    m = zeros(eltype(d[1]),(n, nhits, T))
    for i in 1:n, j in 1:nhits, t in 1:T
        m[i,j,t] = d[t,(i-1)*nhits + j]
    end
    @. m = rescale*m
    return m
end
function gaugegroup_log(file)
    for line in eachline(file)
        if occursin("Gauge group",line)
            pos = findlast(' ',line)
            return strip(line[pos:end])
        end
    end
    return ""
end
function quarkmasses_log(file;pattern=r"\[MAIN\]\[0\](M|m)ass")
    masses = Float64[]
    for line in eachline(file)
        if occursin(pattern,line)
            s = split(line,(","))
            for i in eachindex(s)
                s0 = split(s[i],r"(=|:|;)")[2:end]
                m  = parse.(Float64,s0)
                append!(masses,m)
            end
            return unique(masses)
        end
    end
    return ""
end
function quarkmasses_chimera_log(file)
    mf  = quarkmasses_log(file;pattern="[MAIN][0]mf[0]")
    mas = quarkmasses_log(file;pattern="[MAIN][0]mas[0]")
    return mf, mas
end
function latticesize_log(file)
    for line in eachline(file)
        if occursin("Global size is",line)
            pos  = last(findfirst("Global size is",line))+1
            sizestring  = lstrip(line[pos:end])
            latticesize = parse.(Int,split(sizestring,"x"))
            return latticesize
        end
    end
    return ""
end
function plaquettes_log(file)
    plaquettes = Float64[]
    for line in eachline(file)
        if occursin("Plaquette",line)
            line = replace(line,"="=>" ")
            line = replace(line,":"=>" ")
            p = parse(Float64,split(line)[end])
            append!(plaquettes,p)
        end
    end
    return plaquettes
end
function _match_config_name(filename)
    regex = r".*/(?<run>[^/]*)_(?<T>[0-9]+)x(?<L>[0-9]+)x[0-9]+x[0-9]+nc[0-9]+(?:r[A-Z]+)?(?:nf[0-9]+)?b(?<beta>[0-9]+\.[0-9]+)?(?:m-?[0-9]+\.[0-9]+)?n(?<conf>[0-9]+)"
    return match(regex,filename)
end
function inverse_coupling_log(file)
    try
        l = split(file,"beta")[end]
        β = parse(Float64,split(l,"m")[1])
        return β
    catch
        for line in eachline(file)
            if occursin("Configuration from",line)
                match = _match_config_name(line)
                β = parse(Float64,match[:beta])
                return β
            end
        end
    end
end
function APE_smearing_logfile(file)
    APE_level = Int64[]
    APE_eps   = Float64[]
    # start with reference level for no smearing
    N = -1
    for line in eachline(file)
        if startswith(line,"[APE][0]APE smearing with val")
            eps = parse(Float64,last(split(line,"="))) 
            append!(APE_eps,eps)
        end
        if startswith(line,"[APE][0]N=")
            pos1 = length("[APE][0]N=") + 1
            pos2 = findfirst('<',line) - 1
            N0 = parse(Int,line[pos1:pos2]) 
            N  = max(N,N0)
        end
        if startswith(line,"[APE][0]APE smearing END")
            append!(APE_level,N)
            N = -1 # reset to reference value for no smearing
        end
    end
    return APE_eps, APE_level
end
function Wuppertal_smearing_mixed_logfile(file)
    fundamental_eps = Float64[]
    antisymmetric_eps = Float64[]
    
    #default reference values for no smearing
    epsA = 0.0
    epsF = 0.0

    for line in eachline(file)
        if startswith(line,"[SMEAR][0]source smearing epsilon =")
            pos1 = length("[SMEAR][0]source smearing epsilon =") + 1
            pos2 = first(findfirst("iterations:",line)) - 1
            epsA = parse(Float64,line[pos1:pos2])
        end
        if startswith(line,"[SMEAR][0]Fundamental source smearing epsilon =")
            pos1 = length("[SMEAR][0]Fundamental source smearing epsilon =") + 1
            pos2 = first(findfirst("iterations:",line)) - 1
            epsF = parse(Float64,line[pos1:pos2])
        end
        if occursin("analysed", line)
            append!(antisymmetric_eps,epsA)
            append!(fundamental_eps,epsF)
        end
    end
    return antisymmetric_eps, fundamental_eps
end
function correlators_logfile(file,type,key;kws...)
    corrs = parse_spectrum(file,type;filter_channels=true,channels=[key],kws...)
    return reduce(hcat,getindex.(corrs,key))
end
function parse_spectrum(file,type;disconnected=false,masses=false,mass="",filter_channels=false,channels="",nhits=1,with_progress=false,re_im=true)
    T = latticesize(file)[1]
    corr = zeros(T) # preallocate array for parsing of correlator
    dict = Dict{String,Vector{Float64}}()
    dictarray = Dict{String,Vector{Float64}}[]
    conf0 = 0
    src0  = 0
    # when filtering for specific keys also allow them to end with "_re" and "_im"
    if filter_channels
        all_channels = channels
        if re_im
            if disconnected
                all_channels = hcat(channels,channels.*"_disc_re",channels.*"_disc_im")
            else
                all_channels = hcat(channels,channels.*"_re",    channels.*"_im")
            end
        end
    end
    # keep track of position in file for progress meter
    with_progress && (p = Progress(countlines(file); dt=1, desc="Match $type: Progress:"))
    for line in eachline(file)
        with_progress && next!(p)
        if occursin(type,line)
            if masses
                occursin("mass=$mass",line) || continue
            end
            # get configuration number
            pos_num = findfirst('#',line)
            end_num = findnext(' ',line,pos_num)
            conf = parse(Int,line[pos_num+1:end_num-1])
            # find number of the source if available
            if disconnected
                pos_src = last(findfirst("src",line))+1
                end_src = findnext(' ',line,pos_src+1)
                src = parse(Int,line[pos_src:end_src])
            else
                src = 0
            end
            # find last '=' sign which separates values from Γ structure
            # TODO this does not work for momenta
            pos_eq = findlast('=',line)
            #key_st = findprev(' ',line,pos_eq)
            key_st = last(findfirst(type,line))+1
            key = line[key_st+1:pos_eq-1]
            if filter_channels
                if last(split(key,"/")) ∉ all_channels
                    continue
                end
            end
            if disconnected
                # create new entry if configuration or source number changes
                # if we need to parse more than one source at a time per configuration
                if conf0 != conf || src0 != src
                    if !isempty(dict)
                        push!(dictarray,dict)
                        dict = Dict{String,Vector{Float64}}()
                    end
                end
            end
            # parse corrrelator values
            pos_0 = findnext(' ',line,pos_eq)
            for t in 1:T
                pos_1 = findnext(' ',line,pos_0+1)
                corr[t] = Parsers.parse(Float64,line[pos_0:pos_1])
                pos_0 = pos_1
            end
            dict[key] = copy(corr)
            conf0 = conf
            src0  = src
        end
        if !disconnected
            # If we only have one source at a time and possibly one configuration
            # at a time: the method used to separate distinct measurements fails. 
            # In this case the end of measurement on a given confiuration is 
            # signalled by a line that reads:
            # [MAIN][0]Configuration #N: analysed in [a sec b usec]
            if occursin("analysed",line)
                if !isempty(dict)
                    push!(dictarray,dict)
                    dict = Dict{String,Vector{Float64}}()
                end
            end
        end
    end
    if !isempty(dict)
        push!(dictarray,dict)
    end
    return _reshape_connected(dictarray;disconnected,nhits)
end
function _reshape_connected(dict;disconnected=false,nhits=1)
    corrs = Dict{String,Array{Float64}}()
    for k in keys(dict[1])
        if disconnected
            corrs[k] = reshape_disconnected(reduce(hcat,getindex.(dict,k)),nhits)
            # rename keys for disconnected measurements
            corrs = _rename_disconnected_keys(corrs)
        else
            corrs[k] = permutedims(reduce(hcat,getindex.(dict,k)))
        end
    end
    return corrs
end
function _rename_disconnected_keys(corrs_discon)
    new_corrs_discon = Dict{String,Array{Float64}}()
    for k in keys(corrs_discon)
        val = corrs_discon[k]
        k_new = replace(k,"_re" => "", "_disc" =>"")
        new_corrs_discon[k_new] = val
    end
    return new_corrs_discon
end
function confignames(file)
    fns = AbstractString[]
    for line in eachline(file)
        if occursin("read",line)
            if occursin("Configuration",line)
                pos1 = findlast('/',line)
                pos2 = findnext(']',line,pos1)
                push!(fns,line[pos1+1:pos2-1])
            end
        end
    end
    return fns
end
function confignames_and_plaquette(file)
    fns = AbstractString[]
    plaquettes = Float64[]
    for line in eachline(file)
        if startswith(line,"[IO][0]Configuration")
            # parse filename first
            pos1 = findlast('/',line)
            pos2 = findnext(']',line,pos1)
            push!(fns,line[pos1+1:pos2-1])
            # then parse plaquette
            posP = findlast('=',line) + 1            
            p = parse(Float64,line[posP:end])
            append!(plaquettes,p)
        end
    end
    return fns, plaquettes
end
function permutation_names(names)
    numbers = parse.(Int,last.(split.(names,"n")))
    return sortperm(numbers)
end
unique_indices(v) = unique(i -> v[i], eachindex(v))
function _write_lattice_setup(file,h5file;mixed_rep=false,h5group="",sort=false,smearing=true,deduplicate=false)
    names, plaq = confignames_and_plaquette(file)
    perm  = sort ? permutation_names(names) : eachindex(names)
    inds  = deduplicate ? unique_indices(names[perm]) : eachindex(names[perm]) 
    plaq  = plaq[perm][inds]
    names = names[perm][inds]
    # save other relevant quantities
    h5write(h5file,joinpath(h5group,"plaquette"),plaq)
    h5write(h5file,joinpath(h5group,"configurations"),names)
    h5write(h5file,joinpath(h5group,"gauge group"),gaugegroup(file))
    h5write(h5file,joinpath(h5group,"beta"),inverse_coupling(file))
    h5write(h5file,joinpath(h5group,"lattice"),latticesize(file))
    # write information on the applied sorting and deduplication to file
    h5write(h5file,joinpath(h5group,"sorted"),sort)
    h5write(h5file,joinpath(h5group,"deduplicated"),deduplicate)
    h5write(h5file,joinpath(h5group,"sort_permutation"),collect(perm))
    h5write(h5file,joinpath(h5group,"deduplicated_indices"),collect(inds))
    # get smearing parameters (arrays are empty if no smearing is used)
    if smearing
        APE_eps, APE_level = APE_smearing(file)
        Wuppertal_eps_anti, Wuppertal_eps_fund = Wuppertal_smearing_mixed(file)
        h5write(h5file,joinpath(h5group,"APE_eps"),APE_eps)
        h5write(h5file,joinpath(h5group,"APE_level"),APE_level)
        h5write(h5file,joinpath(h5group,"Wuppertal_eps_anti"),Wuppertal_eps_anti)
        h5write(h5file,joinpath(h5group,"Wuppertal_eps_fund"),Wuppertal_eps_fund)
    end
    # special case fermion masses for mixed representations
    if !mixed_rep
        h5write(h5file,joinpath(h5group,"quarkmasses"),quarkmasses(file))
    else
        mf, mas = quarkmasses_chimera(file)
        h5write(h5file,joinpath(h5group,"quarkmasses_fundamental"),mf)
        h5write(h5file,joinpath(h5group,"quarkmasses_antisymmetric"),mas)
    end
end
function writehdf5_spectrum(file,h5file,type::AbstractString;sort=false,h5group="",setup=true,mixed_rep=false,h5group_setup=h5group,filter_channels=false,channels=nothing,re_im=true,kws...)
    names = confignames(file)
    perm  = sort ? permutation_names(names) :  collect(eachindex(names))
    setup && _write_lattice_setup(file,h5file;mixed_rep,h5group=h5group_setup,sort)
    # read correlator data
    c = parse_spectrum(file,type;disconnected=false,filter_channels,channels,re_im)
    # write matrices to file
    for Γ in keys(c)
        label = joinpath(h5group,type,Γ)
        filter_channels && Γ ∉ channels && continue
        sort && h5write(h5file,label,c[Γ][perm,:];kws...)
        sort || h5write(h5file,label,c[Γ];kws...)
    end
end