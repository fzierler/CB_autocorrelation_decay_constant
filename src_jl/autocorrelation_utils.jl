# Autocorrelation (exponential autocrrelation time)
function autocorrelation(x, lag)
    # (wasteful in terms of allocations but clear)
    z = x .- mean(x)
    norm = sum(abs2,z)
    a = sum(z[1+lag:end].*z[1:end-lag]) / norm
    return a
end
function autocorrelation(x;minlags=0)
    lx   = length(x)
    # lower limit of lags
    nlags = min(lx-1, round(Int,10log10(lx)))
    # apply a minimal number of lags
    nlags  = max(minlags,nlags)
    # create array of all lags considered
    lags  = collect(0:nlags)
    a = zeros(eltype(x),length(lags))
    for i in eachindex(lags)
        a[i] = autocorrelation(x, lags[i])
    end
    return a
end
function exponential_autocorrelation_time(O;minlags=0)
    # discard autocorrelation(t=0)=1 in fitting
    a = autocorrelation(O;minlags)[2:end]
    @. modelτ(x,p) = exp(-x/p[1])
    # we have previosuly discarded the data point at t=0
    x = collect(1:length(a))
    c = curve_fit(modelτ, x, a, ones(1))
    τ = c.param[1]
    return τ
end
function madras_sokal_estimator_fixedt(x, t;biased = false)
    Γ = zero(eltype(x))
    N = length(x)
    m = mean(x)
    v = var(x)
    τ = t - 1
    for i in 1:N-t
        norm = biased ? N : N-τ            
        Γ += (x[i]-m)*(x[i+τ]-m)/ norm / v
    end
    return Γ
end
# With the default maximal window size tmax=length(x)÷10,  
# the Madras-Sokal variance estimate is such that at tmax
# Δτ ÷ τ = 1/sqrt(2.5) ≈ 0.63
function madras_sokal_estimator_windows(x;max_window=length(x)÷10)
    Γ = zeros(eltype(x),max_window)
    for t in 1:max_window
        Γ[t] = madras_sokal_estimator_fixedt(x, t)
    end
    return Γ 
end
function madras_sokal_windows(x;kws...)
    Γ  = madras_sokal_estimator_windows(x;kws...)
    τ  = 1/2 .+ cumsum(Γ)
    Δτ = similar(τ)
    N  = length(x)
    for i in eachindex(τ)
        Δτ[i] = sqrt(τ[i]^2 * (4i+2)/N)
    end
    return τ, Δτ
end
function madras_sokal_time(x,therms;stop=length(x).-therms,kws...)
    τ  = zeros(Float64,size(therms))
    Δτ = zeros(Float64,size(therms))
    for j in eachindex(therms)
        # default value of stop is chosen, such that it is the end of x
        xc = x[therms[j]:therms[j]+stop[j]]
        τ_windows, Δτ_windows = madras_sokal_windows(xc;kws...)
        τ[j], W = findmax(τ_windows)
        Δτ[j]   = Δτ_windows[W]
    end
    return τ, Δτ
end
function madras_sokal_time(x;kws...)
    τ_windows, Δτ_windows = madras_sokal_windows(x;kws...)
    τ, W = findmax(τ_windows)
    Δτ   = Δτ_windows[W]
    return τ, Δτ
end
function errorstring(x,Δx;nsig=2)
    @assert Δx > 0
    sgn = x < 0 ? "-" : ""
    x = abs(x)  
    # round error part to desired number of signficant digits
    # convert to integer if no fractional part exists
    Δx_rounded = round(Δx,sigdigits=nsig) 
    # get number of decimal digits for x  
    floor_log_10 = floor(Int,log10(Δx))
    dec_digits   = (nsig - 1) - floor_log_10
    # round x, to desired number of decimal digits 
    # (standard julia function deals with negative dec_digits) 
    x_rounded = round(x,digits=dec_digits)
    # get decimal and integer part if there is a decimal part
    if dec_digits > 0
        digits_val = Int(round(x_rounded*10.0^(dec_digits)))
        digits_unc = Int(round(Δx_rounded*10.0^(dec_digits)))
        str_val = _insert_decimal(digits_val,dec_digits) 
        str_unc = _insert_decimal(digits_unc,dec_digits)
        str_unc = nsig > dec_digits ? str_unc : string(digits_unc)
        return sgn*"$str_val($str_unc)"
    else
        return sgn*"$(Int(x_rounded))($(Int(Δx_rounded)))"
    end
end
function _insert_decimal(val::Int,digits)
    str = lpad(string(val),digits,"0")
    pos = length(str) - digits
    int = rpad(str[1:pos],1,"0")
    dec = str[pos+1:end]
    inserted = int*"."*dec
    return inserted
end
