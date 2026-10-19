# Numerical helpers that reproduce R's arithmetic and number printing, so that
# results and log files match the R implementation exactly.
#
# R accumulates `sum`, `rowSums` and `mean` of doubles in `long double`
# (80-bit x87 extended precision, 64-bit mantissa, on x86-64). Emulated here
# with 64-bit-precision BigFloat. Integer data are summed exactly.

const _LDOUBLE_PRECISION = 64

"""
    _rsum(xs) -> Float64

Sum with R semantics: sequential accumulation in extended precision.
"""
function _rsum(xs)
    setprecision(BigFloat, _LDOUBLE_PRECISION) do
        s = BigFloat(0)
        for x in xs
            s += x
        end
        Float64(s)
    end
end

_rsum_exact(xs::AbstractVector{<:Integer}) = Float64(sum(Int128, xs; init = Int128(0)))

"""
    _rmean(xs) -> Float64

`mean()` of a double vector as computed by R (`summary.c`): extended-precision
sum divided by `n`, followed by one refinement pass. `xs` is any iterable of
`Float64` that can be iterated twice (e.g. a generator), so no vector needs
to be allocated.
"""
function _rmean(xs)
    setprecision(BigFloat, _LDOUBLE_PRECISION) do
        s = BigFloat(0)
        n = 0
        for x in xs
            s += x
            n += 1
        end
        n == 0 && return NaN
        s /= n
        if isfinite(Float64(s))
            t = BigFloat(0)
            for x in xs
                t += x - s
            end
            s += t / n
        end
        Float64(s)
    end
end

"""
    _rformat(x::Real) -> String

`as.character(x)` for a double in R: up to 15 significant digits, fixed
notation unless scientific notation is shorter (R's `formatReal` with
`scipen = 0`).
"""
function _rformat(x::Real)
    x = Float64(x)
    isnan(x) && return "NaN"
    isinf(x) && return x > 0 ? "Inf" : "-Inf"
    x == 0 && return "0"
    s = @sprintf("%.14e", x)
    mant, ex = split(s, 'e')
    e = parse(Int, ex)
    digits = rstrip(replace(mant, "-" => "", "." => ""), '0')
    nsig = length(digits)
    neg = x < 0 ? 1 : 0
    width_sci = neg + nsig + (nsig > 1 ? 1 : 0) + (abs(e) >= 100 ? 5 : 4)
    if e >= 0
        decimals = max(0, nsig - e - 1)
        width_fixed = neg + e + 1 + (decimals > 0 ? decimals + 1 : 0)
    else
        decimals = nsig - e - 1
        width_fixed = neg + 2 + decimals
    end
    if width_fixed <= width_sci
        # R rounds halves away from zero when no decimals are printed.
        decimals == 0 && return @sprintf("%.0f", round(x, RoundNearestTiesAway))
        return _printf_fixed(decimals, x)
    else
        return Printf.format(Printf.Format(string("%.", nsig - 1, "e")), x)
    end
end

# `@sprintf("%.*f", d, x)` (dynamic precision is not available in all
# supported Julia versions).
_printf_fixed(d::Integer, x::Real) = Printf.format(Printf.Format(string("%.", d, "f")), x)
