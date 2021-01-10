export Finel

"""
Struct defining a finite element

* `pts`: corner points (vector of vectors)
* `domain`: `Domain` struct with the finite element's generic definition
* `abs_coord_regressors`: dictionary mapping `InterpolationFunction` structs to
    regressors
"""
mutable struct Finel
    pts::Vector{
        Vector{Float64}
    }
    domain::Domain
    abs_coord_regressors::Dict{
        InterpolationFunction, Regressor
    }
    _coord_linpars::Vector{Vector{Float64}}
end

"""
Constructor for finite element

* `pts`: corner points (vector of vectors)
* `domain`: `Domain` struct with the finite element's generic definition
"""
function Finel(pts::Vector{Vector{Float64}}, domain::Domain)
    pts_transpose=Vector{Vector{Float64}}(undef, length(pts[1]))
    for i=1:length(pts_transpose)
        pts_transpose[i]=[
            pts[j][i] for j=1:length(pts)
        ]
    end

    Finel(
        pts,
        domain,
        Dict{InterpolationFunction, Regressor}(),
        [
            regress(domain, coords) for coords in pts_transpose
        ]
    )
end

"""
Either create a new regressor for the given corner points or fetch a pre-existing one from the
struct's dictionary
"""
function finel_get_regr(fin::Finel, interp::InterpolationFunction)
    if haskey(fin.abs_coord_regressors, interp)
        return fin.abs_coord_regressors[interp]
    end

    regr=Regressor(
        interp,
        fin.pts
    )

    fin.abs_coord_regressors[interp]=regr

    regr
end

"""
Function to get absolute coordinates from finite element
"""
function get_abs_coords!(absc::Vector{Float64}, fin::Finel, xi::Vector{Float64})
    for (i, lps) in enumerate(fin._coord_linpars)
        absc[i]=fin.domain.interp._funcs[1](xi..., lps...)
    end
end
