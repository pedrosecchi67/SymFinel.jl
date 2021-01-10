using LinearAlgebra

export Domain

# Geometric tolerance
geps=1e-18

"""
Struct to define a finite element domain

* `interp`: `InterpolationFunction` instance
* `restrictions`: definitions of domain limits in terms of auxiliary coordinates (`interp.coords`)
"""
mutable struct Domain
    interp::InterpolationFunction
    restrictions::Vector{
        Tuple{Any, Any}
    }
    _restriction_funcs::Vector{
        Tuple{Function, Function}
    }
    _restriction_derivatives::Vector{
        Tuple{Function, Function}
    }
    _regressor::Regressor
end

"""
Function to obtain restriction function
"""
function _get_restrfunc(re::Union{Expr, Number}, coords::Vector{Symbol}, ncoord::Int64)
    q=quote
        bfun=($(coords...),) -> $(re)
    end

    try
        eval(q)

        return bfun
    catch
        throw(error("_get_restrfunc:ERROR:unable to obtain restriction function"))
    end
end

"""
Function to obtain restriction gradient function
"""
function _get_gradfunc(re::Expr, coords::Vector{Symbol}, ncoord::Int64)
    der_exprs=:([$(differentiate(re, coords)...)])
    
    q=quote
        bfun=($(coords...),) -> $(der_exprs)
    end

    try
        eval(q)

        return bfun
    catch
        throw(error("_get_gradfunc:ERROR:unable to obtain restriction gradient. Please check for correctness with respect to doc. directives"))
    end
end

"""
Alternative implementation for literals
"""
function _get_gradfunc(re::Number, coords::Vector{Symbol}, ncoord::Int64)
    q=quote
        bfun=($(coords...),) -> ( 0.0 )
    end

    try
        eval(q)

        return bfun
    catch
        throw(error("_get_gradfunc:ERROR:error parsing restriction expression of literal type"))
    end
end

"""
Constructor for Domain struct

* `interp`: `InterpolationFunction` instance
* `restrictions`: definitions of domain limits in terms of auxiliary coordinates (`interp.coords`),
    which are negative within the domain and positive within it
"""
function Domain(
    interp::InterpolationFunction,
    restrictions,
    points::Vector{Vector{Float64}}
)
    ncoords=length(restrictions)

    _restriction_derivatives=[
        (_get_gradfunc(:($(vr[1])-$(c)), interp.coords, i), _get_gradfunc(:($(c)-$(vr[2])), interp.coords, i))
            for (i, (c, vr)) in enumerate(zip(interp.coords, restrictions))
    ]

    _restriction_functions=[
        (_get_restrfunc(vr[1], interp.coords, i), _get_restrfunc(vr[2], interp.coords, i))
            for (i, vr) in enumerate(restrictions)
    ]

    _restrictions=[
        Tuple{Any, Any}(vr)
                for vr in restrictions
    ]

    Domain(
        interp,
        _restrictions,
        _restriction_functions,
        _restriction_derivatives,
        Regressor(
            interp,
            points
        )
    )
end

"""
Same as `regress(regr::Regressor, Vector{Float64})`, but using a domain's regressor
"""
function regress(dmn::Domain, vals::Vector{Float64})
    return regress(dmn._regressor, vals)
end

"""
Get projection of a vector onto another (`u` onto `v`)
"""
function _getproj(u::Vector{Float64}, v::Vector{Float64})
    v*(dot(u, v)/(norm(v)^2+geps))
end

"""
Get determinant module by orthonormalizing the matrix at hand
"""
function orthonormal_det!(M::Matrix{Float64})
    d=1.0

    nd=size(M, 2)

    locn=1.0
    for i=1:nd
        for j=1:(i-1)
            M[:, i].-=dot(M[:, j], M[:, i]).*M[:, j]
        end

        locn=norm(M[:, i])+geps

        M[:, i]./=locn

        d*=locn
    end

    d
end

"""
Get component normal to the Jacobian space
"""
function normal_volspace!(M::Matrix{Float64}, nrestr::Int64, restr_derivs, isfirst::Bool)
    d=1.0

    nd=size(M, 2)

    locn=1.0
    for i=1:nd
        if i!=nrestr
            M[:, i].-=M[:, nrestr].*(restr_derivs[i]/restr_derivs[nrestr])

            for j=1:(i-1)
                if j!=nrestr
                    M[:, i].-=dot(M[:, j], M[:, i]).*M[:, j]
                end
            end

            locn=norm(M[:, i])+geps

            M[:, i]./=locn

            d*=locn
        end
    end

    for i=1:nd
        if i!=nrestr
            M[:, nrestr].-=dot(M[:, i], M[:, nrestr]).*M[:, i]
        end
    end

    M[:, nrestr].*=d/(norm(M[:, nrestr])+geps)

    if isfirst
        @. M[:, nrestr]=-M[:, nrestr]
    end
    
    M[:, nrestr]
end
