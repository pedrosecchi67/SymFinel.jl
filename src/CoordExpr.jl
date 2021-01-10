"""
Get a function representing an expression's value, given coordinates and linear parameters
"""
function get_basefun(ex::Expr, coords::Vector{Symbol}, linpars::Vector{Symbol})
    q=quote
        bfun=($(coords...), $(linpars...)) -> ( $(ex) )
    end

    try
        eval(q)

        return bfun
    catch
        throw(error("get_basefun:ERROR:error when parsing basis function. Check inputs"))
    end
end

"""
Get base function derivative of a given order in respect to each coordinate
"""
function get_basefun_derivative(ex::Expr, coords::Vector{Symbol}, linpars::Vector{Symbol}, deriv_inds::Vector{Int64})
    excp=copy(ex)
    dinds=copy(deriv_inds)

    for i=1:size(dinds, 1)
        while dinds[i]>0
            dinds[i]-=1

            excp=differentiate(excp, coords[i])
        end
    end

    q=quote
        bfun=($(coords...), $(linpars...)) -> ( $excp )
    end

    try
        eval(q)

        return bfun
    catch
        throw(error("get_basefun_derivative:ERROR:wrong parsing of basis function derivative"))
    end
end

"""
Get linear contribution of a linear parameter for a basis function, as a function of position
"""
function get_linpar_contribution(ex::Expr, coords::Vector{Symbol}, linpars::Vector{Symbol}, 
    influencer::Symbol, deriv_inds::Vector{Int64})
    excp=copy(ex)
    dinds=copy(deriv_inds)

    for i=1:size(dinds, 1)
        while dinds[i]>0
            dinds[i]-=1

            excp=differentiate(excp, coords[i])
        end
    end

    excp=differentiate(excp, influencer)

    q=quote
        bfun=($(coords...), $(linpars...)) -> ( $excp )
    end

    try
        eval(q)

        return bfun
    catch
        throw(error("get_linpar_contribution:ERROR:error parsing linear transformation coefficient for basis function. Please verify linearity"))
    end
end
