@assert begin
    bfun=_get_gradfunc(:(1.0-xi-eta-zeta), [:xi, :eta, :zeta], 2)

    grad=bfun(zeros(Float64, 3)...)

    all(
        @. isapprox(grad, -1.0)
    )
end

@assert begin
    intp=InterpolationFunction(
        :(a0+xi*a1+eta*a2),
        [
            :xi, :eta
        ],
        [
            :a0, :a1, :a2
        ],
        [
            [0, 0]
        ]
    )

    dmn=Domain(
        intp,
        [
            (
                0.0, 1.0
            ),
            (
                0.0, :(1.0-xi)
            )
        ],
        [
            [0.0, 0.0],
            [1.0, 0.0],
            [0.0, 1.0]
        ]
    )

    vals=[1.0, 0.0, 0.0]
    coefs=regress(dmn, vals)

    isapprox(coefs, [1.0, -1.0, -1.0])
end
