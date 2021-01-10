@assert begin
    local_interp=InterpolationFunction(
        :(a0+a1*xi+a2*eta),
        [:xi, :eta],
        [:a0, :a1, :a2],
        [
            [0, 0]
        ]
    )

    abs_interp=InterpolationFunction(
        :(a0+a1*x+a2*y),
        [:x, :y],
        [:a0, :a1, :a2],
        [
            [0, 0]
        ]
    )

    dmn=Domain(
        local_interp,
        [
            (0.0, 1.0),
            (0.0, :(1.0-xi))
        ],
        [
            [0.0, 0.0],
            [1.0, 0.0],
            [0.0, 1.0]
        ]
    )

    fin=Finel(
        [
            [1.0, 1.0],
            [2.0, 1.0],
            [2.0, 2.0]
        ],
        dmn
    )

    vals=[
        0.0, 0.0, 1.0
    ]

    regr=finel_get_regr(fin, abs_interp)

    coefs=regress(regr, vals)

    isright=isapprox(coefs, [-1.0, 0.0, 1.0])

    Jac=Matrix{Float64}(undef, 2, 2)
    from_finel_jac!(Jac, fin, 0.0, 0.0)
    
    isright && isapprox(
        Jac,
        [1.0 1.0; 0.0 1.0]
    )
end
