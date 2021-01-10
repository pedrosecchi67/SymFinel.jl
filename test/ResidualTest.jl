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

    resd=Residual(
        [
            (
                :u, abs_interp, [
                    (:du!dx, [1, 0]),
                    (:du!dy, [0, 1])
                ]
            ),
            (
                :v, abs_interp, []
            )
        ],
        :(k*(u*(du!dx-du!dy)-v)),
        3;
        extra_args=[:k]
    )

    surfresd=Residual(
        [
            (
                :u, abs_interp, [
                    (:du!dx, [1, 0]),
                    (:du!dy, [0, 1])
                ]
            ),
            (
                :v, abs_interp, []
            )
        ],
        :([k*du!dx, k*du!dy]),
        3;
        extra_args=[:k]
    )

    resdfunc=get_punctual_resfunc(resd)

    isright=isapprox(1.0, resdfunc(1.0, 2.0, 3.0, -2.0, 1.0))

    varvals=[
        [0.0, 0.0, 1.0],
        [0.0, 1.0, 0.0]
    ]

    pls=vals_to_polys(
        fin,
        resd,
        varvals
    )

    isright=isright && isapprox(
        pls[1], [-1.0, 0.0, 1.0]
    )
    isright=isright && isapprox(
        pls[2], [0.0, 1.0, -1.0]
    )

    vs=polys_to_vars(resd, pls, [1.5, 1.5])

    isright=isright && isapprox(vs, [0.5, 0.0, 1.0, 0.0])

    absc=Vector{Float64}(undef, 2)
    get_abs_coords!(absc, fin, [0.5, 0.5])

    isright=isright && isapprox(
        absc,
        [2.0, 1.5]
    )

    funresd=get_volume_residual_function(resd)

    R=funresd(varvals, fin, -2.0)
    @time R=funresd(varvals, fin, -2.0)
    @time R=funresd(varvals, fin, -2.0)
    @time R=funresd(varvals, fin, -2.0)
    @time R=funresd(varvals, fin, -2.0)
    @time R=funresd(varvals, fin, -2.0)
    @time R=funresd(varvals, fin, -2.0)
    @time R=funresd(varvals, fin, -2.0)

    isright=isright && isapprox(R, 2.0/3)

    surf_funresd=get_surface_residual_function(surfresd)

    R=surf_funresd(varvals, fin, -2.0)
    @time R=surf_funresd(varvals, fin, -2.0)
    @time R=surf_funresd(varvals, fin, -2.0)
    @time R=surf_funresd(varvals, fin, -2.0)
    @time R=surf_funresd(varvals, fin, -2.0)
    @time R=surf_funresd(varvals, fin, -2.0)
    @time R=surf_funresd(varvals, fin, -2.0)
    @time R=surf_funresd(varvals, fin, -2.0)

    isright=isright && isapprox(R, 0.0; atol=1e-10)

    isright
end
