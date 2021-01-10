@assert begin
    intp=InterpolationFunction(
        :(a0+a1*x+a2*y+a11*x*y),
        [
            :x, :y
        ],
        [
            :a0, :a1, :a2, :a11
        ],
        [
            [0, 0]
        ]
    )

    mat=get_interp_linmat(
        intp,
        [
            [0.0, 0.0],
            [1.0, 0.0],
            [1.0, 1.0],
            [0.0, 1.0]
        ]
    )

    isapprox(
        [1.0 0.0 0.0 0.0; 1.0 1.0 0.0 0.0; 1.0 1.0 1.0 1.0; 1.0 0.0 1.0 0.0],
        mat
    )
end
