module Utils

    export rel_err

    include("utils/diag_utils.jl")

    export @define, @reexport, @dec128

    include("utils/macro_utils.jl")

    export compensated_summation, L2norm, l2norm

    include("utils/sum_utils.jl")

    export istriustrict, istrilstrict
    export simd_copy_xy_first!, simd_copy_yx_first!, simd_copy_yx_second!,
           simd_copy_yx_first_last!,
           simd_axpy!, simd_aXbpy!, simd_abXpy!, simd_mult!

    include("utils/matrix_utils.jl")

end
