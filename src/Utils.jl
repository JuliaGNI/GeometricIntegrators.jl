__precompile__()

module Utils

    export @define, @reexport

    include("utils/macro_utils.jl")

    export istriustrict, istrilstrict, L2norm, l2norm
    export simd_scale!, simd_copy!, simd_copy_scale!,
           simd_copy_xy_first!, simd_copy_yx_first!, simd_copy_yx_second!,
           simd_copy_yx_first_last!,
           simd_xpy!, simd_axpy!, simd_wxpy!, simd_waxpy!,
           simd_aXbpy!, simd_abXpy!, simd_mult!

    include("utils/matrix_utils.jl")

end
