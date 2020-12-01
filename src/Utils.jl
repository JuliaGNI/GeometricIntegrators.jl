module Utils

    export rel_err, l2_err

    include("utils/diag_utils.jl")

    export @define, @dec128, @big

    include("utils/macro_utils.jl")

    export compensated_summation, L2norm, l2norm

    include("utils/sum_utils.jl")

    export istriustrict, istrilstrict
    
    include("utils/matrix_utils.jl")

end
