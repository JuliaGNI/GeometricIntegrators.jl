
eachstage(int::GeometricIntegrator) = eachstage(method(int))
nstages(int::GeometricIntegrator) = nstages(tableau(method(int)))

"""
Create a vector of S solution vectors of type DT to store the solution of S
internal stages for a problem with `D` dimensions.
"""
function create_internal_stage_vector(DT, D, S)
    [zeros(DT, D) for i in 1:S]
end


"""
Create a vector of S solution matrices of type DT to store the solution of S
internal stages for a problem with `DxM` dimensions.
"""
function create_internal_stage_matrix(DT, D, M, S)
    [zeros(DT, D, M) for i in 1:S]
end


"""
Create a vector of S+1 solution vectors of type DT to store the solution of S
internal stages and the solution of the previous timestep for a problem with `D`
    dimensions.
"""
function create_internal_stage_vector_with_zero(DT, D, S)
    a = OffsetArray{Vector{DT}}(undef, 0:S)

    for i in 0:S
        a[i] = zeros(DT, D)
    end

    return a
end


"""
Create a vector of (S,M+1) solution vectors of type DT to store the solution of S
internal stages and M random processes for a problem with `D` dimensions.
"""
function create_internal_stage_vector_with_zero(DT, D, M, S)
    a = OffsetArray{Vector{DT}}(undef, 0:M, 1:S)

    for i in 0:M
        for j in 1:S
            a[i, j] = zeros(DT, D)
        end
    end

    return a
end
