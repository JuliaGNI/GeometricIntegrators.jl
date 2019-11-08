
@testset "$(rpad("Linear solvers",80))" begin

    A = [[+4.  +5.  -2.]
         [+7.  -1.  +2.]
         [+3.  +1.  +4.]]
    b = [-14., +42., +28.]
    x = [+4., -4., +5.]


    function test_lu_solver(solver, A, b, x)
        # for T in (Float64, ComplexF64, Float32, ComplexF32) # TODO
        for T in (Float64, ComplexF64)
            AT = convert(Array{T,2}, A)
            bT = convert(Array{T,1}, b)
            xT = convert(Array{T,1}, x)

            lu = solver(AT, bT)
            factorize!(lu)
            solve!(lu)
            @test lu.b == xT
        end
    end

    test_lu_solver(LUSolverLAPACK, A, b, x)
    test_lu_solver(LUSolver, A, b, x)

end
