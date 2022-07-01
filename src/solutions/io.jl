
using Test

abstract type SolutionIO end

Base.close(sio::SolutionIO) = error("Base.close(::SolutionIO) not implemented for ", typeof(sio))

save_attributes(sio::SolutionIO, ::AbstractSolution) = error("save_attributes(::SolutionIO, ::Solution) not implemented for ", typeof(sio))
save_solution(sio::SolutionIO, ::AbstractSolution) = error("save_solution(::SolutionIO, ::Solution) not implemented for ", typeof(sio))
load_solution(sio::SolutionIO, ::AbstractSolution) = error("load_solution(::SolutionIO, ::Solution) not implemented for ", typeof(sio))

function test_interface(sio::SolutionIO, sol::AbstractSolution)
AbstractSolution
    @test_nowarn save_attributes(sio, sol)
    @test_nowarn save_solution(sio, sol)
    @test_nowarn load_solution(sio, sol)
    @test_nowarn close(sio)
    
end
