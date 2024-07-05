
function update_vectorfields!(sol, problem::Union{AbstractProblemODE, SODEProblem})
    initialguess(problem).v(sol.v, sol.t, sol.q, parameters(problem))
    return sol
end

function update_vectorfields!(sol, problem::AbstractProblemPODE)
    initialguess(problem).v(sol.v, sol.t, sol.q, sol.p, parameters(problem))
    initialguess(problem).f(sol.f, sol.t, sol.q, sol.p, parameters(problem))
    return sol
end

function update_vectorfields!(sol, problem::AbstractProblemIODE)
    initialguess(problem).v(sol.v, sol.t, sol.q, sol.p, parameters(problem))
    initialguess(problem).f(sol.f, sol.t, sol.q, sol.v, parameters(problem))
    return sol
end
