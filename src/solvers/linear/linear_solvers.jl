
abstract LinearSolver

factorize!(s::LinearSolver) = error("factorize! not implemented for $(typeof(s))")
solve!(s::LinearSolver) = error("solve! not implemented for $(typeof(s))")
