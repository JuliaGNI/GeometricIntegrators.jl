
abstract NonlinearSolver{T}

solve!(s::NonlinearSolver) = error("solve! not implemented for $(typeof(s))")
