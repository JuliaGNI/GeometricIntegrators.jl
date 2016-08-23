
abstract NonlinearSolver

solve!(s::NonlinearSolver) = error("solve! not implemented for $(typeof(s))")
