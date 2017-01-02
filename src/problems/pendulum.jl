
function pendulum_ode_f(t, x, f)
    f[1] = x[2]
    f[2] = sin(x[1])
    nothing
end

function pendulum_ode(x0=[acos(0.4), 0.0])
    ODE(pendulum_ode_f, x0)
end


function pendulum_pode_v(t, q, p, v)
    v[1] = q[1]
    nothing
end

function pendulum_pode_f(t, q, p, f)
    f[1] = sin(q[1])
    nothing
end

function pendulum_pode(q0=[acos(0.4)], p0=[0.0])
    PODE(pendulum_pode_v, pendulum_pode_f, q0, p0)
end
