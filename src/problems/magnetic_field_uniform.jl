
include("magnetic_field.jl")


const A₀ = 1

function A1(t,q)
   return - 0.5 * A₀ * q[2]
end

function A2(t,q)
   return + 0.5 * A₀ * q[1]
end

function A3(t,q)
   return zero(eltype(q))
end


function dA1d1(t,q)
   return zero(eltype(q))
end

function dA1d2(t,q)
   return - 0.5 * A₀
end

function dA1d3(t,q)
   return zero(eltype(q))
end


function dA2d1(t,q)
   return + 0.5 * A₀
end

function dA2d2(t,q)
   return zero(eltype(q))
end

function dA2d3(t,q)
   return zero(eltype(q))
end


function dA3d1(t,q)
   return zero(eltype(q))
end

function dA3d2(t,q)
   return zero(eltype(q))
end

function dA3d3(t,q)
   return zero(eltype(q))
end
