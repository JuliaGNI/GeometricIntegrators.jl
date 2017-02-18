
function B1(t,q)
   return dA3d2(t,q) - dA2d3(t,q)
end

function B2(t,q)
   return dA1d3(t,q) - dA3d1(t,q)
end

function B3(t,q)
   return dA2d1(t,q) - dA1d2(t,q)
end


function B(t,q)
   return sqrt(B1(t,q)^2 + B2(t,q)^2 + B3(t,q)^2)
end


function b1(t,q)
   return B1(t,q) / B(t,q)
end

function b2(t,q)
   return B2(t,q) / B(t,q)
end

function b3(t,q)
   return B3(t,q) / B(t,q)
end
