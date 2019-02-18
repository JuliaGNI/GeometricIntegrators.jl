
using DecFP
using DecFP: Dec128, _parse


macro define(name, definition)
    quote
        macro $(esc(name))()
            esc($(Expr(:quote, definition)))
        end
    end
end


function _dec128(x::Expr)
    y = x
    y.args .= _dec128.(y.args)
    return  y
end

_dec128(x::Int) = Dec128(x)
_dec128(x::String) = _parse(Dec128, x)

_dec128(x) = x

macro dec128(x)
    return esc(_dec128(x))
end
