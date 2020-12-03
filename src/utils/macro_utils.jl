
using DecFP
using DecFP: Dec128, _parse


macro define(name, definition)
    quote
        macro $(esc(name))()
            esc($(Expr(:quote, definition)))
        end
    end
end


_dec128(x) = x
_dec128(x::Number) = Dec128(x)
_dec128(x::String) = _parse(Dec128, x)

function _dec128(x::Expr)
    y = x
    y.args .= _dec128.(y.args)
    return  y
end

macro dec128(x)
    return esc(_dec128(x))
end


_big(x) = x
_big(x::Number) = big(x)
_big(x::String) = parse(BigFloat, x)

function _big(x::Expr)
    y = x
    y.args .= _big.(y.args)
    return  y
end

macro big(x)
    return esc(_big(x))
end
