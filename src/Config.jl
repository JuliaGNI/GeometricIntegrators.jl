__precompile__()

module Config

    export add_config, get_config, set_config

    if !isdefined(:GICONFIG)
        global GICONFIG = Dict()
    end

    function add_config(name, value)
        if !haskey(GICONFIG, name)
            GICONFIG[name] = value
        end
    end

    function set_config(name, value)
        if haskey(GICONFIG, name)
            GICONFIG[name] = value
        else
            println("  WARNING: Unknown parameter name.")
        end
    end

    function get_config(name)
        if haskey(GICONFIG, name)
            return GICONFIG[name]
        else
            println("  WARNING: Unknown parameter name.")
            return nothing
        end
    end

end
