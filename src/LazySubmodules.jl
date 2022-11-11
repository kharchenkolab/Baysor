# import Random.randstring

# mutable struct LazySubmodule
#     _include_path::String
#     _module::Union{Module, Nothing}
# end

# # Later it should check for environment variable to decide whether to import immediately
# # lazy_submodule(include_path::String) = LazySubmodule(include_path, nothing)
# function lazy_submodule(include_path::String)
#     parent = Base.source_path(nothing)
#     include_path = (parent === nothing) ? abspath(include_path) : normpath(joinpath(dirname(parent), include_path))
#     return LazySubmodule(include_path, nothing)
# end

# function rename_mod(expr::Expr, new_name::Symbol)
#     if expr.head == :module
#         expr.args[2] = new_name
#     end
#     return expr
# end

# function Base.getproperty(m::LazySubmodule, s::Symbol)
#     if s in (:_include_path, :_module)
#         return getfield(m, s)
#     end

#     if m._module === nothing
#         mod_name = Symbol("mod_$(randstring(16))")
#         # m._module = include(e -> rename_mod(e, mod_name), m._include_path)
#         m._module = include(m._include_path)
#     end

#     return getfield(m._module, s)
# end

module LazySubmodules

export @lazy_submodule

const _LAZYMODE = Ref(true)

mutable struct LazySubmodule
    _module_name::Symbol
    _include_path::String
    _module::Union{Module, Nothing}
end

macro lazy_submodule(ex)
    name, path = ex.args
    parent = Base.source_path(nothing)
    path = (parent === nothing) ? abspath(path) : normpath(joinpath(dirname(parent), path))

    if !_LAZYMODE[]
        Core.eval(__module__, :($(name) = include($(path))))
    else
        s = LazySubmodule(name, String(path), __module__)
        Core.eval(__module__, :($(name) = $s))
    end
end

function Base.getproperty(m::LazySubmodule, s::Symbol)
    if s in (:_include_path, :_module, :_module_name)
        return getfield(m, s)
    end

    mod = Core.eval(m._module, :($(m._module_name) = include($(m._include_path))))
    return getfield(mod, s)
end

function __init__()
    if get(ENV, "LazyModules_lazyload", "") == "false"
        _LAZYMODE[] = false
    end
end

end