module LazySubmodules

# TODO: use include_dependency for Revise.jl (see https://docs.julialang.org/en/v1/base/base/)

export @lazy_submodule

const _LAZYMODE = Ref(true)
const _LOAD_LOCKER = Threads.ReentrantLock()

mutable struct LazySubmodule
    _alias::Symbol
    _module_name::Union{Symbol, Nothing}
    _include_path::String
    _parent_module::Module
end

macro lazy_submodule(ex)
    name, path = ex.args
    parent = Base.source_path(nothing)
    path = (parent === nothing) ? abspath(path) : normpath(joinpath(dirname(parent), path))

    if !_LAZYMODE[]
        Core.eval(__module__, :($(name) = include($(path))))
    else
        s = LazySubmodule(name, nothing, String(path), __module__)
        Core.eval(__module__, :($(name) = $s))
    end
end

function Base.getproperty(m::LazySubmodule, s::Symbol)
    if s in (:_alias, :_module_name, :_include_path, :_parent_module)
        return getfield(m, s)
    end

    f = nothing
    lock(_LOAD_LOCKER) do
        if m._module_name === nothing
            # Ideally, this function should be called only once per module
            # However in case of multi-threading or for keyword argument functons it can be called multiple times
            # See # https://discourse.julialang.org/t/how-to-properly-implement-lazy-compilation-of-submodules-problem-with-invokelatest-for-kwarg-functions/90459
            mod = Core.eval(m._parent_module, :(include($(m._include_path))))
            m._module_name = Symbol(split("$(mod)", ".")[end])
        else
            mod = getfield(m._parent_module, m._module_name)
        end

        f = getfield(mod, s)
    end

    if !(f isa Function)
        return f
    end

    # invokelatest solves issues with world age and missing functions
    return (args...; kwargs...) -> Base.invokelatest(f, args...; kwargs...)
    # return (args...; kwargs...) -> begin
        # res = Base.invokelatest(f, args...; kwargs...)
        # Substitude the LazySubmodule with an actual compiled module
        # lock(_LOAD_LOCKER) do
        #     Core.eval(m._parent_module, :($(m._alias) = $(m._module_name)))
        # end
        # return res
    # end
end

function __init__()
    if get(ENV, "LazyModules_lazyload", "") == "false"
        _LAZYMODE[] = false
    end
end

end