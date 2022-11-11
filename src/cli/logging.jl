using Logging
import Dates
import Base.CoreLogging: shouldlog, min_enabled_level, catch_exceptions, handle_message
import Base.CoreLogging: Info, Warn

struct DoubleLogger <: AbstractLogger
    cli_stream::IO
    file_stream::Union{IO, Nothing}
    min_level::LogLevel
    message_limits::Dict{Any,Int}
    force_flush::Bool
end
DoubleLogger(cli_stream::IO=stderr, file_stream::Union{IO, Nothing}=nothing, level=Info; force_flush::Bool=false) =
    DoubleLogger(cli_stream, file_stream, level, Dict{Any,Int}(), force_flush)
shouldlog(logger::DoubleLogger, level, _module, group, id) = get(logger.message_limits, id, 1) > 0

min_enabled_level(logger::DoubleLogger) = logger.min_level

catch_exceptions(logger::DoubleLogger) = false

function handle_message(logger::DoubleLogger, level, message, _module, group, id,
                        filepath, line; maxlog=nothing, kwargs...)
    if maxlog !== nothing && maxlog isa Integer
        remaining = get!(logger.message_limits, id, maxlog)
        logger.message_limits[id] = remaining - 1
        remaining > 0 || return
    end

    log_message(logger.cli_stream, level, message, _module, filepath, line; force_flush=logger.force_flush)

    if logger.file_stream !== nothing
        log_message(logger.file_stream, level, message, _module, filepath, line; force_flush=logger.force_flush)
    end

    return nothing
end

function log_message(stream::IO, level, message, _module, filepath, line; force_flush::Bool)
    buf = IOBuffer()
    iob = IOContext(buf, stream)
    levelstr = level == Warn ? "Warning" : string(level)
    msglines = Base.split(chomp(string(message)), '\n')

    prefix = "[$(Dates.format(Dates.now(), "HH:MM:SS"))] $levelstr: "
    println(iob, prefix, msglines[1])

    msglines[2:end] = (" " ^ (length(prefix) - 1)) .* msglines[2:end]
    if level !== Info
        source_address = " " * something("$_module", "nothing") * " " *
            something("$filepath", "nothing") * ":" * something("$line", "nothing")
        push!(msglines, source_address)
    end

    for i in 2:length(msglines)
        println(iob, (i < length(msglines)) ? "|" : "└", msglines[i])
    end

    write(stream, take!(buf))

    if force_flush
        flush(stream)
    end
end