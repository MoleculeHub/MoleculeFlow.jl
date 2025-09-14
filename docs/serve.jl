#!/usr/bin/env julia

# Simple script to serve the documentation locally
using HTTP
using Sockets

function serve_docs()
    port = 8000
    root = joinpath(@__DIR__, "build")

    if !isdir(root)
        println("Documentation not built yet. Run:")
        println("julia --project=. make.jl")
        return nothing
    end

    println("Serving documentation at: http://localhost:$port")
    println("Press Ctrl+C to stop")

    try
        HTTP.serve("localhost", port) do request
            # Simple static file server
            path = HTTP.URIs.unescapeuri(request.target)
            if path == "/"
                path = "/index.html"
            end

            filepath = joinpath(root, path[2:end])  # Remove leading /

            if isfile(filepath)
                return HTTP.Response(200, read(filepath))
            else
                return HTTP.Response(404, "File not found")
            end
        end
    catch e
        if isa(e, InterruptException)
            println("\nStopping server...")
        else
            rethrow(e)
        end
    end
end

serve_docs()
