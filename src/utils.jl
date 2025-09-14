function pil_png_to_rgb(pyimg::Py)
    w, h = pyconvert(Tuple{Int, Int}, pyimg.size)
    img = pyconvert(Vector{UInt8}, pyimg.tobytes())
    img = reinterpret(NTuple{3, UInt8}, img)
    img = permutedims(reshape(img, w, h))
    return map(c -> RGB{N0f8}(c[1] / 255, c[2] / 255, c[3] / 255), img)
end

function extract_address(str::String)
    address = match(r"0x[0-9a-fA-F]+", str)

    if !isnothing(address)
        return address.match
    else
        return missing
    end
end

function julia_to_python_indices(julia_indices::Vector{Int})
    return [i - 1 for i in julia_indices]
end

function python_to_julia_indices(python_indices::Vector{Int})
    return [i + 1 for i in python_indices]
end

function julia_to_python_indices(julia_indices::Nothing)
    return nothing
end

function python_to_julia_indices(python_indices::Nothing)
    return nothing
end
