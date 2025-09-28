#######################################################
# Progress Monitoring for Array Processing
#######################################################

using Dates
using Printf

"""
    ProgressTracker

Tracks progress of array processing operations with timing, throughput, and ETA calculations.

# Fields

  - `total_items::Int`: Total number of items to process
  - `processed_items::Ref{Int}`: Number of items processed so far
  - `start_time::Float64`: Start time of the operation
  - `last_update_time::Ref{Float64}`: Last time progress was updated
  - `update_interval::Float64`: Minimum seconds between progress updates
  - `show_progress::Bool`: Whether to display progress information
"""
mutable struct ProgressTracker
    total_items::Int
    processed_items::Ref{Int}
    start_time::Float64
    last_update_time::Ref{Float64}
    update_interval::Float64
    show_progress::Bool

    function ProgressTracker(
        total_items::Int; update_interval::Float64 = 1.0, show_progress::Bool = true
    )
        @assert total_items > 0 "total_items must be positive"
        @assert update_interval > 0 "update_interval must be positive"

        start_time = time()
        new(
            total_items, Ref(0), start_time, Ref(start_time), update_interval, show_progress
        )
    end
end

"""
    update_progress!(tracker::ProgressTracker, processed::Int; force::Bool=false)

Update the progress tracker with the current number of processed items.

# Arguments

  - `tracker::ProgressTracker`: The progress tracker to update
  - `processed::Int`: Current number of processed items
  - `force::Bool`: Force update even if within update interval
"""
function update_progress!(tracker::ProgressTracker, processed::Int; force::Bool = false)
    tracker.processed_items[] = processed
    current_time = time()

    if !tracker.show_progress
        return nothing
    end

    # Only update display if enough time has passed or forced
    if force || (current_time - tracker.last_update_time[]) >= tracker.update_interval
        tracker.last_update_time[] = current_time
        display_progress(tracker)
    end
end

"""
    display_progress(tracker::ProgressTracker)

Display current progress information.
"""
function display_progress(tracker::ProgressTracker)
    processed = tracker.processed_items[]
    total = tracker.total_items
    elapsed = time() - tracker.start_time

    # Calculate progress percentage
    progress_pct = processed / total * 100

    # Calculate throughput
    throughput = processed / elapsed

    # Estimate time remaining
    if processed > 0
        eta_seconds = (total - processed) / throughput
        eta_str = format_duration(eta_seconds)
    else
        eta_str = "unknown"
    end

    # Create progress bar
    bar_width = 40
    filled = floor(Int, progress_pct * bar_width / 100)
    bar = "█"^filled * "░"^(bar_width - filled)

    # Format output
    elapsed_str = format_duration(elapsed)

    print(
        "\r[$bar] $(round(progress_pct, digits=1))% " *
        "($processed/$total) " *
        "Elapsed: $elapsed_str " *
        "ETA: $eta_str " *
        "Speed: $(round(throughput, digits=1)) items/s",
    )

    if processed >= total
        println()  # New line when complete
    end

    flush(stdout)
end

"""
    format_duration(seconds::Float64)

Format duration in seconds to human-readable string.
"""
function format_duration(seconds::Float64)
    if seconds < 60
        return @sprintf("%.1fs", seconds)
    elseif seconds < 3600
        mins = floor(Int, seconds / 60)
        secs = seconds - mins * 60
        return @sprintf("%dm%.1fs", mins, secs)
    else
        hours = floor(Int, seconds / 3600)
        mins = floor(Int, (seconds - hours * 3600) / 60)
        return @sprintf("%dh%dm", hours, mins)
    end
end

"""
    with_progress(func, items; show_progress=true, desc="Processing")

Execute a function on array items with progress tracking.

# Arguments

  - `func`: Function to apply to each item
  - `items`: Array of items to process
  - `show_progress::Bool`: Whether to show progress bar
  - `desc::String`: Description for the progress bar

# Returns

  - Array of results from applying func to each item
"""
function with_progress(func, items; show_progress = true, desc = "Processing")
    n = length(items)
    if n == 0
        return similar(items, 0)
    end

    tracker = ProgressTracker(n; show_progress = show_progress)

    # Get the return type by calling func on the first item
    first_result = func(first(items))
    results = Vector{typeof(first_result)}(undef, n)
    results[1] = first_result

    show_progress && println("$desc $(n) items...")
    update_progress!(tracker, 1)

    for i in 2:n
        results[i] = func(items[i])
        update_progress!(tracker, i)
    end

    show_progress && update_progress!(tracker, n; force = true)

    return results
end

"""
    map_with_progress(func, items...; show_progress=true, desc="Processing")

Map function over arrays with progress tracking (similar to Base.map).

# Arguments

  - `func`: Function to apply
  - `items...`: Arrays to map over
  - `show_progress::Bool`: Whether to show progress bar
  - `desc::String`: Description for the progress bar

# Returns

  - Array of results from applying func to items
"""
function map_with_progress(func, items...; show_progress = true, desc = "Processing")
    if isempty(items)
        return []
    end

    n = length(first(items))
    @assert all(length(arr) == n for arr in items) "All arrays must have the same length"

    if n == 0
        return []
    end

    tracker = ProgressTracker(n; show_progress = show_progress)

    # Get the return type by calling func on the first elements
    first_args = tuple((arr[1] for arr in items)...)
    first_result = func(first_args...)
    results = Vector{typeof(first_result)}(undef, n)
    results[1] = first_result

    show_progress && println("$desc $(n) items...")
    update_progress!(tracker, 1)

    for i in 2:n
        args = tuple((arr[i] for arr in items)...)
        results[i] = func(args...)
        update_progress!(tracker, i)
    end

    show_progress && update_progress!(tracker, n; force = true)

    return results
end
