module Utils

"""
`find_value` finds the range of `x` for which `f(x) = y0`.

It assumes that the function is continuous and monotonically increasing.
Under these constraints it guarantees (or at least tries to) that the value
is between the limits.

Keyword arguments:

    - `xeps` -- tolerance in x
    - `yeps` -- tolerance in y
    - `maxiters`
    - `verbose`
"""
function find_value(f::Function, y0::Number; xeps=1e-10, yeps=xeps, maxiters=100, verbose=false)
    if verbose println("[find_value: start] xeps=$xeps, yeps=$yeps") end

    xmin,xmax = _fv_find_range_double(f,y0,xeps,yeps,maxiters,verbose)
    if verbose println("[find_value: doubling] xmin=$xmin, xmax=$xmax") end
    converged,xmin,xmid,xmax = _fv_binary(xmin,xmax, f,y0,xeps,yeps,maxiters,verbose)
    if verbose println("[find_value: binary] xmin=$xmin, xmid=$xmid, xmax=$xmax ($(converged?"CONVERGED":"NOT CONVERGED"))") end
    xmin,xmax = if converged
        xmid, xmid
    else
        xmin = _fv_binary_lower(xmin,xmid, f,y0,xeps,yeps,maxiters,verbose)
        xmax = _fv_binary_upper(xmid,xmax, f,y0,xeps,yeps,maxiters,verbose)
        xmin, xmax
    end
    if verbose
        println("[find_value: final] xmin=$xmin, xmax=$xmax")
        println("[find_value: final] f(xmin)=$(f(xmin))")
        println("[find_value: final] f(xmax)=$(f(xmax))")
    end
    @assert abs(f(xmin) - y0) < yeps
    @assert abs(f(xmax) - y0) < yeps
    xmin,xmax
end

function _fv_binary_upper(xmin,xmax,f,y0,xeps,yeps,maxiters,verbose)
    xs = xmin,xmax
    y::Float64 = 0.0
    while abs(xs[2]-xs[1]) > xeps
        x = 0.5*(xs[2]+xs[1])
        y = f(x)

        @assert (y-y0) > -yeps
        xs = (abs(y-y0)<yeps) ? (x, xs[2]) : (xs[1], x)

        maxiters -= 1
        if maxiters == 0 error("Max iterations reached (binary_upper)") end
    end

    xs[1]
end

function _fv_binary_lower(xmin,xmax,f,y0,xeps,yeps,maxiters,verbose)
    xs = xmin,xmax
    y::Float64 = 0.0
    while abs(xs[2]-xs[1]) > xeps
        x = 0.5*(xs[2]+xs[1])
        y = f(x)

        @assert (y-y0) < yeps
        xs = (abs(y-y0)<yeps) ? (xs[1], x) : (x, xs[2])

        maxiters -= 1
        if maxiters == 0 error("Max iterations reached (binary_lower)") end
    end

    xs[2]
end

# returns: (converged::Bool, xmin, xmid, xmax)
# Convergence here means that the limit on the x-axis (xeps) has been reached.
# In that sense it is a bad thing in this context, since it also implies that
# the y value has not converged.
function _fv_binary(xmin,xmax,f,y0,xeps,yeps,maxiters,verbose)
    xs = xmin,xmax
    y::Float64 = 0.0
    while abs(xs[2]-xs[1]) > xeps
        x = 0.5*(xs[2]+xs[1])
        y = f(x)

        if abs(y-y0)<yeps
            return false, xs[1], x, xs[2]
        end
        xs = (y > y0) ? (xs[1], x) : (x, xs[2])

        maxiters -= 1
        if maxiters == 0 error("Max iterations reached (binary)") end
    end

    true, xs[1], 0.5*(xs[2]+xs[1]), xs[2]
end

function _fv_find_range_double(f,y0,xeps,yeps,maxiters,verbose)
    # The first step in the iteration to find the appropriate lambda
    # is to establish the initial search range. This is done by starting
    # from zero, ±1 and then doubling, until the target value A is somewhere
    # in that range.
    xmin,xmax = -1.0,1.0

    while true
        y = f(xmin)
        if y < y0 && abs(y-y0) > yeps
            break
        end
        xmin *= 2

        maxiters -= 1
        if maxiters == 0 error("Max iterations reached (xmin doubling)") end
    end

    while true
        y = f(xmax)
        if y > y0 && abs(y-y0) > yeps
            break
        end
        xmax *= 2

        maxiters -= 1
        if maxiters == 0 error("Max iterations reached (xmax doubling)") end
    end

    xmin,xmax
end

end # module
