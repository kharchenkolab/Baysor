# TODO: replace it with existing triangulation polygons
convex_hull(a::Matrix{<: Real}) = convex_hull(collect.(eachcol(a)))
function convex_hull(a::Array{Vector{T},1} where T<: Real)
    cw(a, b, c) = (a[1]*(b[2]-c[2])+b[1]*(c[2]-a[2])+c[1]*(a[2]-b[2]) < 0);
    ccw(a, b, c) = (a[1]*(b[2]-c[2])+b[1]*(c[2]-a[2])+c[1]*(a[2]-b[2]) > 0);

    if length(a) == 1
        return a
    end

    a = sort(a, lt = (a,b) -> (a[1] < b[1] || a[1] == b[1] && a[2] < b[2]))
    p1 = a[1];
    p2 = a[end];
    up = [p1]
    down = [p1]

    for i in 2:length(a)
        if (i==length(a) || cw(p1, a[i], p2))
            while (length(up) >= 2 && !cw(up[end-1], up[end], a[i]))
                up = up[1:end-1]
            end

            push!(up, a[i])
        end
        if (i==length(a) || ccw(p1, a[i], p2))
            while (length(down) >= 2 && !ccw(down[end-1], down[end], a[i]))
                down = down[1:end-1]
            end

            push!(down, a[i])
        end
    end

    return hcat(vcat(up, reverse(down[1:end-1]))...);
end

function area(polygon::Matrix{<:Real})
    area = 0.0;

    for i in 1:size(polygon, 2)
        j = i % size(polygon, 2) + 1;
        area += polygon[1, i] * polygon[2, j] - polygon[2, i] * polygon[1, j];
    end

   return abs(area) / 2;
end
