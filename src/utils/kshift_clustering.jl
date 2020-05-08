# Copyright (c) 2015: Rene Donner.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of
# this software and associated documentation files (the "Software"), to deal in the Software
# without restriction, including without limitation the rights to use, copy, modify, merge,
# publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons
# to whom the Software is furnished to do so, subject to the following conditions:
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
# NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
# WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
# OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

using NearestNeighbors
using StatsBase: denserank

# DEPRECATED
function kshiftmedoids(data, k)
    centers = kshifts(data,k)
    tree = KDTree(float(data))
    ids = vec(mapslices(a -> knn(tree, a, 1)[1][1], centers, dims=1))

    return data[:,ids], ids
end

# DEPRECATED
kshifts(data::Array{T,2} where T, k::Int) = kshifts!(deepcopy(data[:,rand(1:size(data, 2), k)]), data)

# DEPRECATED
@inbounds function kshifts!(centers::Matrix{T}, data::Matrix{T}; f1::Float64=29/30) where T
    @assert size(data, 1) == size(centers, 1) == 2
    f2 = 1 - f1
    for i = 1:size(data, 2)
        c_data = data[:,i]
        v = typemax(eltype(centers))
        ind = 0
        for j = 1:size(centers, 2)
            d1 = (c_data[1]-centers[1,j])
            d1 *= d1

            if d1 > v
                continue
            end

            d2 = (c_data[2]-centers[2,j])
            d = d1 + d2 * d2

            if d < v
                v = d
                ind = j
            end
        end

        centers[:,ind] .= f1 .* centers[:,ind] .+ f2 .* c_data
    end

    return centers
end

kshiftlabels(data::Matrix{T}, centers::Matrix{T}) where T = vcat(knn(KDTree(centers), data, 1)[1]...)
