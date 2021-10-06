using LinearAlgebra: diag, norm, Symmetric
using JuMP
import Gurobi

function compute_distances(coordinates)::Symmetric{Float64, Array{Float64, 2}}
    n = size(coordinates, 1)
    distances = zeros(n, n)
    for i=2:n
        for j=1:i-1
            distances[i, j] = floor(norm(coordinates[i, :] - coordinates[j, :]))
        end
    end
    return Symmetric(distances, :L)
end