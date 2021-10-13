using LinearAlgebra: diag, norm, Symmetric
using JuMP
import Gurobi

function compute_distances(coordinates)
    n = size(coordinates, 1)
    distances = zeros(n, n)
    for i=2:n
        for j=1:i-1
            distances[i, j] = floor(norm(coordinates[i, :] - coordinates[j, :]))
        end
    end
    return Symmetric(distances, :L)
end


function sommets_visites(x,d,f)
    path = [d]
    next_airport = d
    while next_airport != f
        next_airport = argmax(x[next_airport, :])
        append!(path, next_airport)
    end
    return path
end
