#import des packages

using LinearAlgebra: diag, norm, Symmetric
using JuMP
import Gurobi


include("read.jl")
include("functions.jl")


#function resolution(n,d,f,Amin,Nr,R,regions,coords,mod)
# ---------------------------------------------------------------
# Function that compiutes the optimum
# Parameter meanings :
# n : number of airports
# d : "depart"/start
# f : "fin"/end
# Amin : minimum number of airports to go through
# Nr : number of regions (all must be visited !)
# R : maximum range of a plane
# regions : array such that regions[region number] =
# list of airports in that region
# coords : coordinates of all airports
# mod : exponential or polynomial formulation to solve the problem.
# ---------------------------------------------------------------

n,d,f,Amin,Nr,R,regions,coords=readInstance("/Users/antoine/Desktop/3A_ENSTA/SOD321_ELLOUMI/PROJET/Instances-20211005/instance_6_1.txt")
mod="polynomial"

	D = compute_distances(coords)

	#declaration du modele
	model = Model(Gurobi.Optimizer)
	
	#declaration des variables
	@variable(model, x[1:n,1:n], Bin)  #xij= 1 if the pilot travels from airport i to airport j during the flight, 0 otherwise
	

	#declaration de l objectif
	@objective(model, Min, sum(x .* D))

	#declaration des contraintes
	@constraint(model, diag(x) .== 0)  # removes x_ii
	
	#On n’atterit au plus une fois dans chaque aeroport
	@constraint(model, [i in 1:n], sum(x[i,j] for j in 1:n)<=1)
	
	#On décolle au plus une fois dans chaque aeroport
	@constraint(model, [i in 1:n], sum(x[j,i] for j in 1:n)<=1) 

	#contrainte sur le depart
	@constraint(model, sum(x[d, :]) == 1)

	#contrainte sur l'arrivée
	@constraint(model, sum(x[:, f]) == 1)
	
	#si départ différent d'arrivée, on ne peut arriver au départ, ou partir de l'arrivée
	if d != f
		@constraint(model, sum(x[:, d]) == 0)
		@constraint(model, sum(x[f, :]) == 0)
	end

	# contrainte sur la distance que peut parcourir l'avion
    @constraint(model, x[D .> R] .== 0) 
   	
    #On décolle et atterit d'un aeroport si celui ci n'est pas le départ ou l'arrivée
   	depart_arrive = []
	for i in 1:n
    	if i != d && i != f 
    		append!(depart_arrive,i)
    	end
    end   	

    # on decolle et on atterit des aeroports visités en dehors du depart et de l'arrivee
    @constraint(model, [i in depart_arrive],sum(x[i,j] - x[j,i] for j in 1:n) == 0) 

    #toutes les regions sont visitées:

    reg = [i for i in 1:Nr]
    @constraint(model, [i in reg], sum(sum(x[k,j] + x[j,k] for j in 1:n) for k in regions[i]) >= 1)

	#nombre minimum d'aeroports à visiter   
	if d == f
		indicatrice_depart_arrivee = 1
	end
	if d !=f
		indicatrice_depart_arrivee = 0
	end
	@constraint(model, min_visits, sum(x) >= Amin -1 + indicatrice_depart_arrivee) 

	# Subtour
	if mod == "polynomial"
		@variable(model, u[i = 1:n], Int)
		@constraint(model, subtour[i=1:n, j=1:n],
					u[j] >= u[i] + 1 - n * (1 - x[i, j]))
	end

	#return(model)
#end


#function optimization(model)
#	JuMP.optimize!(model)
#    
#	#affichage des resultats
#	obj_value = JuMP.objective_value(model)
#	println("Objective value: ", obj_value)
#	value.(x)
#end


