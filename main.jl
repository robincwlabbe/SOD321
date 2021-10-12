#import des packages

using LinearAlgebra: diag, norm, Symmetric
using JuMP
import Gurobi


include("read.jl")
include("functions.jl")


#function resolution(path,mod)
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

n,d,f,Amin,Nr,R,regions,coords = readInstance("C:/Users/User/Desktop/3A/M2/SOD321/SOD321/instance_70_1.txt")
mod="polynomial"
	#n,d,f,Amin,Nr,R,regions,coords = readInstance(path)	
	D = compute_distances(coords)

	#declaration du modele
	model = Model(Gurobi.Optimizer)
	
	#declaration des variables
	@variable(model, x[1:n,1:n], Bin)  #xij= 1 if the pilot travels from airport i to airport j during the flight, 0 otherwise
	

	#declaration de l objectif
	@objective(model, Min, sum(x .* D))

	#declaration des contraintes
	@constraint(model, diag(x) .== 0)  # removes x_ii
	
	#On décolle au plus une fois dans chaque aeroport
	@constraint(model, [i in 1:n], sum(x[i,j] for j in 1:n)<=1)
	
	#On  au plus une fois dans chaque aeroport
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
		JuMP.optimize!(model)
				
		println("fin")
		println(JuMP.objective_value(model))
	else
		JuMP.optimize!(model)
	#definition variable pour boucle while
		global obj_value = 1
	
		while obj_value >0
			#creation sous probleme
			sub_problem = Model(with_optimizer(Gurobi.Optimizer))
			@variable(sub_problem, a[1:n], Bin)  #ai= 1 si l'aeroport i est visité
			y = model[:x]
			@objective(sub_problem, Max, sum(sum(JuMP.value(y[i,j])*a[i]*a[j] for j in 1:n) for i in 1:n ) - sum(a[i] for i in 1:n) + 1)
			@constraint(sub_problem, sum(a[i] for i in 1:n)>=1)
			#resolution sous probleme
			JuMP.optimize!(sub_problem)
			#ajout de la contrainte qui viole le plus et resolution du sous probleme
			@constraint(model, sum(sum(x[i,j]*JuMP.value(a[i])*JuMP.value(a[j]) for j in 1:n) for i in 1:n ) <= sum(JuMP.value(a[i]) for i in 1:n)-1)
			JuMP.optimize!(model)
			#mise a jour pour voir si l'objectif est positif
			global obj_value = JuMP.objective_value(sub_problem)

		end

		println("fin")
		println(JuMP.objective_value(model))


	end

	#return(model)
#end

#### Separation ####

#On optimise le modele maitre



	#declaration des variables
#    
#	#affichage des resultats
#obj_value = JuMP.objective_value(model)
#println("Objective value: ", obj_value)
#	value.(x)
#end


