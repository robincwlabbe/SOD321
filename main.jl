
#import des packages
using LinearAlgebra: diag, norm, Symmetric, Tables, CSV
using JuMP
import Gurobi
include("read.jl")
include("functions.jl")


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
function resolution(path,mod)
	n,d,f,Amin,Nr,R,regions,coords = readInstance(path)
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
	if (d == f)
		indicatrice_depart_arrivee = 1
	else
		indicatrice_depart_arrivee = 0
	end

	@constraint(model, min_visits, sum(x) >= Amin -1 + indicatrice_depart_arrivee)

	solving_start = time()

	# Subtour
	if (mod == "polynomial")
		@variable(model, u[i = 1:n], Int)
		@constraint(model, subtour[i=1:n, j=1:n],
			u[j] >= u[i] + 1 - n * (1 - x[i, j]))
		JuMP.optimize!(model)
	elseif (mod == "exponential_1") #resolution SEP
	#### Separation ####
	#### on enleve le caractere binaires ####
		for v in all_variables(model)
				if is_binary(v)
				    unset_binary(v)
				    lb = has_lower_bound(v) ? lower_bound(v) : -Inf
				    ub = has_upper_bound(v) ? upper_bound(v) : Inf
				    set_lower_bound(v, max(0.0, lb))
				    set_upper_bound(v, min(1.0, ub))
				end
		end
		isBinary = false
		set_silent(model)
		isOptimal = false
		while !(isOptimal)

			optimize!(model)
			println("Current objective value : ",objective_value(model))

			#creation sous probleme
			oldstd = stdout
			redirect_stdout(open("null", "w")) #rediriger temporairement l'output
			sub_problem = Model(Gurobi.Optimizer)
			redirect_stdout(oldstd) #rediriger l'output vers sa sortie normale

			set_silent(sub_problem)

			@variable(sub_problem, a[1:n], Bin)  #ai= 1 si l'aeroport i est visité
			@objective(sub_problem, Max, sum(sum(JuMP.value(x[i,j])*a[i]*a[j] for j in 1:n) for i in 1:n ) - sum(a[i] for i in 1:n) + 1)
			@constraint(sub_problem, sum(a[i] for i in 1:n)>=1)
			#resolution sous probleme
			JuMP.optimize!(sub_problem)

			if (objective_value(sub_problem)>0)
				println("Violated subtour constraint found.")
				#ajout de la contrainte qui viole le plus et resolution du sous probleme
				@constraint(model, sum(sum(x[i,j]*JuMP.value(a[i])*JuMP.value(a[j]) for j in 1:n) for i in 1:n )<=sum(JuMP.value(a[i]) for i in 1:n)-1)
			else
				#println("No violated subtour constraints.")
				#isOptimal = true
				if !(isBinary)
					println("No violated subtour constraints for root relaxation.")
					isOptimal = false
					isBinary = true
					for v in all_variables(model)
						set_binary(v)
					end
				else
					println("No violated subtour constraints.")
					isOptimal=true
				end
								
			end
		end
	elseif (mod == "exponential_2") #resolution SEP generalisé
		#### Separation ####
		### ajout de la variable y ###
		@variable(model, y[1:n], Bin)
	    @constraint(model, [i in filter(x->!(x in [f]), 1:n)], sum(x[i,:]) == y[i])
	    @constraint(model, y[f] == 1)
	    @constraint(model, y[d] == 1)

	    for v in all_variables(model)
				if is_binary(v)
				    unset_binary(v)
				    lb = has_lower_bound(v) ? lower_bound(v) : -Inf
				    ub = has_upper_bound(v) ? upper_bound(v) : Inf
				    set_lower_bound(v, max(0.0, lb))
				    set_upper_bound(v, min(1.0, ub))
				end
		end
		isBinary = false

	    set_silent(model)
		isOptimal = false
		while !(isOptimal)
			optimize!(model)
			println("Current objective value : ",objective_value(model))

			#creation sous probleme
			oldstd = stdout
			redirect_stdout(open("null", "w")) #rediriger temporairement l'output
			sub_problem = Model(Gurobi.Optimizer)
			redirect_stdout(oldstd) #rediriger l'output vers sa sortie normale

			set_silent(sub_problem)

			@variable(sub_problem, a[1:n], Bin)  #ai= 1 si l'aeroport i est visité
			@variable(sub_problem, h[1:n], Bin)

			#objectifs et contraintes du sous probleme
			@objective(sub_problem, Max, sum(sum(JuMP.value(x[i,j])*a[i]*a[j] for j in 1:n) for i in 1:n )
			- sum(JuMP.value(y[i])*a[i] for i in 1:n)
			+ sum(JuMP.value(y[i])*h[i] for i in 1:n))

			@constraint(sub_problem, sum(h[i] for i in 1:n)==1)
			@constraint(sub_problem, [i in 1:n],h[i]<=a[i])

			optimize!(sub_problem)

			if (objective_value(sub_problem)>0)
				println("Violated subtour constraint found.")
				#ajout de la contrainte qui viole le plus et resolution du sous probleme
				
				S = BitArray(value.(a))
	        	i0 = BitArray(value.(h))
	        	@constraint(model, sum(x[S,S]) <= sum(y[S]) - sum(y[i0]))


	        	#=
				@constraint(model, sum(x[i,j]*value(a[i])*value(a[j]) for i in 1:n for j in 1:n)
				<= sum(sum(x[i,j] for j in 1:n)*(value(a[i]) - value(h[i]))
				for i in filter(e->!(e in [d,f]),1:n)) + value(a[d]) + value(a[f])
				- value(h[d]) - value(h[f]))
				=#
				
			else
				if !(isBinary)
					println("No violated subtour constraints for root relaxation.")
					isOptimal = false
					isBinary = true
					for v in all_variables(model)
						set_binary(v)
					end
				else
					println("No violated subtour constraints.")
					isOptimal=true
				end
			end
		end
	end

	println("="^50)
	println("Optimal value : ", objective_value(model))
	println("Total resolution time : ", time() - solving_start)
	return(value.(x),objective_value(model),time() - solving_start)
end
