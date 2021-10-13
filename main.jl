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

n,d,f,Amin,Nr,R,regions,coords = readInstance("/Users/antoine/Desktop/3A_ENSTA/SOD321_ELLOUMI/PROJET/Instances-20211005/instance_40_1.txt")
mod = "exponential_1"
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
		global solving_time = MOI.get(model, MOI.SolveTime())
		println("aaa")
	end
	
	if mod == "exponential_1" #resolution SEP
	#### Separation ####
		#JuMP.unset_binary(x::Matrix{VariableRef})

		#On optimise le modele maitre
		JuMP.optimize!(model)
		global solving_time = MOI.get(model, MOI.SolveTime())
		#definition variable pour boucle while
		sub_problem = Model(with_optimizer(Gurobi.Optimizer))
		@variable(sub_problem, a[1:n], Bin)  #ai= 1 si l'aeroport i est visité
		y = model[:x]
		@objective(sub_problem, Max, sum(sum(JuMP.value(y[i,j])*a[i]*a[j] for j in 1:n) for i in 1:n ) - sum(a[i] for i in 1:n) + 1)
		@constraint(sub_problem, sum(a[i] for i in 1:n)>=1)
		#resolution sous probleme
		JuMP.optimize!(sub_problem)
		#mise a jour pour voir si l'objectif est positif
		global obj_value = JuMP.objective_value(sub_problem)
		global solving_time = solving_time + MOI.get(sub_problem, MOI.SolveTime())
		
		if obj_value>0
			#ajout de la contrainte qui viole le plus et resolution du sous probleme
			@constraint(model, sum(sum(x[i,j]*JuMP.value(a[i])*JuMP.value(a[j]) for j in 1:n) for i in 1:n ) <= sum(JuMP.value(a[i]) for i in 1:n)-1)
			JuMP.optimize!(model)
			global solving_time = solving_time + MOI.get(model, MOI.SolveTime())
		end

		while obj_value >0 
			#creation sous probleme
			oldstd = stdout
			redirect_stdout(open("null", "w")) #rediriger temporairement l'output
			sub_problem = Model(Gurobi.Optimizer)
			
			@variable(sub_problem, a[1:n], Bin)  #ai= 1 si l'aeroport i est visité
			y = model[:x]
			@objective(sub_problem, Max, sum(sum(JuMP.value(y[i,j])*a[i]*a[j] for j in 1:n) for i in 1:n ) - sum(a[i] for i in 1:n) + 1)
			@constraint(sub_problem, sum(a[i] for i in 1:n)>=1)
			#resolution sous probleme
			JuMP.optimize!(sub_problem)
			global solving_time = solving_time + MOI.get(sub_problem, MOI.SolveTime())
			#mise a jour pour voir si l'objectif est positif
			global obj_value = JuMP.objective_value(sub_problem)

			if obj_value>0
				#ajout de la contrainte qui viole le plus et resolution du sous probleme
				@constraint(model, sum(sum(x[i,j]*JuMP.value(a[i])*JuMP.value(a[j]) for j in 1:n) for i in 1:n ) <= sum(JuMP.value(a[i]) for i in 1:n)-1)
				JuMP.optimize!(model)
				global solving_time = solving_time + MOI.get(model, MOI.SolveTime())
			end
			redirect_stdout(oldstd) #rediriger l'output vers sa sortie normale
		end
		#JuMP.unset_binary(x::Matrix{VariableRef})
		#JuMP.optimize!(model)
		#global solving_time = solving_time + MOI.get(model, MOI.SolveTime())
		JuMP.optimize!(model)
		global solving_time = solving_time + MOI.get(model, MOI.SolveTime())
		println("Temps de résolution : ", solving_time)
		println(JuMP.objective_value(model))
		println("fin1")
	end


	if mod == "exponential_2" #resolution SEP generalisé
	#### Separation ####
		#JuMP.unset_binary(x::Matrix{VariableRef})

		#On optimise le modele maitre
		JuMP.optimize!(model)
		global solving_time = MOI.get(model, MOI.SolveTime())

		#creation du sous probleme
		sub_problem = Model(with_optimizer(Gurobi.Optimizer))
		@variable(sub_problem, a[1:n], Bin)  #ai= 1 si l'aeroport i est visité
		@variable(sub_problem, h[1:n], Bin)
		
		#definition des variables selections d'arcs et sommets pour la solutiondu modele
		y = model[:x]
		sommets_visites = []
		for i in 1:n
			for j in 1:n
				if JuMP.value(x[i,j])!=0
					append!(sommets_visites,[i,j])
				end	
			end
		end

		x_tilde = [0 for i in 1:n]
		for i in sommets_visites
			x_tilde[i] = 1
		end
		#objectifs et contraintes du sous probleme
		@objective(sub_problem, Max, sum(sum(JuMP.value(y[i,j])*a[i]*a[j] for j in 1:n) for i in 1:n ) - sum(JuMP.value(x_tilde[i])*a[i] for i in 1:n) + sum(JuMP.value(x_tilde[i])*h[i] for i in 1:n))
		@constraint(sub_problem, sum(h[i] for i in 1:n)==1)
		@constraint(sub_problem, [i in 1:n],h[i]<=a[i])

		#resolution sous probleme
		JuMP.optimize!(sub_problem)
		global solving_time = solving_time + MOI.get(sub_problem, MOI.SolveTime())

		#mise a jour pour voir si l'objectif est positif
		global obj_value = JuMP.objective_value(sub_problem)


		if obj_value>0
			S = BitArray(value.(a))
        	i0 = BitArray(value.(h))
        	@constraint(model, sum(x[S,S]) <= sum(x_tilde[S]) - sum(x_tilde[i0]))
        	JuMP.optimize!(model)
			global solving_time = solving_time + MOI.get(model, MOI.SolveTime())
		end

		while obj_value >0 && solving_time<60

			#creation sous probleme
			oldstd = stdout
			redirect_stdout(open("null", "w")) #rediriger temporairement l'output
			
			
			sub_problem = Model(with_optimizer(Gurobi.Optimizer))
			@variable(sub_problem, a[1:n], Bin)  #ai= 1 si l'aeroport i est visité
			@variable(sub_problem, h[1:n], Bin)
			
			y = model[:x]
			sommets_visites = []
			for i in 1:n
				for j in 1:n
					if JuMP.value(x[i,j])!=0
						append!(sommets_visites,[i,j])
					end	
				end
			end
			x_tilde = [0 for i in 1:n]
			for i in sommets_visites
				x_tilde[i] = 1
			end

			@objective(sub_problem, Max, sum(sum(JuMP.value(y[i,j])*a[i]*a[j] for j in 1:n) for i in 1:n ) - sum(JuMP.value(x_tilde[i])*a[i] for i in 1:n) + sum(JuMP.value(x_tilde[i])*h[i] for i in 1:n))
			@constraint(sub_problem, sum(h[i] for i in 1:n)==1)
			@constraint(sub_problem, [i in 1:n],h[i]<=a[i])
			#resolution sous probleme
			JuMP.optimize!(sub_problem)
			global solving_time = solving_time + MOI.get(sub_problem, MOI.SolveTime())
			
			#mise a jour pour voir si l'objectif est positif
			global obj_value = JuMP.objective_value(sub_problem)

			if obj_value>0
				#ajout de la contrainte qui viole le plus et resolution du sous probleme
				S = BitArray(value.(a))
	        	i0 = BitArray(value.(h))
	        	@constraint(model, sum(x[S,S]) <= sum(x_tilde[S]) - sum(x_tilde[i0]))

				#@constraint(model, sum(sum(x[i,j]*JuMP.value(a[i])*JuMP.value(a[j]) for j in 1:n) for i in 1:n ) <= sum(JuMP.value(a[i]) for i in 1:n)-sum(JuMP.value(a[i])*JuMP.value(h[i]) for i in 1:n))
				JuMP.optimize!(model)
				global solving_time = solving_time + MOI.get(model, MOI.SolveTime())
			end
			if floor(solving_time)% 30 ==0
				println("temps passé : ",solving_time)
			end
			redirect_stdout(oldstd) #rediriger l'output vers sa sortie normale
		end
		#JuMP.unset_binary(x::Matrix{VariableRef})
		#JuMP.optimize!(model)
		#global solving_time = solving_time + MOI.get(model, MOI.SolveTime())
		JuMP.optimize!(model)
		global solving_time = solving_time + MOI.get(model, MOI.SolveTime())
		println("Temps de résolution : ", solving_time)
		println(JuMP.objective_value(model))
		println("fin2")
	end




