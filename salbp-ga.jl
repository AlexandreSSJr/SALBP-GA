### General Utility Functions ###

# Data Read: Reads the structured data from the given instance file
# Input:  - filename: A string containing the file name with its relative path.
# Output: - n: An Integer containing the number of tasks;
#         - tasks_times: An (n) Array of Integers with each task's time;
#         - prec_relations: A (n,n) 2-dimensional Binary Array representing
#                           each task's Precedence Relation.
#         - times_sum: An Integer representing the sum of all the tasks' times.
function data_read(filename::String)
    # Gets every line in the instance file for processing
    lines = []
    open(filename,"r") do f
        lines = readlines(filename)
    end
    
    # The number n of tasks is given in the first line
    n = parse(Int64, lines[1])
    
    # Gets each task's time in the next n lines and the sum of all times
    task_times = Int64[]
    times_sum = 0
    for i = 2:(n+1)
        curr_time = parse(Int64, lines[i])
        append!(task_times, curr_time)
        times_sum += curr_time
    end
    
    # Gets the i,j Direct Precedence Relations until the "-1,-1" end mark
    # and transforms it into a Binary Matrix of the Relations
    prec_relations = zeros(Int64, (n,n))
    curr_line = n+2
    # Julia has no "Do While" so a first verification is necessary to avoid a break
    rs = split(lines[curr_line], ",")
    r = [parse(Int64, rs[1]), parse(Int64, rs[2])]
    while (r[1] > 0)
        prec_relations[r[1], r[2]] = 1
        
        curr_line += 1        
        rs = split(lines[curr_line], ",")
        r = [parse(Int64, rs[1]), parse(Int64, rs[2])]
    end
    
    println("Data successfully read from file \"", filename, "\" containing ", n, " tasks.")
    println("")
    
    # Returns the structured data
    return n, task_times, prec_relations, times_sum
end

# Solution Write: Writes the final solution and time in a file with stdout format
# Input:  - filename: A string with the name for the file that will be written;
#         - solution: The final best solution encountered during execution;
#         - time: The time it took to find the given solution.
# Output: None.
function solution_write(filename::String, solution, time)
    open(filename, "w") do f
        if (time == 0)
            println(f, "No feasible solution found.\n")
            
            println("- FINAL RESULT: No feasible solution found.")
        else
            println(f, "Best solution: ", solution[2], "\r\n")
            println(f, "With score: ", solution[1], "\r\n")
            println(f, "Found in: ", time, " seconds.\r\n")
            
            println("- FINAL RESULT: ")
            println("- Best solution: ", solution[2])
            println("- With score: ", solution[1])
            println("- Found in: ", time, " seconds.")
        end
    end
    
    println("Solution saved in file: ", filename)
    println("")
end

# Args Read: Reads any (or provides default) Argments passed through command line
# Input:  - args: An Array of Strings containing the argments from the command line.
# Output: - file: A string with the name of the file for the instance given;
#         - stations: An Integer with the number of stations available, defaults to 1;
#         - population: An Integer with the size of the population, defaults to 100.
function args_read(args = ["instances/HAHN.IN2", "1", "500"]) 
    if (length(args) == 0)
        println("WARNING: No argument provided, using default instance Hahn and stations and population sizes of 1 and 100, respectively")
        file = "instances/HAHN.IN2"
        stations = 1
        population = 500
    elseif (length(args) == 1)
        println("WARNING: Only one argument provided, using default stations and population sizes of 1 and 100, respectively")
        file = args[1]
        stations = 1
        population = 500
    elseif (length(args) == 2)
        println("WARNING: Population size not provided, using default of 500")
        file = args[1]
        stations = parse(Int64, args[2])
        population = 500
    else
        file = args[1]
        stations = parse(Int64, args[2])
        population = parse(Int64, args[3])
    end
    
    println("Arguments successfully read. Using these parameters:")
    println("File Path and Name: \"", file, "\".")
    println("Number of Stations: ", stations, ".")
    println("Size of Population: ", population, ".")
    println("")
    
    return file, stations, population
end

### Genetic Algorithm Functions ###

# Rand Gene: Generates a random gene containing a solution for the problem
# Input:  - s: An Integer representing the given number of stations available;
#         - n: An Integer representing the number of tasks to be assigned.
# Output: - gene: A (s,n) 2-dimensional Integer Array representing the sequence
#                 of tasks n assigned to each station s.
function rand_gene(s::Int64, n::Int64)
    gene = zeros(Int64, (s,n))
    task_counter = ones(Int64, s)
    
    for i = 1:n
        rnds = rand(1:s)
        gene[rnds, task_counter[rnds]] = i
        task_counter[rnds] += 1
    end
    
    # 0.1% chance of just generating a linear gene
    rnd_linear = rand(1:1000)
    if (rnd_linear == 500)
        gene = zeros(Int64, (s,n))
        task_counter = 1
        threshold = ceil(Int64, (n/s))
        
        for i = 1:s
            for j = 1:threshold
                if (task_counter <= n)
                    gene[i, j] = task_counter
                    task_counter += 1
                end
            end
        end
    end
    
    return gene
end

# Rand Population: 
# Input:  - size: An Integer representing the size of the generated population;
#         - s, n, task_times, times_sum, prec_relations: Needed from gene and fitness functions.
# Output: - population: A random population containing "size"s [fitness, gene] pairs.
function rand_population(size::Int64, s::Int64, n::Int64, task_times, times_sum::Int64, prec_relations)
    population = []
    for i = 1:size
        gene = rand_gene(s, n)
        fit = fitness(gene, task_times, times_sum, prec_relations)
        curr_individual = [fit, gene]
        append!(population, [curr_individual])
    end
    
    return population
end

# Crossover: Creates a new Gene by doing a crossover between two parents
# Input:  - parent_1: A [fitness, gene] individual representing the first parent;
#         - parent_2: A [fitness, gene] individual representing the second parent.
#         - task_times, times_sum, prec_relations: Needed from fitness function.
# Output: - new_individual: A new [fitness, gene] individual.
function crossover(parent_1, parent_2, task_times, times_sum::Int64, prec_relations)
    gene1 = parent_1[2]
    gene2 = parent_2[2]
    
    s = size(gene1)[1]
    n = size(gene1)[2]
    
    genes = [gene1, gene2]
    task_counter = ones(Int64, s)
    new_gene = zeros(Int64, (s,n))
    
    for i = 1:n
        rnd_parent = rand(1:2)
        for j = 1:s
            discover_task = find(genes[rnd_parent][j, :] .== i)
            if (!isempty(discover_task))
                new_gene[j, task_counter[j]] = i
                task_counter[j] += 1
            end
        end
    end
    
    fit = fitness(new_gene, task_times, times_sum, prec_relations)
    new_individual = [fit, new_gene]
    
    return new_individual
end

# Fitness: Calculates the fitness of the given gene
#          In SALBP, it's the longest cycle among all stations
#          A cycle of a station is the sum the times of all the tasks given to it
# Input:  - gene: A (s,n) 2-dimensional Integer Array representing a gene;
#         - task_times: An (n) Array of Integers with each task's time.
#         - times_sum: An Integer with the sum of all tasks' times to use as
#                      "punishment" for each infeasibility in the gene.
#         - prec_relations: A (n,n) 2-dimensional Binary Array representing
#                           each task's Precedence Relation.
# Output: - max_cycle: An Integer representing the fitness of the given gene.
function fitness(gene, task_times, times_sum::Int64, prec_relations)
    max_cycle = 0
    s = size(gene)[1]
    n = size(gene)[2]
    tasks_not_exec = ones(Int64, n)
    
    # For each station
    for i = 1:s
        current_cycle = 0
        # And each task assigned to it
        for j = 1:n
            curr_task_id = gene[i, j]
            if (curr_task_id > 0)
                # Adds the current task time to the current cycle total
                current_cycle += task_times[curr_task_id]
                # This is for penalities addition on infeasibilities
                for t = 1:n
                    # If a task has not yet been executed
                    if (tasks_not_exec[t] == 1)
                        # And it's in the precedence list for the current task
                        if(prec_relations[t, curr_task_id] == 1)
                            # Adds the penalty
                            current_cycle += times_sum
                        end
                    end
                end
                # Removes the current task from the not executed yet list
                tasks_not_exec[curr_task_id] = 0
            end
        end
        # Keeps the maximum cycle found
        if (current_cycle > max_cycle)
            max_cycle = current_cycle
        end
    end
        
    return max_cycle
end

# Genetic Algorithm: 
# Input:  - init_pop: An array containing the initial [fitness, gene] pairs population
#                     of the previously defined size;
#         - pop_size, s, n, task_times, times_sum, prec_relations: Needed from gene, population
#                                                                  and fitness functions;
#         - max_gen: An integer containing the stop criteria of the GA, max generations.
# Output: - best_solution: The [fitness, gene] pair of the best solution found;
#         - time_taken: The time it took to find the best solution found;
function genetic_algorithm(init_pop, pop_size::Int64, s::Int64, n::Int64, task_times, times_sum::Int64, max_gen::Int64, prec_relations)
    tic()
    population = init_pop
    best_found_fit = times_sum + 1
    best_solution = [0, 0]
    time_taken = 0
    
    println("Genetic Algorithm started.")
        
    # Until the end criteria is hit
    for i = 1:max_gen
        population = sort(population, lt=lexless)
        curr_best = population[1][1]
        
        if (curr_best <= times_sum)
            if (curr_best < best_found_fit)
                best_found_fit = curr_best
                best_solution = population[1]
                time_taken += toc()
                tic()
                println("IT [", i, "] - ", "Found a better solution with score: ", best_solution[1])
            end
        end
        
        # Creates a new population for the next generation
        class_A = ceil(Int64, pop_size * 0.2)
        class_C = ceil(Int64, pop_size * 0.5)
        
        # Class A stays the same, changes start in Class B
        for i = class_A:class_C-1
            # Class B crossover
            curr_parents = i - class_A + 1
            population[i] = crossover(population[curr_parents], population[curr_parents + 1], task_times, times_sum, prec_relations)
        end
        # Class C just replaces the worst half of the population for new random genes
        population[(class_C+1):end,:] = rand_population(class_C, s, n, task_times, times_sum, prec_relations)
    end
        
    return best_solution, time_taken
end

### Main Function ###
function main()
    max_gen = 500
    out_filename = "main_test.txt"
    
    println("- Reading arguments...")
    println("")
    file_name, num_stations, population_size = args_read(ARGS)
    
    # TEST #
    #file_name = "instances/TEST1.IN2"
    #file_name = "instances/HAHN.IN2"
    #num_stations = 10
    #population_size = 100
    # TEST #
    
    println("- Setting up data from given file...")
    println("")
    num_tasks, task_times, prec_relations, times_sum = data_read(file_name)
    
    println("- Starting Genetic Algorithm execution...")
    println("")
    initial_population = rand_population(population_size, num_stations, num_tasks, task_times, times_sum, prec_relations)
    best_solution, time_taken = genetic_algorithm(initial_population, population_size, num_stations, num_tasks, task_times, times_sum, max_gen, prec_relations)
        
    solution_write(out_filename, best_solution, time_taken)
    
    println("- Execution ended.")
    println("")
end

main()

########## TESTS ##########

#=
file1 = "instances/HAHN.IN2"
file2 = "instances/LUTZ3.IN2"
file3 = "instances/WEE-MAG.IN2"

ns , ts, rs, bigm = data_read(file1)
println(ns)
println(ts)
for i = 1:ns
    println(i, rs[i,:])
end
println(bigm)

fs, ms, ps = args_read(["3"])

genes0 = rand_gene(ms, ns)

fs, ms, ps = args_read(["1", "2"])

genes = rand_gene(ms, ns)
println(genes)

testfits = [0, 0]
for i = 1:size(genes)[2]
    if (genes[1, i] > 0)
        testfits[1] += ts[genes[1, i]]
    end
    if (genes[2, i] > 0)
        testfits[2] += ts[genes[2, i]]
    end
end
println(testfits)

fits = fitness(genes, ts, bigm)
println(fits)

new_genes = crossover(genes0, genes)
println(new_genes)
=#
