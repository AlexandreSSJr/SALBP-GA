# This is just an experiment with linear only tasks #

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
    gene = Int[]
    task_counter = ones(Int64, s)
    
    if (s > 1)
        cut = rand(1:n)
        for i = 1:(s-1)
            append!(gene, cut)
            cut = rand(cut:n)
        end
    end
    
    return gene
end

# Rand Population: 
# Input:  - size: An Integer representing the size of the generated population;
#         - s, n, task_times: Needed from gene and fitness functions.
# Output: - population: A random population containing "size"s [fitness, gene] pairs.
function rand_population(size::Int64, s::Int64, n::Int64, task_times)
    population = []
    for i = 1:size
        gene = rand_gene(s, n)
        fit = fitness(s, n, gene, task_times)
        curr_individual = [fit, gene]
        append!(population, [curr_individual])
    end
    
    return population
end

# Crossover: Creates a new Gene by doing a crossover between two parents
# Input:  - parent_1: A [fitness, gene] individual representing the first parent;
#         - parent_2: A [fitness, gene] individual representing the second parent.
#         - s, n, task_times: Needed from fitness function.
# Output: - new_individual: A new [fitness, gene] individual.
function crossover(s, n, parent1, parent2, task_times)
    gene1 = parent1[2]
    gene2 = parent2[2]
    
    new_gene = gene1
    rnd = rand(1:2)
    
    if (rnd == 1)
        for i = 1:(s-1)
            new_gene[i] = min(new_gene[i], gene2[i])
        end
    else
        for i = 1:(s-1)
            it = s - i
            new_gene[it] = max(new_gene[it], gene2[it])
        end
    end
    
    fit = fitness(s, n, new_gene, task_times)
    new_individual = [fit, new_gene]
    
    return new_individual
end

# Fitness: Calculates the fitness of the given gene
#          In SALBP, it's the longest cycle among all stations
#          A cycle of a station is the sum the times of all the tasks given to it
# Input:  - gene: A (s,n) 2-dimensional Integer Array representing a gene;
#         - task_times: An (n) Array of Integers with each task's time.
#         - s, n: Needed from gene function.
# Output: - max_cycle: An Integer representing the fitness of the given gene.
function fitness(s, n, gene, task_times)
    max_cycle = 0
    
    if (s == 1)
        max_cycle = sum(task_times)
    elseif (s == 2)
        max_cycle = max(task_times[1:gene[1]], task_times[(gene[1]+1):n])
    else
        curr_cut = 0
        for i = 1:(s-1)
            if(curr_cut != gene[i])
                curr_cut += 1
                curr_cycle = sum(task_times[curr_cut:gene[i]])
                curr_cut = gene[i]
                max_cycle = max(max_cycle, curr_cycle)
            end
        end
        curr_cut += 1
        if (curr_cut < n)
            curr_cycle = sum(task_times[curr_cut:n])
            max_cycle = max(max_cycle, curr_cycle)
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
        
        # Now to create a new population for the next generation
        class_A = ceil(Int64, pop_size * 0.2)
        class_B = ceil(Int64, pop_size * 0.3)
        class_C = ceil(Int64, pop_size * 0.5)
        
        # Class A stays the same, changes start in Class B
        for i = class_A:class_C-1
            # Class B crossover
            curr_parents = i - class_A + 1
            population[i] = crossover(s, n, population[curr_parents], population[curr_parents + 1], task_times)
        end
        # Class C just replaces the worst half of the population for new random genes
        population[(class_C+1):end,:] = rand_population(class_C, s, n, task_times)
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
    
    println("- Setting up data from given file...")
    println("")
    num_tasks, task_times, prec_relations, times_sum = data_read(file_name)
    
    println("- Starting Genetic Algorithm execution...")
    println("")
    initial_population = rand_population(population_size, num_stations, num_tasks, task_times)
    best_solution, time_taken = genetic_algorithm(initial_population, population_size, num_stations, num_tasks, task_times, times_sum, max_gen, prec_relations)
        
    solution_write(out_filename, best_solution, time_taken)
    
    println("- Execution ended.")
    println("")
end

main()
