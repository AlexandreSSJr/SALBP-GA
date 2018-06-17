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
    for i = 1:n
        curr_time = parse(Int64, lines[(i+1)])
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
            println(f, "No feasible solution found.")
            
            println("- FINAL RESULT: No feasible solution found.")
        else
            println(f, solution)
            println(f, time)
            
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
function args_read(args = ["Hahn", "1", "100"])    
    f = lowercase(args[1])
    
    if (length(args) == 1)
        println("WARNING: Only one argument provided, using default stations and population sizes of 1 and 100, respectively")
        println("")
        stations = 1
        population = 100
    elseif (length(args) == 2)
        println("WARNING: Population size not provided, using default of 100")
        println("")
        stations = parse(Int64, args[2])
        population = 100
    else
        stations = parse(Int64, args[2])
        population = parse(Int64, args[3])
    end
    
    # Defines the name of the instance file based on the first argument information
    if (f == "hahn" || f == "1")
        file = "instances/HAHN.IN2"
    elseif (f == "lutz3" || f == "2")
        file = "instances/LUTZ3.IN2"
    elseif (f == "wee-mag" || f == "3")
        file = "instances/WEE-MAG.IN2"
    else
        println("WARNING: Instance name not recognized, using default instance Hahn")
        println("")
        file = "instances/HAHN.IN2"
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
    
    return gene
end

# Rand Population: 
# Input:  - size: An Integer representing the size of the generated population;
#         - s, n, task_times, times_sum: Needed from gene and fitness functions.
# Output: - population: 
function rand_population(size::Int64, s::Int64, n::Int64, task_times, times_sum::Int64)
    population = []
    for i = 1:size
        gene = rand_gene(s, n)
        fit = fitness(gene, task_times, times_sum)
        curr_gene = [fit, gene]
        append!(population, [curr_gene])
    end
    
    return population
end

# Crossover: Creates a new Gene by doing a crossover between two existing ones
# Input:  - gene1: A (s,n) 2-dimensional Integer Array representing the first gene
#           gene2: A (s,n) 2-dimensional Integer Array representing the second gene
# Output: - new_gene: A (s,n) 2-dimensional Integer Array representing the new gene
function crossover(gene1, gene2)
    #TODO
    return gene1
end

# TODO: Evaluate use of sum and dict to remove both fors
# Fitness: Calculates the fitness of the given gene
#          In SALBP, it's the longest cycle among all stations
#          A cycle of a station is the sum the times of all the tasks given to it
# Input:  - gene: A (s,n) 2-dimensional Integer Array representing a gene;
#         - task_times: An (n) Array of Integers with each task's time.
#         - times_sum: An Integer with the sum of all tasks' times to use as
#                      "punishment" for each infeasibility in the gene.
# Output: - max_cycle: An Integer representing the fitness of the given gene.
function fitness(gene, task_times, times_sum::Int64)
    max_cycle = 0
    
    for i = 1:(size(gene)[1])
        current_cycle = 0
        for j = 1:(size(gene)[2])
            if (gene[i, j] > 0)
                current_cycle += task_times[gene[i, j]]
            end
        end
        if (current_cycle > max_cycle)
            max_cycle = current_cycle
        end
    end
        
    return max_cycle
end

# Genetic Algorithm: 
# Input:  - init_pop: An array containing the initial [fitness, gene] pairs population
#                     of the previously defined size;
#         - pop_size, s, n, task_times, times_sum: Needed from gene, population
#                                                  and fitness functions;
#         - max_it: An integer containing the stop criteria of the GA, max iterations.
# Output: - best_solution: The [fitness, gene] pair of the best solution found;
#         - time_taken: The time it took to find the best solution found;
function genetic_algorithm(init_pop, pop_size::Int64, s::Int64, n::Int64, task_times, times_sum::Int64, max_it::Int64)
    tic()
    population = init_pop
    best_found_fit = times_sum + 1
    best_solution = [0, 0]
    time_taken = 0
        
    for i = 1:max_it
        #sort(population, lt=lexless)
        curr_best = population[1][1]
        
        if (curr_best <= times_sum)
            if (curr_best < best_found_fit)
                best_found_fit = curr_best
                best_solution = population[1]
                time_taken = toc()
            end
        end
        
        population = rand_population(pop_size, s, n, task_times, times_sum)
    end
        
    return best_solution, time_taken
end

### Main Function ###
function main()
    # TODO: Possibly move max iterations to parameters
    max_it = 10
    out_filename = "main_test.txt"
    
    println("- Reading arguments...")
    println("")
    file_name, num_stations, population_size = args_read(ARGS)
    
    println("- Setting up data from given file...")
    println("")
    num_tasks, task_times, pred_relations, times_sum = data_read(file_name)
    
    println("- Starting Genetic Algorithm execution...")
    println("")
    initial_population = rand_population(population_size, num_stations, num_tasks, task_times, times_sum)
    best_solution, time_taken = genetic_algorithm(initial_population, population_size, num_stations, num_tasks, task_times, times_sum, max_it)
        
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
