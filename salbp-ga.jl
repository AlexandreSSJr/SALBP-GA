### General Utility Functions ###

# Data Read: Reads the structured data from the given instance file
# Input:  - filename: A string containing the file name with its relative path.
# Output: - n: An Integer containing the number of tasks;
#         - tasks_times: An (n) Array of Integers with each task's time;
#         - prec_relations: A (n,n) 2-dimensional Binary Array representing
#                           each task's Precedence Relation.
function data_read(filename)
    # Gets every line in the instance file for processing
    lines = []
    open(filename,"r") do f
        lines = readlines(filename)
    end
    
    # The number n of tasks is given in the first line
    n = parse(Int64, lines[1])
    
    # Gets each task's time in the next n lines
    task_times = Int64[]
    for i = 1:n
        append!(task_times, parse(Int64, lines[(i+1)]))
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
    
    # Returns the structured data
    return n, task_times, prec_relations
end

# Solution Write: Writes the final solution and time in a file with stdout format
# Input:  - filename: A string with the name for the file that will be written;
#         - solution: The final best solution encountered during execution;
#         - time: The time it took to find the given solution.
# Output: None.
function solution_write(filename, solution, time)
    open(filename, "w") do f
        println(f, solution)
        println(f, time)
    end
end

# Args Read: Reads any (or provides default) Argments passed through command line
# Input:  - args: An Array of Strings containing the argments from the command line.
# Output: - file: A string with the name of the file for the instance given;
#         - stations: An Integer with the number of stations available, defaults to 1;
#         - population: An Integer with the size of the population, defaults to 100.
function args_read(args = ["Hahn", "1", "100"])    
    f = lowercase(args[1])
    
    if (length(args) == 1)
        warn("Only one argument provided, using default stations and population sizes of 1 and 100, respectively")
        stations = 1
        population = 100
    elseif (length(args) == 2)
        warn("Population size not provided, using default of 100")
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
        warn("Instance name not recognized, using default instance Hahn")
        file = "instances/HAHN.IN2"
    end
    
    return file, stations, population
end

### Genetic Algorithm Functions ###

# TODO: Evaluate if this generated solution should already be feasible
#       This changes not only generation but also crossover and mutation
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

# TODO: If there are infeasible solution genes this should return -1
#       Also, evaluate use of sum and dict to remove both fors
# Fitness: Calculates the fitness of the given gene
#          In SALBP, it's the longest cycle among all stations
#          A cycle of a station is the sum the times of all the tasks given to it
# Input:  - gene: A (s,n) 2-dimensional Integer Array representing a gene;
#         - task_times: An (n) Array of Integers with each task's time.
# Output: - max_cycle: An Integer representing the fitness of the given gene.
function fitness(gene, task_times)
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

### Main Function ###
function main()
    fname, m, pop = args_read(ARGS)
    n, times, relations = data_read(fname)
    
    #TODO
    
    #solution_write(out_fname, best_solution, time_taken)
end

main()

########## TESTS ##########

file1 = "instances/HAHN.IN2"
file2 = "instances/LUTZ3.IN2"
file3 = "instances/WEE-MAG.IN2"

ns , ts, rs = data_read(file1)
println(ns)
println(ts)
for i = 1:ns
    println(i, rs[i,:])
end
    
#solution_write("test.txt", 36, 20132.12302)

fs, ms, ps = args_read(["3"])
println(fs)
println(ms)
println(ps)

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

fits = fitness(genes, ts)
println(fits)
