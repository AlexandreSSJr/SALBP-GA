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
for i = 1:n
    println(i, rs[i,:])
end
    
#solution_write("test.txt", 36, 20132.12302)

fs, ms, ps = args_read(["3"])
println(fs)
println(ms)
println(ps)
