### Data Parser ###

# Data Read: Reads the structured data from the given instance file
# Input:  - filename: A string containing the file name with its relative path.
# Output: - n: An Integer containing the number of tasks;
#         - tasks_times: An (n) Array of Integers with each task's time;
#         - prec_relations: A (n,n) 2-dimensional Binary Array representing
#                           each task's Precedence Relation.
function data_read(filename::String)
    # Gets every line in the instance file for processing
    lines = []
    open(filename,"r") do f
        lines = readlines(filename)
    end
    
    # The number n of tasks is given in the first line
    n = parse(Int64, lines[1])
    
    # Gets each task's time in the next n lines
    task_times = Int64[]
    for i = 2:(n+1)
        append!(task_times, parse(Int64, lines[i]))
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
    return n, task_times, prec_relations
end

# Parser Write: Parses the Data from the given Instance to fit GLPK pattern.
# Input:  - filename: A string with the name for the parser file that will be written;
#         - n, s, task_times, times_sum, prec_relations: Needed to write the Data.
# Output: None.
function parser_write(filename::String, n, s, task_times, prec_relations)
    open(filename, "w") do f
        print(f, "data;\r\n")
        print(f, "\r\n")
        print(f, "param n := ", n, ";\r\n")
        print(f, "param m := ", s, ";\r\n")
        print(f, "param times := ")
        for i = 1:n
            print(f, i, " ", task_times[i], "\r\n")
        end
        print(f, ";\r\n")  
        print(f, "precedences : ")
        for i = 1:n
            print(f, i, " ")
        end
        print(f, ":=\r\n")
        for i = 1:n
            prec_line = string(prec_relations[i, :])
            prec_line = replace(prec_line, ",", "")
            prec_line = replace(prec_line, "[", "")
            prec_line = replace(prec_line, "]", "")
            print(f, i, " ", prec_line, "\r\n")
        end      
        print(f, ";\r\n")
    end
    
    println("Solution saved in file: ", filename)
    println("")
end

### Main Function ###
function main()
    s = 3
    
    filename = ["instances/HAHN.IN2", "instances/LUTZ3.IN2", "instances/WEE-MAG.IN2"]
    out_file = ["instances/HAHN.DAT", "instances/LUTZ3.DAT", "instances/WEE-MAG.DAT"]
    
    for i = 1:(length(filename))
        n, task_times, prec_relations = data_read(filename[i])
        parser_write(out_file[i], n, s, task_times, prec_relations)
    end
    
    println("Parsing finished.")
end

main()
