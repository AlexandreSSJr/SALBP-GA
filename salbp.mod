# Simple Assembly Line Balancing Problem - Type 2
#
# The problem consists of finding an assignment of tasks
# to an arrengement of linear stations where these tasks
# are executed in a particular order that respects its
# precedence relations. The type 2 SALB Problem
# implemented here aims to minimize the maximum
# cycle time among the stations for a given
# number of stations and tasks' times.
#
# References: Marcus Ritt and Alysson M. Costa. 
#             A comparison of formulations for the 
#             simple assembly line balancing problem. 
#             Int. Trans. Oper. Res., 2015.

# The set of Tasks
param n;
set TASKS := 1..n;

# The set of Stations
param m;
set STATIONS:= 1..m;

# Each task's execution time
param times {i in TASKS};

# The tasks precedence relation binary matrix
param precedences {i in TASKS, p in TASKS};

# The cycle to be minimized
var c, integer >= 0;

# A variable X to denote if the Task i is assigned to the Station j
var X {j in STATIONS, i in TASKS}, binary;

# Objective function to minimize cycle time
minimize cycle_time: c;

# Constraint to allow a Task to be assigned only once in the solution
s.t. assigned_once {i in TASKS}: sum{j in STATIONS} x[j][i] = 1;

# Constraint to assure the precedence relations are respected
s.t. respect_precedences {i in TASKS, p in TASKS, s in STATIONS: precedences[i,p] = 1}: X[s,p] <= sum{j in STATIONS: j <= s} X[j,i];

# Sets the final cycle to the maximum sum of times among stations
s.t. set_cycle {j in STATIONS}: sum{i in TASKS} times[i] * x[j][i] <= c;

solve;

display_obj;

end;
