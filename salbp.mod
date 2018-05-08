# Simple Assembly Line Balancing Problem
#
# TODO: PROBLEM DESCRIPTION
#
# References: Marcus Ritt and Alysson M. Costa. 
#             A comparison of formulations for the 
#             simple assembly line balancing problem. 
#             Int. Trans. Oper. Res., 2015.

set TASKS;
set STATIONS;

param times {i in TASKS};
param precedences {i in TASKS, j in STATIONS};

var c, integer >= 0;
var X {i in TASKS, j in STATIONS}, binary;

minimize cycle_time: c;

s.t. assigned_once {i in TASKS}: sum{j in STATIONS} x[i][j] = 1;

# TODO
# s.t. respect_precedences :

s.t. set_cycle {j in STATIONS}: sum{i in TASKS} ttimes[i] * x[i][j] <= c;