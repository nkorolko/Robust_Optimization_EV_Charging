Robust Optimization of EV Charging Schedules in Unregulated Electricity Markets

Nikita Korolko & Zafer Sahinoglu

In this project, we address the problem of optimal electric vehicle charging in an unregulated electricity market. This problem is known to be highly nonlinear even in the case of fixed electricity prices due to a nonlinear state-of-charge curve representing physical battery limitations. We design tractable formulations for single and multiple EV charging frameworks. In the first part of the paper, we develop a new efficient cutting plane method, that can be used for solving charging optimization problem for both scenarios of known and uncertain electricity prices. The latter scenario with real-time electricity rates is considered in the second part of the paper. We obtain robust optimization counterparts of the nominal charging problems that are particularly important from an economic perspective when budget constraints are strictly enforced. 

New robust formulations are proven to be tractable. Moreover, computational experiments
illustrate that a decision maker can find solutions that are close to optimal in terms of the corresponding objective values, and
robust with respect to uncertain electricity prices.

File "Article_Robust.pdf" is a paper published in IEEE Transactions on Smart Grid. 

Julia file "MINLP_1.jl" solves a mixed-integer nonlinear optimization problem (9) on page 7. The solving algorithm uses lazy constraints technique described in Section 2 of the paper and also calls Gurobi. This algorithm is more than two orders of magnitude faster than off-the-shelf nonlinear solvers like Bonmin. 

File "Figure3.m" is an example of visualization for the paper (Fig.3 on page 6) coded in MATLAB. 

[1] http://www.gurobi.com/

[2] http://julialang.org/

[3] https://projects.coin-or.org/Bonmin

