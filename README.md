# CASTANET (Computational package for deriving Analytical STAtionary distribution of biochemical reaction networks with NEtwork Translation)

This is a matlab code to analytically derive stationary distributions for a given stochastic biochemical reaction networks. Detailed step-by-step manual will be uploaded soon. 

## Code Description
1. CRN_main.m
> Main function for this package. Users need to work only with this function. If one specify the source complexes, product complexes, and propensity functions of reactions then the code provides a symbolic expression for a stationary distribution. 

Since all the below functions are automatically run by CRN_main.m function, users do not need to run nor edit the below functions separately.

2. CRN_translation.m
> This function generates all possible translated networks whose reaction orders are at most 'max_order', which is specified by users.

3. CRN_find_elementary_path.m
> This function finds an elementary path. Here, the elementary path means that the sequence of a pair of the source complexes, which form a basis for the lattice.

4. CRN_find_elementary_function.m
> This function finds elementary functions corresponding to the elementary basis, identified by the function CRN_find_elementary_path.m. The elementary function is the term to be multiplied when the state moves from n-e_j to n. Here, e_j is the jth element of the elementary basis, which corresponds to the jth elementary function.

5. CRN_solve_sym_linear.m
> This functions compute the coordinate of n-n_0 where n_0 is the start point with respect to the elementary basis.

6. CRN_check_factorization_condition.m
> This function tests whether the construced candidate theta_c satisfies the factorization condition.

7. CRN_compute_cbe.m
> This functions compute a complex balanced equilibrium of the deterministic mass action model for the translated network. 


## Limitations
1. Solve the system of linear equation with symbolic variables.
 - When the CRN_solve_sym_linear.m is performed, the code may not fully use underlying assumptions for the variables. It sometimes leads to failure of the code. For instance, if one sovle the system of linera equations x = 5 - a and y = a - 5 under the assumption x == -y, then the solution is x=5-a itself because the assumption x==-y is consistent to the system of equations. However, the code cannot solve this system. This occurs especially when there is a conservation law in a given biochemical reaction network.

