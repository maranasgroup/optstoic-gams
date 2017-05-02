# optstoic-gams
optStoic-minFlux/minRxn GAMS code. If you use this code in your work, please cite the following paper: Anupam Chowdhury and Costas D. Maranas,
"Designing overall stoichiometric conversions and intervening metabolic reactions", Scientific Reports, 2015 (PMID: 26530953).

Requirement:
===========
1. GAMS (https://www.gams.com/)
2. BARON solver
3. CPLEX solver


Description:
===========
1. optStoic 
- An MINLP-based optimization procedure for identifying the optimal stoichiometry of a reaction subject to thermodynamic constraints, charge and elemental balance. 

2. minFlux or minRxn
- An MILP procedure to identify pathway(s) that carry minimal flux or consist of minimal number of reaction that conform to the given reaction stoichiometry.

The reaction and metabolite IDs are derived from the KEGG database (http://www.genome.jp/kegg/).  
