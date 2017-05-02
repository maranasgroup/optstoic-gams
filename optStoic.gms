$Ontext
optStoic

Authors: Anupam Chowdhury, Chiam Yu Ng

DESCRIPTION
The objective of this code (optStoic.gms) is to identify the overall stoichiometry of
reactants and products. The code identifies the maximum achievable yields of the
target product subject to thermodynamic restrictions.

In this sample code, we use glucose as substrate and acetate/1,3-propanediol as product.

OUTPUT FILE
optStoic_output.csv - a comma-separated list of model status, metabolites and stoichiometric
                        coefficients, the standard Gibbs free energy of reaction.

For more information, please refer to the paper below:
Anupam Chowdhury and Costas D. Maranas
"Designing overall stoichiometric conversions and intervening metabolic reactions"
Scientific Reports, 2015
$Offtext

$INLINECOM /*  */
$set fpath data/
$set fsuffix _modified
$set substrate C00267

Options
    limrow = 1000
    optCR = 1E-9
    optCA = 0.0
    iterlim = 1000000
    decimals = 8
    reslim = 1000000
    work = 50000000
	mip = cplex
	minlp = baron;

***********************
Sets

i       metabolites
/
$include "%fpath%metabolites%fsuffix%.txt"
/

j       "elements, charge and dGf"
/
'C'
'H'
'N'
'O'
'P'
'S'
'F'
'Cl'
'Mg'
'Fe'
'Se'
'Co'
'As'
'Br'
'I'
'R'
'charge'
'dGf'
/

elem(j) subset of j
/
'C'
'H'
'O'
'charge'
/

dG(j)   change in Gibbs free energy of formation / 'dGf'/

pdt(i)  list of target products
/
'C00033'    /*acetate*/
'C02457'	/*13pdo*/
/

allow(i) metabolites allowed to participate in a reaction
/
'C00007'    /*o2*/
'C00267'    /*glc*/
'C00011'    /*co2*/
'C00080'    /*h+*/
'C00001'    /*h2o*/
'C02457'    /*13pdo*/
'C00033'    /*acetate*/
/
;

***********************
Parameters

m(i,j)  metabolite details
/
$include "%fpath%metabolite_details.txt"
/

c       "standard Gibbs free energy of a reaction" /0/

dGmax   "maximum allowed value of free energy change" /-5/
;

***********************
Variables
zp      objective function
s(i)    stoichiometric coefficient of metabolite i
;

***********************
Equations

obj     define objective function
stoic   elemental and charge balance
bal     the overall conversion should have a free energy change lesser or equal to dGmax
;

obj.. 			        zp =e= sum(i$pdt(i),s(i))/(-s('%substrate%'));
stoic(j)$(elem(j))..   sum(i,s(i)*m(i,j)) =e= 0;
bal(j)$(dG(j))..	    sum(i,s(i)*m(i,j)) =l= dGmax*(-s('%substrate%'));

***********************
Model optstoic
/
obj
stoic
bal
/
;


file f1 /optStoic_output.csv/;
f1.pc=5;
put f1;

*CSV file header
put 'product', 'modelstat', 'obj_value', 'dGr_overall', 'stoichiometric_coefficient', 'metabolite'/;


********************************************************************************
*                       Example 1: glucose -> acetate                          *
********************************************************************************

*set the upper and lower bound of stoichoimetric coefficient
s.lo(i) = -15;
s.up(i) = 15;

*set stoichiometric coefficient glucose to -1 (reactant)
s.up('%substrate%') = -1;

*allow water to be either co-reactant or co-product
s.up('C00001') = 15;
s.lo('C00001') = -15;

*turn off oxygen
s.up('C00007') = 0;
s.lo('C00007') = 0;

*turn off other product
s.fx(i)$pdt(i)= 0;

*set lb of target product's stoichiometry coefficient to 1
s.lo('C00033') = 1; /*acetate*/

*set ub of target product's stoichiometry coefficient to 15
s.up('C00033') = 15; /*acetate*/

*set stoichiometric coefficient of metabolites not allowed to participate in final reaction to 0
s.fx(i)$(not allow(i)) = 0;

optstoic.optfile = 1;
optstoic.holdfixed = 1;

Solve optstoic using minlp maximizing zp;

*get dG of the entire reaction
c = sum(i, s.l(i)*m(i,'dGf'))/(-s.l('%substrate%'));

*output result
put 'C00033',optstoic.modelstat:0:0,zp.l:0:4,c:0:4,;

loop(i$(s.l(i) ne 0),
    put s.l(i):0:4,i.tl,;
);
put /;

if(optstoic.modelstat ne 1,
    put optstoic.modelstat:0:0,"infeasible"/;
);


********************************************************************************
*                       Example 2: glucose -> 1,3-propanediol                  *
********************************************************************************

*set the upper and lower bound of stoichoimetric coefficient
s.lo(i) = -15;
s.up(i) = 15;

*set ub of stoichiometric coefficient of substrate(e.g. glucose) to -1 (reactant)
s.up('%substrate%') = -1;

*make oxygen a reactant
s.up('C00007') = 0;

*turn off other products
s.fx(i)$pdt(i) = 0;

*set lb of target product's stoichiometry coefficient to 1
s.lo('C02457') = 1;

*set ub of target product's stoichiometry coefficient to 15
s.up('C02457') = 15;

*set stoichiometric coefficient of metabolites not allowed to participate in final reaction to 0
s.fx(i)$(not allow(i)) = 0;

optstoic.optfile = 1;
optstoic.holdfixed = 1;

Solve optstoic using minlp maximizing zp;

*get dG of the entire reaction
c = sum(i, s.l(i)*m(i,'dGf'))/(-s.l('%substrate%'));

*output result
put 'C02457', optstoic.modelstat:0:0,zp.l:0:4,c:0:4,;
loop(i$(s.l(i) ne 0),
    put s.l(i):0:4,i.tl,;
);
put /;

if(optstoic.modelstat ne 1,
    put optstoic.modelstat:0:0,"infeasible"/;
);
