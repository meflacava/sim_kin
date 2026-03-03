# sim_kin
Code to simulate hatchery spawning and calculate kinship trends for Delta Smelt

## Definitions:

Kinship coefficient = probability that two alleles, one from individual i and one from j, are identical by descent. Calculated between 2 potential parents
Relatedness = 2x kinship coefficient to get at how related 2 fish are to each other. Calculated between 2 potential parents
Inbreeding coefficient = Probability that the two alleles within individual i are IBD. Calculated for individuals. Equivalent to parents' kinship coefficient

Example: 2 parents that are siblings
 - relatedness between parents = 0.5
 - kinship coefficient between parents = 0.25
 - inbreeding coefficient for their offspring = 0.25


Relatedness  Kinship   Pair type
1            0.5       Self
0.5          0.25      Full sib
0.25         0.125     Half sib
0.125        0.0625    Cousins
0.0625       0.03125   Second cousins
