clear UT

## X TRANSLATION ##

# 2mm scale
setscale 2
setoption smoothing 2
setoption boundguess 2
setoption paramsubset 1  0 0 0 1 0 0 0 0 0 0 0 0
clear U
clear UA
setrow UA 1 0 0 0  0 1 0 0  0 0 1 0  0 0 0 1
optimise 12 UA:1   0.0   0.0   0.0    0.0   0.0   0.0   0.0   abs 4
optimise 12 UA:1   0.0   0.0   0.0    2.0   0.0   0.0   0.0   abs 4
optimise 12 UA:1   0.0   0.0   0.0   -2.0   0.0   0.0   0.0   abs 4
optimise 12 UA:1   0.0   0.0   0.0    4.0   0.0   0.0   0.0   abs 4
optimise 12 UA:1   0.0   0.0   0.0   -4.0   0.0   0.0   0.0   abs 4
clear UA
copy U UA

# 1mm scale
setscale 1
setoption smoothing 1
setoption boundguess 1
setoption paramsubset 1  0 0 0 1 0 0 0 0 0 0 0 0
clear UB
clear UL
clear UM
# remeasure costs at this scale
clear U
measurecost 12 UA 0 0 0 0 0 0 abs
sort U
copy U UL
# optimise best 3 candidates
clear U
optimise 12 UL:1-3  0.0   0.0   0.0   0.0   0.0   0.0   0.0  abs 4
# also try the identity transform as a starting point at this resolution
clear UQ
setrow UQ  1 0 0 0  0 1 0 0  0 0 1 0  0 0 0 1
optimise 12 UQ  0.0   0.0   0.0   0.0   0.0   0.0   0.0  rel 4
clear UB
copy U UB

# 0.5mm scale
setscale 0.5
setoption smoothing 0.5
setoption boundguess 0.5
setoption paramsubset 1  0 0 0 1 0 0 0 0 0 0 0 0
clear U
clear UC
clear UF
# remeasure costs at this scale
measurecost 12 UB 0 0 0 0 0 0 abs
sort U
copy U UC
clear U
optimise 12  UC:1  0.0   0.0   0.0   0.0   0.0   0.0   0.0  abs 4
copy U UF


## Y TRANSLATION ##

# 2mm scale
setscale 2
setoption smoothing 2
setoption boundguess 2
setoption paramsubset 1  0 0 0 0 1 0 0 0 0 0 0 0
clear U
clear UA
setrow UA 1 0 0 0  0 1 0 0  0 0 1 0  0 0 0 1
optimise 12 UA:1   0.0   0.0   0.0   0.0    0.0   0.0   0.0   abs 4
optimise 12 UA:1   0.0   0.0   0.0   0.0    2.0   0.0   0.0   abs 4
optimise 12 UA:1   0.0   0.0   0.0   0.0   -2.0   0.0   0.0   abs 4
optimise 12 UA:1   0.0   0.0   0.0   0.0    4.0   0.0   0.0   abs 4
optimise 12 UA:1   0.0   0.0   0.0   0.0   -4.0   0.0   0.0   abs 4
clear UA
copy U UA

# 1mm scale
setscale 1
setoption smoothing 1
setoption boundguess 1
setoption paramsubset 1  0 0 0 0 1 0 0 0 0 0 0 0
clear UB
clear UL
clear UM
# remeasure costs at this scale
clear U
measurecost 12 UA 0 0 0 0 0 0 abs
sort U
copy U UL
# optimise best 3 candidates
clear U
optimise 12 UL:1-3  0.0   0.0   0.0   0.0   0.0   0.0   0.0  abs 4
# also try the identity transform as a starting point at this resolution
clear UQ
setrow UQ  1 0 0 0  0 1 0 0  0 0 1 0  0 0 0 1
optimise 7 UQ  0.0   0.0   0.0   0.0   0.0   0.0   0.0  rel 4
clear UB
copy U UB

# 0.5mm scale
setscale 0.5
setoption smoothing 0.5
setoption boundguess 0.5
setoption paramsubset 1  0 0 0 0 1 0 0 0 0 0 0 0
clear U
clear UC
clear UF
# remeasure costs at this scale
measurecost 12 UB 0 0 0 0 0 0 abs
sort U
copy U UC
clear U
optimise 12  UC:1  0.0   0.0   0.0   0.0   0.0   0.0   0.0  abs 4
copy U UF

## sort the 2 results to pick the best
clear U
copy UT U
sort U

# now do a general 2 DOF translation to refine this
clear UA
copy U UA

# 2mm scale
setscale 2
setoption smoothing 2
setoption paramsubset 2  0 0 0 1 0 0 0 0 0 0 0 0  0 0 0 0 1 0 0 0 0 0 0 0
clear U
optimise 12 UA:1  0.0   0.0   0.0   0.0   0.0   0.0   0.0  rel 4

# 1mm scale
setscale 1
setoption smoothing 1
setoption paramsubset 2  0 0 0 1 0 0 0 0 0 0 0 0  0 0 0 0 1 0 0 0 0 0 0 0
clear UB
clear UL
clear UM
# remeasure costs at this scale
clear U
measurecost 12 UA 0 0 0 0 0 0 rel
sort U
copy U UL
# optimise best 3 candidates
clear U
optimise 12 UL:1-3  0.0   0.0   0.0   0.0   0.0   0.0   0.0  rel 4
# also try the identity transform as a starting point at this resolution
clear UQ
setrow UQ  1 0 0 0  0 1 0 0  0 0 1 0  0 0 0 1
optimise 7 UQ  0.0   0.0   0.0   0.0   0.0   0.0   0.0  rel 4
copy U UB

# 0.5mm scale
setscale 0.5
setoption smoothing 0.5
setoption paramsubset 2  0 0 0 1 0 0 0 0 0 0 0 0  0 0 0 0 1 0 0 0 0 0 0 0
clear U
clear UC
clear UD
clear UE
clear UF
# remeasure costs at this scale
measurecost 12 UB 0 0 0 0 0 0 rel
sort U
copy U UC
clear U
optimise 12  UC:1  0.0   0.0   0.0   0.0   0.0   0.0   0.0  rel 4
sort U
copy U UF

