mol representation NewCartoon 0.300000 6.000000 4.100000 0
mol color ColorID 0
mol selection {residue 0 1 2 3 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 }
mol material Opaque
mol addrep top
mol representation NewCartoon 0.300000 6.000000 4.100000 0
mol color ColorID 1
mol selection {residue 4 }
mol material Opaque
mol addrep top
graphics top color black
mol representation VDW 1.000000 8.000000
mol color ColorID 1
mol selection {residue 4 and name CA P}
mol material Opaque
mol addrep top
mol representation VDW 1.000000 8.000000
mol color ColorID 0
mol selection {residue 23 and name CA P}
mol material Opaque
mol addrep top
set sel1 [atomselect top "residue 4 and name CA P"]
set sel2 [atomselect top "residue 23 and name CA P"]
set coord1 [lindex [$sel1 get {x y z}] 0]
set coord2 [lindex [$sel2 get {x y z}] 0]
graphics top cylinder $coord1 $coord2 radius 1.000000
$sel1 delete
$sel2 delete
mol representation VDW 1.000000 8.000000
mol color ColorID 1
mol selection {residue 4 and name CA P}
mol material Opaque
mol addrep top
mol representation VDW 1.000000 8.000000
mol color ColorID 0
mol selection {residue 27 and name CA P}
mol material Opaque
mol addrep top
set sel1 [atomselect top "residue 4 and name CA P"]
set sel2 [atomselect top "residue 27 and name CA P"]
set coord1 [lindex [$sel1 get {x y z}] 0]
set coord2 [lindex [$sel2 get {x y z}] 0]
graphics top cylinder $coord1 $coord2 radius 1.000000
$sel1 delete
$sel2 delete
