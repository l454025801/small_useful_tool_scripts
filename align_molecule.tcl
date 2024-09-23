  package require Orient
  namespace import Orient::orient

  set all [atomselect top "all"]
  set I [draw principalaxes $all]
  set A [orient $all [lindex $I 2] {0 0 1}]
  $all move $A
  set I [draw principalaxes $all]
  set A [orient $all [lindex $I 1] {0 1 0}]
  $all move $A
  set I [draw principalaxes $all]

