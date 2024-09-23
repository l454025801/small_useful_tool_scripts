  package require Orient
  namespace import Orient::orient

  set sel [atomselect top "name CB CE2"]
  set all [atomselect top "all"]
  set I [draw principalaxes $sel]
  set A [orient $sel [lindex $I 2] {0 0 1}]
  $all move $A
  set I [draw principalaxes $sel]
  set A [orient $sel [lindex $I 1] {0 1 0}]
  $all move $A
  set I [draw principalaxes $sel]
  set begin [atomselect top "name CB"]
  set coor [lindex [$begin get {x y z}] 0]
  $all moveby [vecscale -1 $coor]
