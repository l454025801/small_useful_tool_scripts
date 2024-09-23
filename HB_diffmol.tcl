#draw_hbonds: draw cylinders between hydrogen and acceptor atoms in hydrogen bonds formed between two selections. 
#Cutoff and angle are also specified.
#1. source HB_diffmol.tcl
#2. draw_hhbonds 3.5 30 {chain A} {chain B}

proc draw_hhbonds { cutoff angle sel1text sel2text} {
  set a [atomselect top $sel1text]
  set b [atomselect top $sel2text]
  set coords [[atomselect top all] get {x y z}]
  
  # get hbonds with 1 as donor and 2 as acceptor
  foreach { d1 a1 h1 } [measure hbonds $cutoff $angle $a $b] break
  # get hbonds wiith 2 as donor and 1 as acceptor
  foreach { d2 a2 h2 } [measure hbonds $cutoff $angle $b $a] break
  # combine the lists
  set hyd [concat $h1 $h2]
  set acc [concat $a1 $a2]
  # draw them!
  draw delete all
  foreach hi $hyd ai $acc {
    set hpos [lindex $coords $hi]
    set apos [lindex $coords $ai]
    draw color green
    draw cylinder $hpos $apos radius 0.3
  }
}

proc draw_hhbridges { cutoff angle sel1text sel2text sel3text} {
  set a [atomselect top $sel1text]
  set b [atomselect top $sel2text]
  set c [atomselect top $sel3text]
  set coords [[atomselect top all] get {x y z}]

  # get hbonds with 1 as donor and 2 as acceptor
  foreach { d1 a1 h1 } [measure hbonds $cutoff $angle $a $b] break
  # get hbonds wiith 2 as donor and 1 as acceptor
  foreach { d2 a2 h2 } [measure hbonds $cutoff $angle $b $a] break

  # get hbonds with 2 as donor and 3 as acceptor
  foreach { d1 a3 h3 } [measure hbonds $cutoff $angle $b $c] break
  # get hbonds wiith 3 as donor and 2 as acceptor
  foreach { d2 a4 h4 } [measure hbonds $cutoff $angle $c $b] break

  # combine the lists
  set hyd1 [concat $h1 $h2]
  set acc1 [concat $a1 $a2]

  set hyd2 [concat $h3 $h4]
  set acc2 [concat $a3 $a4]

  # take the overlap 
  set hyd [intersect $hyd1 $hyd2]
  set acc [intersect $acc1 $acc2]

  # draw them!
  draw delete all
  foreach hi $hyd ai $acc {
    set hpos [lindex $coords $hi]
    set apos [lindex $coords $ai]
    draw color green
    draw cylinder $hpos $apos radius 0.3
  }
}

