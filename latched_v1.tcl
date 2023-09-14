mol addfile "../p00-noWater.parm7" type parm7 
mol addfile "DISTANCES_filtered_2to6nm.nc" type netcdf waitfor all

set outputFile1 "output_intersection_points_and_penetration"
set file_out1 $outputFile1.txt
set outfile1 [open $file_out1 w]
set outputFile2 "output_prelasso_vs_pretadpole"
set file_out2 $outputFile2.txt
set outfile2 [open $file_out2 w]

################################ PROCEDURES ####################################

proc triNorm {comLoop com1 com2} {	
	set cP [veccross [vecsub $com1 $comLoop] [vecsub $com2 $comLoop]] 
	return $cP
}

proc findIntersection {pointOnLine pointOnPlane lineDirection planeNorm} {
    set denominator [vecdot $lineDirection $planeNorm]
    set numerator [vecdot [vecsub $pointOnPlane $pointOnLine] $planeNorm]
    if {$denominator == 0} {
        return -1
    } else {
        set d [expr {double($numerator) / $denominator}]
    }
   if {$d > 1 || $d < 0} { ;# Line segment : 0 <= d <= 1
       return -2
   } else {
       set intersection [vecadd $pointOnLine [vecscale $d $lineDirection]]
       return $intersection
   } 
}

proc barycentricTransform {comLoop comRes1 comRes2 intersectionPoint} {
    set w0 [vecsub $comRes2 $comLoop]
    set w1 [vecsub $comRes1 $comLoop]
    set w2 [vecsub $intersectionPoint $comLoop]

    set dot00 [vecdot $w0 $w0]
    set dot01 [vecdot $w0 $w1]
    set dot02 [vecdot $w0 $w2]
    set dot11 [vecdot $w1 $w1]
    set dot12 [vecdot $w1 $w2]

    set denom [expr ($dot00 * $dot11) - ($dot01 * $dot01)]
    set u [expr (($dot11 * $dot02) - ($dot01 * $dot12)) / $denom]
    set v [expr (($dot00 * $dot12) - ($dot01 * $dot02)) / $denom]
    set sumuv [expr $u + $v]
    set insideTriangleBinary [expr ($u >= 0) && ($v >= 0) && ($sumuv <= 1)]
    return $insideTriangleBinary ;# 0 for no, 1 for yes
}

proc vmd_draw_arrow {mol start end} { ;# start = comLoop for the normal vectors
    draw color red
    draw sphere $start radius 0.5 resolution 500
    draw color purple
	set middle [vecadd $start [vecscale 0.9 [vecsub $end $start]]] ;# arrow made of a cylinder and a cone
	graphics $mol cylinder $start $middle radius 0.15 resolution 500
	graphics $mol cone $middle $end radius 0.25 resolution 500
}
################################################################################

set nf [molinfo top get numframes]

puts $outfile1 "intersection at -1 = either no intersection or parallel"
puts $outfile1 "intersection at -2 = no intersection"
puts $outfile1 "(Note that we are interested in the intersection between line segment (not vector) and triangle plane (infinite))"
puts $outfile1 "loop penetrating 0 = no penetration, loop penetrating 1 = yes"
puts $outfile1 "\n"

for {set i 0} {$i < $nf} {incr i} {
    puts $outfile1 "Time Frame $i ..."
    # CA atoms in the loop-forming sequence
    set resLoop [atomselect top "name CA and (resid 1 to 8)"]
	$resLoop frame $i
	$resLoop update
    set comLoop [measure center $resLoop weight mass]
    $resLoop delete
    for {set p 1} {$p <= 8} {incr p} {
        set res($p) [atomselect top "resid $p and name CA"]
        $res($p) frame $i
	    $res($p) update
        set com($p) [measure center $res($p) weight mass]
        $res($p) delete
    }
    # CA atoms in the tail sequence
    for {set q 14} {$q <= 21} {incr q} {
        set res($q) [atomselect top "resid $q and name CA"]
        $res($q) frame $i
	    $res($q) update
        set com($q) [measure center $res($q) weight mass]
        $res($q) delete
    }

    # Normal vectors of each triangle
    for {set r 1} {$r <= 8} {incr r} {
        if {$r < 8} {
            set rNext [expr $r + 1]
            set norm($r) [triNorm $comLoop $com($r) $com($rNext)]
            unset rNext
        } else {
            set rNext 1
            set norm($r) [triNorm $comLoop $com($r) $com($rNext)]
            unset rNext
        }
    }

    # Draw normal vectors (comment out if not needed)
    #for {set w 7} {$w <= 7} {incr w} {
    #    set molID [molinfo top get id]
    #    set scaledNorm [vecscale 0.4 $norm($w)]
    #    set end [vecadd $comLoop $scaledNorm]
    #    vmd_draw_arrow $molID $comLoop $end
    #    unset molID
    #    unset scaledNorm
    #    unset end
    #}

    # Tail vectors
    for {set s 14} {$s < 21} {incr s} {
        set sNext [expr $s + 1]
        set tail($s) [vecsub $com($sNext) $com($s)]
        unset sNext
    }

    # Find the intersecting point between each tail vector and each triangle plane
    set counter 0
    for {set u 1} {$u <= 8} {incr u} {
        for {set v 14} {$v < 21} {incr v} {
            set intersectionPoint [findIntersection $com($v) $comLoop $tail($v) $norm($u)]
            if {$intersectionPoint == -2 || $intersectionPoint == -1} {
                puts $outfile1 "\t intersection between plane $u and tail $v at $intersectionPoint ... loop penetrating N/A"
            } else {
                if {$u == 8} {
                    set uNext 1
                    set insideTriangle [barycentricTransform $comLoop $com($u) $com($uNext) $intersectionPoint]
                    puts $outfile1 "\t intersection between plane $u and tail $v at $intersectionPoint ... loop penetrating $insideTriangle"
                    if {$insideTriangle > 0} {
                        set counter [expr $counter + 1]
                    }
                    unset uNext
                    unset insideTriangle
                } else {
                    set uNext [expr $u + 1]
                    set insideTriangle [barycentricTransform $comLoop $com($u) $com($uNext) $intersectionPoint]
                    puts $outfile1 "\t intersection between plane $u and tail $v at $intersectionPoint ... loop penetrating $insideTriangle"
                    if {$insideTriangle > 0} {
                        set counter [expr $counter + 1]
                    }
                    unset uNext
                    unset insideTriangle
                }
            }
            unset intersectionPoint
        }
        puts -nonewline  $outfile1 "\n"
    }

    if {$counter > 0} {
        puts $outfile2 "Time $i ... Penetrating (Pre-Lasso)"
    } else {
        puts $outfile2 "Time $i ... Non-penetrating (Pre-Tadpole)"
    }

    unset comLoop
	unset com
    unset norm
    unset tail
    unset counter
}

close $outfile1
close $outfile2