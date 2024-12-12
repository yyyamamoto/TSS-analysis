#! /usr/bin/perl

#  2017.01.31
#  Yoshiharu Y. Yamamoto, Faculty of Applied Biological Sciences, Gifu University
#  yamamoto.yoshiharu.h5@f.gifu-u.ac.jp
#  Solexa TSS analysis
#  Preparation of TSS clusers from mapped data
#  rule 1: cut if distance is eq or more than 20 bp (=$gapdistance1)
#  rule 2: cut if distance between secondary peaks is eq or more than 100 bp (=$gapdistance2)
#  cutting site: bottom of P1 peaks, longest space among TSS in the bottom
#  ref: Tokizawa M, Kusunoki K, Koyama H, Kurotani A, Sakurai T, Suzuki Y, Kurata T, Yamamoto YY. "Identification of Arabidopsis genic and non-genic promoters by pair-end sequencing of TSS tags". Plant J, 2017.

#  argument 1: input file name
#  input file: output file from TSSClusteringFromMappedData.pl
#  1st column: chromosome number
#  2nd column: chromosomal position
#  3rd column: direction (F or R)
#  4th column: number of TSS tag
#  5th column: cluster ID (temporal, rename recommended)
#  6th column: serial number of TSS
#  7th column: u/s/e/d mark (ignore this)
#  8th column: P1 peak mark
#  9th column: P1 serial number
#  10th column: P2 peak mark
#  11th column: P2 serial number
#  12th column: CP peak mark (=cluster peak)
#  12th column: P3 peak mark                  <= 180122YYY
#  13th column: P3 serial number              <= 180122YYY
#  14th column: P4 peak mark                  <= 180122YYY
#  15th column: P4 serial number              <= 180122YYY
#  16th column: CP peak mark (=cluster peak)  <= 180122YYY

#  argument 2: output file name
#  output file:
#  1st column: Cluster ID (serial)
#  2nd column: chromosome number
#  3rd column: direction (F or R)
#  4th column: North end of the cluster
#  5th column: South end of the cluster
#  6th column: cluster width
#  7th column: position of peak TSS
#  8th column: TSS tags of the cluster (total)
#  9th column: TSS tags of the peak TSS
# 10th column: peak ratio



# open and read table file
$tssfile1 = $ARGV[0];
open (TABLE1, "<$tssfile1") or die "cannnot find table file";
@table1 = <TABLE1>;
close (TABLE1);
chomp @table1;
$p0tablesize = @table1;

# output file
$outputfile = $ARGV[1];
open (OUT, ">$outputfile");

#$cl_name = "CL0A_";

#make 2D table
for ($i = 0; $i < $p0tablesize ; $i++){
    @{$p0_table[$i]} = split (/\t/, $table1[$i]);
}


#summerize for each cluster
$cluster_num = 1;
$p1 = 0; # start number/position of the focused cluster
$p2 = 1; # number of TSS in the cluster
$out_table = ""; $totaltag =0;

while ($p1 < $p0tablesize){
    while ($p0_table[$p1][4] eq $p0_table[$p1 + $p2][4]){ #determine p2
        $p2++;
    }
    for ($i = 0;$i < $p2; $i++){
        $totaltag = $totaltag + $p0_table[$p1 + $i][3];
        if ($p0_table[$p1 + $i][15] eq CP){
            $peaktag = $p0_table[$p1 + $i][3];
            $peakposi = $p0_table[$p1 + $i][1];
        }
    }
    #$clusterID = sprintf("%06d",$cluster_num);
    $cluster_width = $p0_table[$p1 + $p2 - 1][1] - $p0_table[$p1][1] +1;
    $peakratio = $peaktag / $totaltag;
    $out_table = $out_table."$p0_table[$p1][4]"."\t"."$p0_table[$p1][0]"."\t"."$p0_table[$p1][2]"."\t"."$p0_table[$p1][1]"."\t"."$p0_table[$p1+$p2 - 1][1]"."\t"."$cluster_width"."\t"."$peakposi"."\t"."$totaltag"."\t"."$peaktag"."\t"."$peakratio"."\n";
    
    for ($i1 = $p1; $i1 < p1 + $p2 - 1; $i1++){ #try to release memory
        for ($i2 = 0; $i2 < 16; $i2++){
            $p0_table[$i1][$i2] = "";
        }
    }
    
    $cluster_num++; $p1 = $p1 + $p2;$p2 = 1;$totaltag = 0; #reset for the next round
    print "$cluster_num\n";
}

print OUT $out_table;

close (OUT);

#output file = @tablesplit
