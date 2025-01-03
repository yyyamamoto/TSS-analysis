#! /usr/bin/perl

#  2017.01.31
#  2018.01.22 modified for detection of promoter switching
#  2024.10.   corrected CP marking of chopped (new) clusters  (problem = no CP marks for them)
#  Yoshiharu Y. Yamamoto, Faculty of Applied Biological Sciences, Gifu University
#  yamamoto.yoshiharu.h5@f.gifu-u.ac.jp
#  Solexa TSS analysis
#  Preparation of TSS clusers from mapped data
#  rule 1: cut if distance is eq or more than 10 bp (=$gapdistance1)
#  rule 2: cut if distance between secondary/tertiary/quaternary peaks is eq or more than 100 bp (=$gapdistance2)
#  cutting site: bottom of P1 peaks, longest space among TSS in the bottom
#  ref: Tokizawa M, Kusunoki K, Koyama H, Kurotani A, Sakurai T, Suzuki Y, Kurata T, Yamamoto YY. "Identification of Arabidopsis genic and non-genic promoters by pair-end sequencing of TSS tags". Plant J, 2017.

#  argument 1: input file name
#  input file: mapped TSS data (tab delimited table)
#  1st column: chromosome number (Chr01~)
#  2nd column: chromosomal position
#  3rd column: direction (F or R)
#  4th column: number of TSS tag
## Table should be sorted according to 1)chromosome number, 2)direction, and 3)chromosomal position in this order.

#  argument 2: output file name
#  output file:
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
#  12th column: P3 peak mark                  <= 180122YYY
#  13th column: P3 serial number              <= 180122YYY
#  14th column: P4 peak mark                  <= 180122YYY
#  15th column: P4 serial number              <= 180122YYY
#  16th column: CP peak mark (=cluster peak)  <= 180122YYY


#conditions for cutting ### change them if needed
#$gapdistance1 = 20;                          #original setting
$gapdistance1 = 10;                           #<= 180122YYY modified for detection of promoter switch
$gapdistance2 = 100;

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

$cl_name_fin = "CL_";

#make 2D table
for ($i = 0; $i < $p0tablesize ; $i++){
    @{$p0_table[$i]} = split (/\t/, $table1[$i]);
}


#First chop, label $p0_table[][10] & mark $p0_table[][6] with u/d/e, mark s after marking P1
# u: up, d: down, e: even, s: start of a cluster
#write p0 serial# at [][5]
$ClusterNumber = 1;
for ($i = 0; $i < $p0tablesize - 1; $i++){
    $cl_name = "CL_$p0_table[$i][0]_";
    $p0_table[$i][4] = "$cl_name"."$ClusterNumber";
    $p0_table[$i][5] = $i;
    if ($p0_table[$i+1][0] eq $p0_table[$i][0]){ #chr
        if ($p0_table[$i+1][2] eq $p0_table[$i][2]){ #direct
            if($p0_table[$i][3] > $p0_table[$i+1][3]){
                $p0_table[$i+1][6] = "d";
            }
            elsif($p0_table[$i][3] < $p0_table[$i+1][3]){
                $p0_table[$i+1][6] = "u";
            }
            elsif($p0_table[$i][3] == $p0_table[$i+1][3]){
                $p0_table[$i+1][6] = "e";
            }
            
            if($p0_table[$i+1][1] - $p0_table[$i][1] > $gapdistance1){ #chop cluster, evaluated distance
                $ClusterNumber++;
            }
            
        }
        else {$ClusterNumber++;}
    }
    else {$ClusterNumber++;}
}
$p0_table[0][6] = "s";$p0_table[$p0tablesize-1][5] = $p0tablesize-1;# exception
if ($p0_table[$p0tablesize - 1][1] - $p0_table[$p0tablesize - 2][1] > $gapdistance1){# exception
    $cl_name ="CL_$p0_table[$p0tablesize - 1][0]_";
    $p0_table[$p0tablesize - 1][4] = "$cl_name"."$ClusterNumber";$p0_table[$p0tablesize - 1][6] = "s";#ClusterNumber: same as one-north
}
else {
    $cl_name ="CL_$p0_table[$p0tablesize - 1][0]_";
    $p0_table[$p0tablesize - 1][4] = "$cl_name"."$ClusterNumber";# $ClusterNumber: add one from one-north (L77)
}

print "first chop & u/d/e labeling done\n";

#mark P1 at $p0_table[][7], ud >> P1
for ($k = 0; $k < $p0tablesize - 1; $k++){
    if ($p0_table[$k][6] eq "u" && $p0_table[$k + 1][6] eq "d"){
        $p0_table[$k][7] ="P1";
    }
}


#mark $p0_table[][6] with s
for ($j = 0; $j < $p0tablesize -1; $j++){
    if ($p0_table[$j][4] ne $p0_table[$j + 1][4]){
        $p0_table[$j + 1][6] = "s";
    }
}

#mark P1 at sd & us
for ($i1 = 0; $i1 < $p0tablesize; $i1++){
    if ($p0_table[$i1][6] eq "s"){# 1 TSS in a cluster
        if ($p0_table[$i1][4] ne $p0_table[$i1+1][4]){# single TSS in a cluster, SS >> P1,P2,CP
            $p0_table[$i1][7] = "P1"; $p0_table[$i1][9] = "P2"; #$p0_table[$i1][15] = "CP";
        }
        elsif ($p0_table[$i1+1][6] eq "d"){# sd >> P1
            $p0_table[$i1][7] = "P1";
        }
        if ($p0_table[$i1-1][6] eq "u"){# us  >> P1
            $p0_table[$i1-1][7] = "P1";
        }
    }
}

#mark P1 at e
$j1 = 0;
while ($j1 < $p0tablesize){
    if ($p0_table[$j1][6] ne "e"){
        $j1++;
    }
    elsif ($p0_table[$j1][6] eq "e"){
        $j2 = 0;
        while ($p0_table[$j1+$j2][6] eq "e"){
            $j2++;
        }# $j2 = length of "e"; e = 1 >> $j2 = 1
        if ($p0_table[$j1 - 1][6] eq "d"){#deee
            #nothing to do
        }
        elsif ($p0_table[$j1 + $j2 - 1][6] eq "u"){# eeeu
            #nothing to do
        }
        elsif ($p0_table[$j1 - 1][6] eq "s" or $p0_table[$j1 - 1][6] eq "u"){# s/u_eee
            if ($p0_table[$j1 + $j2][6] eq "d" or $p0_table[$j1 + $j2][6] eq "s"){# s/u_eee  s/d
                $mid = 0;
                $mid = $j1 - 1 + int($j2 / 2);# cut in the middle if $j2 is odd; cut at two-north of the middle if even
                $p0_table[$mid][7] = "P1";
            }
        }
        $j1 = $j1 + $j2 + 1;
    }
}

##mark P2


#prep P1 table
%p1_table = "";    $p1n = 0;
for ($k = 0; $k < $p0tablesize; $k++){
    if ($p0_table[$k][7] eq "P1"){
        $p1_table[$p1n][0] = $p0_table[$k][0];$p1_table[$p1n][1] = $p0_table[$k][1];$p1_table[$p1n][2] = $p0_table[$k][2];$p1_table[$p1n][3] = $p0_table[$k][3];$p1_table[$p1n][4] = $p0_table[$k][4];$p1_table[$p1n][5] = $p0_table[$k][5];$p1_table[$p1n][6] = "";$p1_table[$p1n][7] = $p0_table[$k][7];$p1_table[$p1n][8] = $p1n; $p1_table[$p1n][9] = $p0_table[$k][9];$p1_table[$p1n][11] = $p0_table[$k][11];$p1_table[$p1n][10] = "";
        $p0_table[$k][8] = $p1n;
        $p1n++;
    }
}
$p1tablesize = $p1n;

#mark P2

for ($i = 0; $i < $p1tablesize; $i++){# clear $p1_table[][6]
    $p1_table[$i][6] = "";
}
#mark u/d/e at $p1_table[][6]
for ($i = 0;$i < $p1tablesize - 1; $i++){
    if ($p1_table[$i][3] < $p1_table[$i + 1][3]){
        $p1_table[$i + 1][6] = "u";
    }
    elsif ($p1_table[$i][3] > $p1_table[$i + 1][3]){
        $p1_table[$i + 1][6] = "d";
    }
    elsif ($p1_table[$i][3] == $p1_table[$i + 1][3]){
        $p1_table[$i + 1][6] = "e";
    }
}
#mark P2 1st time
for ($i = 0; $i < $p1tablesize -1; $i++){#single peak, ud >> P2
    if ($p1_table[$i][6] eq "u" && $p1_table[$i + 1][6] eq "d"){
        $p0serialtemp = $p1_table[$i][5];
        $p0_table[$p0serialtemp][9] = "P2";
    }
}
#mark s at $p1_table[][6]
for ($i = 0; $i < $p1tablesize; $i++){
    if ($p1_table[$i][4] ne $p1_table[$i + 1][4]){
        $p1_table[$i + 1][6] = "s";
    }
}
$p1_table[0][6] = "s"; #treat exceptions

#mark P2 at sd & us
for ($i1 = 0; $i1 < $p1tablesize; $i1++){
    if ($p1_table[$i1][6] eq "s"){# 1 TSS in a cluster
        if ($p1_table[$i1][4] ne $p1_table[$i1 + 1][4]){# single TSS in a cluster, SS >> P1,P2,CP
            $p0serialtemp = $p1_table[$i1][5];
            $p0_table[$p0serialtemp][9] = "P2"; #$p0_table[$p0serialtemp][15] = "CP";
        }
        elsif ($p1_table[$i1 + 1][6] eq "d"){# sd >> P1
            $p0serialtemp = $p1_table[$i1][5];
            $p0_table[$p0serialtemp][9] = "P2";
        }
        if ($i1 - 1 >= 0 && $p1_table[$i1 - 1][6] eq "u"){# us  >> P1
            $p0serialtemp = $p1_table[$i1 - 1][5];
            if ($p0serialtemp >= 0){
                $p0_table[$p0serialtemp][9] = "P2";
            }
            
        }
    }
}

#mark P2 at e
$j1 = 0;
while ($j1 < $p1tablesize){
    if ($p1_table[$j1][6] ne "e"){
        $j1++;
    }
    elsif ($p1_table[$j1][6] eq "e"){
        $j2 = 0;
        while ($p1_table[$j1 + $j2][6] eq "e"){
            $j2++;
        }# $j2 = length of "e"
        if ($p1_table[$j1 - 1][6] eq "d"){#deee
            #nothing done
        }
        elsif ($p1_table[$j1 + $j2][6] eq "u"){# eeeu
            #nothing done
        }
        elsif ($p1_table[$j1 - 1][6] eq "s" or $p1_table[$j1 - 1][6] eq "u"){# s/u_eee
            if ($p1_table[$j1 + $j2][6] eq "d" or $p1_table[$j1 + $j2][6] eq "s"){# s/u_eee  s/d
                $mid = 0;
                $mid = $j1 - 1 + int ($j2 / 2);# middle if $j2 is odd, two-north of the middle if even
                $p0serialtemp = $p1_table[$mid][5];
                $p0_table[$p0serialtemp][9] = "P2";
            }
        }
        $j1 = $j1 + $j2 + 1;
    }
}

#mark P2 for single P1 in a cluster
for ($i = 0; $i < $p1tablesize - 2; $i++){
    if ($p1_table[$i][4] ne $p1_table[$i + 1][4] && $p1_table[$i + 1][4] ne $p1_table[$i + 2][4]){
        $p0serialtemp = $p1_table[$i + 1][5];
        $p0_table[$p0serialtemp][9] = "P2";
        #$p0_table[$p0serialtemp][15] = "CP"; #<=180122YYY
    }
}

print "P2 marking done\n";

##Second chop at P2 valley

#extract P2 peaks
$p2serialtemp = 0;
for ($i = 0; $i < $p0tablesize; $i++){
    if ($p0_table[$i][9] eq "P2"){
        $p2_table[$p2serialtemp][0] = $p0_table[$i][0]; $p2_table[$p2serialtemp][1] = $p0_table[$i][1]; $p2_table[$p2serialtemp][2] = $p0_table[$i][2]; $p2_table[$p2serialtemp][3] = $p0_table[$i][3]; $p2_table[$p2serialtemp][4] = $p0_table[$i][4]; $p2_table[$p2serialtemp][5] = $p0_table[$i][5]; $p2_table[$p2serialtemp][6] = ""; $p2_table[$p2serialtemp][7] = $p0_table[$i][7]; $p2_table[$p2serialtemp][8] = $p0_table[$i][8]; $p2_table[$p2serialtemp][9] = $p0_table[$i][9]; $p2_table[$p2serialtemp][15] = $p0_table[$i][15];
        $p0_table[$i][10] = $p2serialtemp;
        $p2serialtemp++;
    }
}
$p2tablesize = $p2serialtemp;


#chop at P1 valley, rename clusterID: G#N, G#SN, G#SS
#$p2, $p1n, $p1s, $p1v, p0last_d, $p0first_u,$p0cn, $p0cs
for ($p2 = 0; $p2 < $p2tablesize - 1; $p2++){
    if ($p2_table[$p2][4] eq $p2_table[$p2 + 1][4]){ #same cluster
        $p1n = $p2_table[$p2][8]; $p1s = $p2_table[$p2 + 1][8];
        $p1first_u = 0; $p1last_d = 0; $p1cn = 0; $p1cs = 0;#initialize,
        if ($p2_table[$p2 + 1][1] - $p2_table[$p2][1] > $gapdistance2){ # distance between P2 peaks
            print "cut at P2_$p2_table[$p2][4]\n";
            #remake, determine valley after examination of $p1_table[][12]のs/u/d/e.
            # set $p0last_d & $p0first_u
            for ($p1v = $p1n; $p1v < $p1s + 1; $p1v++){
                if ($p1_table[$p1v][6] eq "d"){
                    $p1last_d = $p1_table[$p1v][8];
                }
            }
            for ($p1v = $p1s; $p1v > $p1n; $p1v--){
                if ($p1_table[$p1v][6] eq "u"){
                    $p1first_u = $p1_table[$p1v][8];
                }
            } # selected P1 valley = $p1last_d - 1 ~ $p1first_u
            
            #select largest gap in the P1 valley
            @distance = ""; %dis2p1serial = "";
            for ($p1vv = $p1last_d - 1; $p1vv < $p1first_u; $p1vv++){
                $dis = $p1_table[$p1vv + 1][1] - $p1_table[$p1vv][1];
                push (@distance, $dis);
                $dis2p1serial{$dis} = $p1vv;
            }
            @distance = sort {$b <=> $a} @distance;
            $max_dis = $distance[0];
            $p1cutn = $dis2p1serial{$max_dis}; $p1cuts = $p1cutn + 1;# cutting site in P1 valley
            $p0n = $p1_table[$p1cutn][5]; $p0s = $p1_table[$p1cuts][5]; #
            
            #search largest gap between $p0cn and $p0cs
            @distance = ""; %dis2p0serial = "";
            for ($p0vv = $p0n; $p0vv < $p0s; $p0vv++){
                $dis = $p0_table[$p0vv + 1][1] - $p0_table[$p0vv][1];
                push (@distance, $dis);
                $dis2p0serial{$dis} = $p0vv;
            }
            @distance = sort {$b <=> $a} @distance;
            $max_dis = $distance[0];
            $p0c = $dis2p0serial{$max_dis};
            $clusterID = $p0_table[$p0c][4];
            $p0_c1 = $p0c;
            
            while ($p0_table[$p0_c1][4] eq $clusterID){#rename cluterID, North of the border
                $p0_table[$p0_c1][4] = "$clusterID"."N";
                $p0_c1--;
            }
            
            $p0_c2 = $p0c + 1;
            while ($p0_table[$p0_c2][4] eq $clusterID){#rename cluterID, South of the border
                $p0_table[$p0_c2][4] = "$clusterID"."S";
                $p0_c2++;
            }
        }
    }
}

print "second chop done\n";


# renew p1_table[][4] (=clusterID) after chopping
for ($i = 0; $i < $p0tablesize; $i++){
    if ($p0_table[$i][7] eq "P1"){
        $p1_table[$p0_table[$i][8]][4] = $p0_table[$i][4];
    }
}
# renew p2_table[][4] (=clusterID) after chopping
for ($i = 0; $i < $p0tablesize; $i++){
    if ($p0_table[$i][9] eq "P2"){
        $p2_table[$p0_table[$i][10]][4] = $p0_table[$i][4];
    }
}
# renew p3_table[][4] (=clusterID) after chopping
for ($i = 0; $i < $p0tablesize; $i++){
    if ($p0_table[$i][11] eq "P3"){
        $p3_table[$p0_table[$i][12]][4] = $p0_table[$i][4];
    }
}

#mark CP
for ($i = 0; $i < $p2tablesize; $i++){#clear p2_table[][]
    @{$p2_table[$i]} = "";
}
$p2serialtemp = 0;
for ($i = 0; $i < $p0tablesize; $i++){#remake p2_table[][] after modification of clusterID.
    if ($p0_table[$i][9] eq "P2"){
        $p2_table[$p2serialtemp][0] = $p0_table[$i][0]; $p2_table[$p2serialtemp][1] = $p0_table[$i][1]; $p2_table[$p2serialtemp][2] = $p0_table[$i][2]; $p2_table[$p2serialtemp][3] = $p0_table[$i][3]; $p2_table[$p2serialtemp][4] = $p0_table[$i][4]; $p2_table[$p2serialtemp][5] = $p0_table[$i][5]; $p2_table[$p2serialtemp][6] = ""; $p2_table[$p2serialtemp][7] = $p0_table[$i][7]; $p2_table[$p2serialtemp][8] = $p0_table[$i][8]; $p2_table[$p2serialtemp][9] = $p0_table[$i][9]; $p2_table[$p2serialtemp][10] = $p0_table[$i][10]; $p2_table[$p2serialtemp][15] = $p0_table[$i][15];#<= 180122YYY
        $p2serialtemp++;
    }
}

### 180122YYY=>
### mark P3
## mark u/d/e at $p2_table[][6]
for ($i = 0; $i < $p2tablesize; $i++){
    $p2_table[$i][6] = "";
}

#mark u/d/e at $p2_table[][6]
for ($i = 0;$i < $p2tablesize - 1; $i++){
    if ($p2_table[$i][3] < $p2_table[$i + 1][3]){
        $p2_table[$i + 1][6] = "u";
    }
    elsif ($p2_table[$i][3] > $p2_table[$i + 1][3]){
        $p2_table[$i + 1][6] = "d";
    }
    elsif ($p2_table[$i][3] == $p2_table[$i + 1][3]){
        $p2_table[$i + 1][6] = "e";
    }
}

#mark s at $p2_table[][6]   overwrite
for ($i = 0; $i < $p2tablesize; $i++){
    if ($p2_table[$i][4] ne $p2_table[$i + 1][4]){
        $p2_table[$i + 1][6] = "s";
    }
}
$p2_table[0][6] = "s"; #treat exceptions

#mark P3 1st time  at ud
for ($i = 0; $i < $p2tablesize -1; $i++){#single peak, ud >> P2
    if ($p2_table[$i][6] eq "u" && $p2_table[$i + 1][6] eq "d"){
        $p0serialtemp = $p2_table[$i][5];
        $p0_table[$p0serialtemp][11] = "P3";
    }
}

#mark P3 at sd & us
for ($i1 = 0; $i1 < $p2tablesize; $i1++){
    if ($p2_table[$i1][6] eq "s"){# 1 TSS in a cluster
        if ($p2_table[$i1][4] ne $p2_table[$i1 + 1][4]){# single TSS in a cluster, SS >> P1,P2,CP
            $p0serialtemp = $p2_table[$i1][5];
            $p0_table[$p0serialtemp][11] = "P3";
        }
        elsif ($p2_table[$i1 + 1][6] eq "d"){# sd >> P1
            $p0serialtemp = $p2_table[$i1][5];
            $p0_table[$p0serialtemp][11] = "P3";
        }
        if ($i1 - 1 >= 0 && $p2_table[$i1 - 1][6] eq "u"){# us  >> P1
            $p0serialtemp = $p2_table[$i1 - 1][5];
            if ($p0serialtemp >= 0){
                $p0_table[$p0serialtemp][11] = "P3";
            }
        }
    }
}

#mark P3 at e
$j1 = 0;
while ($j1 < $p2tablesize){
    if ($p2_table[$j1][6] ne "e"){
        $j1++;
    }
    elsif ($p2_table[$j1][6] eq "e"){
        $j2 = 0;
        while ($p2_table[$j1 + $j2][6] eq "e"){
            $j2++;
        }# $j2 = length of "e"
        if ($p2_table[$j1 - 1][6] eq "d"){#deee
            #nothing done
        }
        elsif ($p2_table[$j1 + $j2][6] eq "u"){# eeeu
            #nothing done
        }
        elsif ($p2_table[$j1 - 1][6] eq "s" or $p2_table[$j1 - 1][6] eq "u"){# s/u_eee
            if ($p2_table[$j1 + $j2][6] eq "d" or $p2_table[$j1 + $j2][6] eq "s"){# s/u_eee  s/d
                $mid = 0;
                $mid = $j1 - 1 + int ($j2 / 2);# middle if $j2 is odd, two-north of the middle if even
                $p0serialtemp = $p2_table[$mid][5];
                $p0_table[$p0serialtemp][11] = "P3";
            }
        }
        $j1 = $j1 + $j2 + 1;
    }
}

#mark P3 for single P2 in a cluster
for ($i = 0; $i < $p2tablesize - 2; $i++){
    if ($p2_table[$i][4] ne $p2_table[$i + 1][4] && $p2_table[$i + 1][4] ne $p2_table[$i + 2][4]){
        $p0serialtemp = $p2_table[$i + 1][5];
        $p0_table[$p0serialtemp][11] = "P3";
    }
}

print "P3 marking done\n";


###  mark P4

# prep p3_table
$p3serialtemp = 0;
for ($i = 0; $i < $p0tablesize; $i++){
    if ($p0_table[$i][11] eq "P3"){
        $p3_table[$p3serialtemp][0] = $p0_table[$i][0]; $p3_table[$p3serialtemp][1] = $p0_table[$i][1]; $p3_table[$p3serialtemp][2] = $p0_table[$i][2]; $p3_table[$p3serialtemp][3] = $p0_table[$i][3]; $p3_table[$p3serialtemp][4] = $p0_table[$i][4]; $p3_table[$p3serialtemp][5] = $p0_table[$i][5]; $p3_table[$p3serialtemp][6] = ""; $p3_table[$p3serialtemp][7] = $p0_table[$i][7]; $p3_table[$p3serialtemp][8] = $p0_table[$i][8]; $p3_table[$p3serialtemp][9] = $p0_table[$i][9]; $p3_table[$p3serialtemp][10] = $p0_table[$i][10]; $p3_table[$p3serialtemp][11] = $p0_table[$i][11];  $p3_table[$p3serialtemp][12] = $p3serialtemp; $p0_table[$i][12] = $p3serialtemp; $p3_table[$p3serialtemp][15] = $p0_table[$i][15];#<= 180122YYY
        $p3serialtemp++;
    }
}

$p3tablesize = $p3serialtemp;

#mark u/d/e at $p3_table[][6]
for ($i = 0;$i < $p3tablesize - 1; $i++){
    if ($p3_table[$i][3] < $p3_table[$i + 1][3]){
        $p3_table[$i + 1][6] = "u";
    }
    elsif ($p3_table[$i][3] > $p3_table[$i + 1][3]){
        $p3_table[$i + 1][6] = "d";
    }
    elsif ($p3_table[$i][3] == $p3_table[$i + 1][3]){
        $p3_table[$i + 1][6] = "e";
    }
}

#mark s at $p3_table[][6]   overwrite
for ($i = 0; $i < $p3tablesize; $i++){
    if ($p3_table[$i][4] ne $p3_table[$i + 1][4]){
        $p3_table[$i + 1][6] = "s";
    }
}
$p3_table[0][6] = "s"; #treat exceptions


#mark P4 1st time
for ($i = 0; $i < $p2tablesize -1; $i++){#single peak, ud >> P2
    if ($p3_table[$i][6] eq "u" && $p3_table[$i + 1][6] eq "d"){
        $p0serialtemp = $p3_table[$i][5];
        $p0_table[$p0serialtemp][13] = "P4";
    }
}

#mark P4 at sd & us
for ($i1 = 0; $i1 < $p3tablesize; $i1++){
    if ($p3_table[$i1][6] eq "s"){# 1 TSS in a cluster
        if ($p3_table[$i1][4] ne $p3_table[$i1 + 1][4]){# single TSS in a cluster, SS >> P1,P2,CP
            $p0serialtemp = $p3_table[$i1][5];
            $p0_table[$p0serialtemp][13] = "P4"; #$p0_table[$p0serialtemp][15] = "CP";
        }
        elsif ($p3_table[$i1 + 1][6] eq "d"){# sd >> P1
            $p0serialtemp = $p3_table[$i1][5];
            $p0_table[$p0serialtemp][13] = "P4";
        }
        if ($i1 - 1 >= 0 && $p3_table[$i1 - 1][6] eq "u"){# us  >> P1
            $p0serialtemp = $p3_table[$i1 - 1][5];
            if ($p0serialtemp >= 0){
                $p0_table[$p0serialtemp][13] = "P4";
            }
        }
    }
}

#mark P4 at e
$j1 = 0;
while ($j1 < $p3tablesize){
    if ($p3_table[$j1][6] ne "e"){
        $j1++;
    }
    elsif ($p3_table[$j1][6] eq "e"){
        $j2 = 0;
        while ($p3_table[$j1 + $j2][6] eq "e"){
            $j2++;
        }# $j2 = length of "e"
        if ($p3_table[$j1 - 1][6] eq "d"){#deee
            #nothing done
        }
        elsif ($p3_table[$j1 + $j2][6] eq "u"){# eeeu
            #nothing done
        }
        elsif ($p3_table[$j1 - 1][6] eq "s" or $p3_table[$j1 - 1][6] eq "u"){# s/u_eee
            if ($p3_table[$j1 + $j2][6] eq "d" or $p3_table[$j1 + $j2][6] eq "s"){# s/u_eee  s/d
                $mid = 0;
                $mid = $j1 - 1 + int ($j2 / 2);# middle if $j2 is odd, two-north of the middle if even
                $p0serialtemp = $p3_table[$mid][5];
                $p0_table[$p0serialtemp][13] = "P4";
            }
        }
        $j1 = $j1 + $j2 + 1;
    }
}

#mark P4 for single P3 in a cluster
for ($i = 0; $i < $p2tablesize - 2; $i++){
    if ($p3_table[$i][4] ne $p3_table[$i + 1][4] && $p3_table[$i + 1][4] ne $p3_table[$i + 2][4]){
        $p0serialtemp = $p3_table[$i + 1][5];
        $p0_table[$p0serialtemp][13] = "P4";
        #$p0_table[$p0serialtemp][15] = "CP";
    }
}

$p4serialtemp = 0;
for ($i = 0; $i < $p0tablesize; $i++){
    if ($p0_table[$i][13] eq "P4"){
        $p0_table[$i][14] = $p4serialtemp;
        $p4serialtemp++;
    }
}
$p4tablesize = $p4serialtemp;


## chop at P3 ~ P3 gap, rename clusterID: G#N, G#SN, G#SS
#$p3, $p2n, $p2s, $p2v, p2last_d, $p2first_u,$p1cn, $p1cs
for ($p3 = 0; $p3 < $p3tablesize - 1; $p3++){
    if ($p3_table[$p3][4] eq $p3_table[$p3 + 1][4]){ #same cluster
        if ($p3_table[$p3 + 1][1] - $p3_table[$p3][1] > $gapdistance2){ # distance between P2 peaks
            print "cut at P3_$p3_table[$p3][4]\n";
            $p0n = ""; $p0s = ""; $p2n = ""; $p2s = ""; $p2first_u = ""; $p2last_d = ""; $p2cn = ""; $p2cs = ""; $p2n = $p3_table[$p3][10]; $p2s = $p3_table[$p3 + 1][10]; #initialize,
            #remake, determine valley after examination of $p2_table[][6]のs/u/d/e.
            # set $p2last_d & $p2first_u
            for ($p2v = $p2n + 1; $p2v < $p2s + 1; $p2v++){
                if ($p2_table[$p2v][6] eq "d"){
                    $p2last_d = $p2_table[$p2v][10];
                }
            }
            for ($p2v = $p2s; $p2v > $p2n; $p2v--){
                if ($p2_table[$p2v][6] eq "u"){
                    $p2first_u = $p2_table[$p2v][10];
                }
            }
            # selected P2 valley (P2~P2) = $p2last_d ~ $p2first_u - 1
            #search largest gap between $p0n(= p0serial# of $p2last_d) and $p0s (= p0serial# of $p2first_u - 1)
            $p0_c1 = ""; $p0_c2 = ""; $p0c = ""; $max_dis = "";$clusterID = "";
            @distance = ""; %dis2p0serial = "";
            $p0n = $p2_table[$p2last_d - 1][5]; $p0s = $p2_table[$p2first_u][5];
            
            for ($p0vv = $p0n; $p0vv < $p0s; $p0vv++){
                $dis = $p0_table[$p0vv + 1][1] - $p0_table[$p0vv][1];
                push (@distance, $dis);
                $dis2p0serial{$dis} = $p0vv;
            }
            @distance = sort {$b <=> $a} @distance;
            $max_dis = $distance[0];
            $p0c = $dis2p0serial{$max_dis};
            $clusterID = $p0_table[$p0c][4];
            $p0_c1 = $p0c;
            
            while ($p0_table[$p0_c1][4] eq $clusterID){#rename cluterID, North of the border
                $p0_table[$p0_c1][4] = "$clusterID"."N";
                $p0_c1--;
            }
            
            $p0_c2 = $p0c + 1;
            while ($p0_table[$p0_c2][4] eq $clusterID){#rename cluterID, South of the border
                $p0_table[$p0_c2][4] = "$clusterID"."S";
                $p0_c2++;
            }
        }
    }
}

print "third chop between P3 peaks done\n";


# prep p4_table
$p4serialtemp = 0;
for ($i = 0; $i < $p0tablesize; $i++){
    if ($p0_table[$i][13] eq "P4"){
        $p4_table[$p4serialtemp][0] = $p0_table[$i][0]; $p4_table[$p4serialtemp][1] = $p0_table[$i][1]; $p4_table[$p4serialtemp][2] = $p0_table[$i][2]; $p4_table[$p4serialtemp][3] = $p0_table[$i][3]; $p4_table[$p4serialtemp][4] = $p0_table[$i][4]; $p4_table[$p4serialtemp][5] = $p0_table[$i][5]; $p4_table[$p4serialtemp][6] = ""; $p4_table[$p4serialtemp][7] = $p0_table[$i][7]; $p4_table[$p4serialtemp][8] = $p0_table[$i][8]; $p4_table[$p4serialtemp][9] = $p0_table[$i][9]; $p4_table[$p4serialtemp][10] = $p0_table[$i][10]; $p4_table[$p4serialtemp][11] = $p0_table[$i][11];  $p4_table[$p4serialtemp][12] = $p0_table[$i][12]; $p4_table[$p4serialtemp][13] = $p0_table[$i][13]; $p4_table[$p4serialtemp][14] = $p4serialtemp; $p0_table[$i][14] = $p4serialtemp; $p4_table[$p4serialtemp][15] = $p0_table[$i][15];#<= 180122YYY
        $p4serialtemp++;
    }
}

# renew p3_table[][4] (=clusterID) after chopping
$p3serialtemp = 0;
for ($i = 0; $i < $p0tablesize; $i++){
    if ($p0_table[$i][11] eq "P3"){
        $p3_table[$p3serialtemp][4] = $p0_table[$i][4];
        $p3serialtemp++;
    }
}

## chop at P4 ~ P4 gap, rename clusterID: G#N, G#SN, G#SS
#$p4, $p3n, $p3s, $p3v, p3last_d, $p3first_u,$p0cn, $p0cs
## chop at P3 ~ P3 gap, rename clusterID: G#N, G#SN, G#SS
#$p3, $p2n, $p2s, $p2v, p2last_d, $p2first_u,$p1cn, $p1cs

for ($p4 = 0; $p4 < $p4tablesize - 1; $p4++){
    if ($p4_table[$p4][4] eq $p4_table[$p4 + 1][4]){ #same cluster

        if ($p4_table[$p4 + 1][1] - $p4_table[$p4][1] > $gapdistance2){ # distance between P2 peaks
            print "cut at P4_$p4_table[$p4][4]\n";
            $p3n = ""; $p3s = ""; $p3first_u = ""; $p3last_d = ""; $p3cn = ""; $p3cs = ""; $p3n = $p4_table[$p4][12]; $p3s = $p4_table[$p4 + 1][12]; #initialize,
            #remake, determine valley after examination of $p3_table[][6]のs/u/d/e.
            # set $p2last_d & $p2first_u
            for ($p3v = $p3n + 1; $p3v < $p3s + 1; $p3v++){
                if ($p3_table[$p3v][6] eq "d"){
                    $p3last_d = $p3_table[$p3v][12];
                }
            }
            for ($p3v = $p3s; $p3v > $p3n; $p3v--){
                if ($p3_table[$p3v][6] eq "u"){
                    $p3first_u = $p3_table[$p3v][12];
                }
            }
            # selected P3 valley (P3~P3) = $32last_d ~ $p3first_u - 1
            #search largest gap between $p0n(= p0serial# of $p3last_d) and $p0s (= p0serial# of $p3first_u - 1)
            
            
            $p0n = ""; $p0s = ""; $p0_c1 = ""; $p0_c2 = ""; $p0c = ""; $max_dis = "";$clusterID = "";
            @distance = ""; %dis2p0serial = "";
            $p0n = $p3_table[$p3last_d - 1][5]; $p0s = $p3_table[$p3first_u][5];
            
            for ($p0vv = $p0n; $p0vv < $p0s; $p0vv++){
                $dis = $p0_table[$p0vv + 1][1] - $p0_table[$p0vv][1];
                push (@distance, $dis);
                $dis2p0serial{$dis} = $p0vv;
            }
            @distance = sort {$b <=> $a} @distance;
            $max_dis = $distance[0];
            $p0c = $dis2p0serial{$max_dis};
            $clusterID = $p0_table[$p0c][4];
            $p0_c1 = $p0c;
            
            while ($p0_table[$p0_c1][4] eq $clusterID){#rename cluterID, North of the border
                $p0_table[$p0_c1][4] = "$clusterID"."N";
                $p0_c1--;
            }
            
            $p0_c2 = $p0c + 1;
            while ($p0_table[$p0_c2][4] eq $clusterID){#rename cluterID, South of the border
                $p0_table[$p0_c2][4] = "$clusterID"."S";
                $p0_c2++;
            }
        }
    }
}

print "forth chop between P4 peaks done\n";


#rename ClusterID, mark CP

$new_clusternumber = 1;
$i1 = 0; $i2 = 1;
while ($i1 + $i2 <= $p0tablesize){
    
    while ($p0_table[$i1 + $i2][4] eq $p0_table[$i1][4]){
        $i2++;
    }
    # $i2 = number of TSSs in the cluster
    # $i1 = start position of the cluster in p0
    
    # rename ClusterID, mark CP
    # select CP from P4
    $cluterID_new = sprintf("%07d",$new_clusternumber);
    @cp_search = ""; %tag2p0serial = ""; $maxtag = "";
    for ($j = 0; $j < $i2; $j++){
        $p0_table[$i1 + $j][4] = "$cl_name_fin"."$cluterID_new"; #rename ClusterID
        if ($p0_table[$i1 + $j][13] eq "P4"){
            push (@cp_search, $p0_table[$i1 + $j][3]);
            $tag2p0serial{$p0_table[$i1 + $j][3]} = $i1 + $j;
        }
    }
    @cp_search = sort {$b <=> $a} @cp_search;
    $maxtag = $cp_search[0];
    $p0_table[$tag2p0serial{$maxtag}][15] = "CP";

    $i1 = $i1 + $i2; $i2 = 1; #next cluster
    $new_clusternumber++;
}

#############################################################################

##Corrected in 2024 Oct
##CP addition to new clsters (=last chopped clusters)
# find cluster with no CP in p0_table[][15]
# prep $clustertable[][], size [0~40k][0,1]
# $clustertable[][0]=clusterID,$clustertable[][1]=CP or none
%clustertable ="";
$ct_i = 0;# number of cluster in clustertable
$cID_prev = "some";
$cID_cur = "";
for ($i_p0 = 0; $i_p0 < $p0tablesize; $i_p0++ ){
    $cID_cur = $p0_table[$i_p0][4];
    if ($cID_cur ne $cID_prev){
        $clustertable[$ct_i][0] = $cID_cur;
        if ($p0_table[$i_p0][15] eq "CP"){
            $clustertable[$ct_i][1] = "CP";
        }
        $ct_i++;
    }
    else {
        if ($p0_table[$i_p0][15] eq "CP"){
            $clustertable[$ct_i - 1][1] = "CP";
        }
    }
    $cID_prev = $cID_cur;
    $cID_cur = "";
}
$ ct_i++;  #$ci_i = number of clusters in clustertable

#clustertable[][] prepared, size[$ct_i][0: clusterID,1: CP or none]

#extract CP-negative clusters, and mark CP1-3
for ($ii = 0; $ii < $ct_i; $ii++){
    $curr_cluster ="";
    if ($clustertable[$ii][1] ne "CP"){
        $curr_cluster = $clustertable[$ii][0];
        #no CP in this cluster
        #check p0_table[][], find peak, label p0_table[][15]=CP
        #copy TSS data of the cluster from p0_table[][]
        %selectedcluster = "";$v = 0; $selectedclustersize = "";
        for ($iii = 0; $iii < $p0tablesize; $iii++){
            if ($p0_table[$iii][4] eq $curr_cluster){
                for ($iv = 0; $iv < 16; $iv++){
                    $selectedcluster[$v][$iv] = $p0_table[$iii][$iv];#ここまだ
                }
                $v++;
            }
        }
        $selectedclustersize = $v; # problem when single TSS
        
        #when single TSS
        if ($selectedclustersize == 1){
            $posi_p0table = "";
            $posi_p0table = $selectedcluster[0][5];
            $p0_table[$posi_p0table][15] = "CP";
        }

        #when multi TSS
        else {
            $checkP1 = ""; @p1list =""; %height2serialno = "";@sortedp1 = "";$p_posi =""; $p1listsize = "";@p0list = "";$checkP3 = ""; $posi_P3 = "";
            for ($vi = 0; $vi < $selectedclustersize; $vi++){
                if ($selectedcluster[$vi][7] eq "P1"){
                    $ckeckP1 = "P1";
                    if ($selectedcluster[$vi][11] eq "P3"){
                        $checkP3 = "P3";
                        $posi_P3 = $selectedcluster[$vi][5];
                    }
                    push @p1list, $selectedcluster[$vi][3];
                    $height2serialno{$selectedcluster[$vi][3]} = $selectedcluster[$vi][5];
                    
                }
            }
            
            #when P1 exists,   select peak from P1 as CP2
            #if P4 or P3 exists, mark it as CP2
            # prefer P3 for CP2  241011
            if ($ckeckP1 eq "P1"){
                if ($checkP3 eq "P3"){
                    $p0_table[$posi_P3][15] = "CP";
                }
                else {
                    @sortedp1= sort {$b <=> $a} @p1list;
                    $p_posi = $height2serialno{$sortedp1[0]};
                    $p0_table[$p_posi][15] = "CP";
                }
            }
    
            #when P1 does not exist, pickup the highest tags in the cluster for CP3
            else {
                %height2serialno = ""; @p0list = ""; $p_posi = "";# reuse vars
                @sortedp0 ="";
                for ($viii = 0; $viii < $selectedclustersize; $viii++){
                    push @p0list, $selectedcluster[$viii][3];
                    $height2serialno{$selectedcluster[$viii][3]} = $selectedcluster[$viii][5];
                }
                @sortedp0 = sort {$b <=> $a} @p0list;
                $p_posi = $height2serialno{$sortedp0[0]};
                $p0_table[$p_posi][15] = "CP";
            }
        }
    }
}

#   (corrected in 2024 Oct//)
####################################################


#output file = @tablesplit
for ($i = 0; $i < $p0tablesize; $i++){
    for ($j = 0; $j < 16; $j++){ #<= 180122YYY
        print OUT "$p0_table[$i][$j]\t";
    }
    print OUT "\n";
}
##print p1_table
#for ($i = 0; $i < $p1tablesize; $i++){
#   for ($j = 0; $j < 16; $j++){ #<= 180122YYY
#       print OUT "$p1_table[$i][$j]\t";
#   }
#   print OUT "\n";
#}
##print p2_table
#for ($i = 0; $i < $p2tablesize; $i++){
#    for ($j = 0; $j < 16; $j++){ #<= 180122YYY
#        print OUT "$p2_table[$i][$j]\t";
#    }
#    print OUT "\n";
#}
##print p3_table
#for ($i = 0; $i < $p3tablesize; $i++){
#    for ($j = 0; $j < 16; $j++){ #<= 180122YYY
#        print OUT "$p3_table[$i][$j]\t";
#    }
#   print OUT "\n";
#}
##print p4_table
#for ($i = 0; $i < $p4tablesize; $i++){
#    for ($j = 0; $j < 16; $j++){ #<= 180122YYY
#        print OUT "$p4_table[$i][$j]\t";
#    }
#    print OUT "\n";
#}

close (OUT);



#output file = @tablesplit
