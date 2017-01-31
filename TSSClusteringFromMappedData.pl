#! /usr/bin/perl

#  2017.01.31
#  Yoshiharu Y. Yamamoto, Faculty of Applied Biological Sciences, Gifu University
#  yyy@gifu-u.ac.jp
#  Solexa TSS analysis
#  Preparation of TSS clusers from mapped data
#  rule 1: cut if distance is eq or more than 20 bp (=$gapdistance1)
#  rule 2: cut if distance between secondary peaks is eq or more than 100 bp (=$gapdistance2)
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
#  12th column: CP peak mark (=cluster peak)


#conditions for cutting ### change them if needed
$gapdistance1 = 20;
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

$cl_name = "CL0A_";

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

#mark P1 at $p0_table[][13], ud >> P1
for ($k = 0; $k < $p0tablesize - 1; $k++){
    if ($p0_table[$k][6] eq "u" && $p0_table[$k + 1][6] eq "d"){
        $p0_table[$k][7] ="P1";
    }
}


#mark $p0_table[][12] with s
for ($j = 0; $j < $p0tablesize -1; $j++){
    if ($p0_table[$j][4] ne $p0_table[$j + 1][4]){
        $p0_table[$j + 1][6] = "s";
    }
}

#mark P1 at sd & us
for ($i1 = 0; $i1 < $p0tablesize; $i1++){
    if ($p0_table[$i1][6] eq "s"){# 1 TSS in a cluster
        if ($p0_table[$i1][4] ne $p0_table[$i1+1][4]){# single TSS in a cluster, SS >> P1,P2,CP
            $p0_table[$i1][7] = "P1"; $p0_table[$i1][9] = "P2"; $p0_table[$i1][11] = "CP";
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
#mark P2 一回目
for ($i = 0; $i < $p1tablesize -1; $i++){#single peak, ud >> P2
    if ($p1_table[$i][6] eq "u" && $p1_table[$i + 1][6] eq "d"){
        $p0serialtemp = $p1_table[$i][5];
        $p0_table[$p0serialtemp][9] = "P2";
    }
}
#mark s at $p1_table[][12]
for ($i = 0; $i < $p1tablesize; $i++){
    if ($p1_table[$i][4] ne $p1_table[$i + 1][4]){
        $p1_table[$i + 1][6] = "s";
    }
}
$p1_table[0][6] = "s"; #例外処理

#mark P2 at sd & us
for ($i1 = 0; $i1 < $p1tablesize; $i1++){
    if ($p1_table[$i1][6] eq "s"){# 1 TSS in a cluster
        if ($p1_table[$i1][4] ne $p1_table[$i1 + 1][4]){# single TSS in a cluster, SS >> P1,P2,CP
            $p0serialtemp = $p1_table[$i1][5];
            $p0_table[$p0serialtemp][9] = "P2"; $p0_table[$p0serialtemp][11] = "CP";
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
            #何もしない
        }
        elsif ($p1_table[$j1 + $j2][6] eq "u"){# eeeu
            #何もしない
        }
        elsif ($p1_table[$j1 - 1][6] eq "s" or $p1_table[$j1 - 1][6] eq "u"){# s/u_eee
            if ($p1_table[$j1 + $j2][6] eq "d" or $p1_table[$j1 + $j2][6] eq "s"){# s/u_eee  s/d
                $mid = 0;
                $mid = $j1 - 1 + int ($j2 / 2);# $j2が奇数なら真ん中、偶数なら真ん中二つの北側
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
        $p0_table[$p0serialtemp][11] = "CP";
    }
}

print "P2 marking done\n";

##Second chop 2

#extract P2 peaks
$p2serialtemp = 0;
for ($i = 0; $i < $p0tablesize; $i++){
    if ($p0_table[$i][9] eq "P2"){
        $p2_table[$p2serialtemp][0] = $p0_table[$i][0]; $p2_table[$p2serialtemp][1] = $p0_table[$i][1]; $p2_table[$p2serialtemp][2] = $p0_table[$i][2]; $p2_table[$p2serialtemp][3] = $p0_table[$i][3]; $p2_table[$p2serialtemp][4] = $p0_table[$i][4]; $p2_table[$p2serialtemp][5] = $p0_table[$i][5]; $p2_table[$p2serialtemp][6] = $p0_table[$i][6]; $p2_table[$p2serialtemp][7] = $p0_table[$i][7]; $p2_table[$p2serialtemp][8] = $p0_table[$i][8]; $p2_table[$p2serialtemp][9] = $p0_table[$i][9]; $p2_table[$p2serialtemp][11] = $p0_table[$i][11];
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
        $p1first_u = 0; $p1last_d = 0; $p1cn = 0; $p1cs = 0;#初期化、スコープを意識
        if ($p2_table[$p2 + 1][1] - $p2_table[$p2][1] > $gapdistance2){ # distance between P2 peaks
            print "cut at P2_$p2\n";
            #作り直し　$p1_table[][12]のs/u/d/eを見て谷を決める。
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


#mark CP
for ($i = 0; $i < $p2tablesize; $i++){#clear p2_table[][]
    @{$p2_table[$i]} = "";
}
$p2serialtemp = 0;
for ($i = 0; $i < $p0tablesize; $i++){#clusterID変わったのでp2_table[][]作り直し
    if ($p0_table[$i][9] eq "P2"){
        $p2_table[$p2serialtemp][0] = $p0_table[$i][0]; $p2_table[$p2serialtemp][1] = $p0_table[$i][1]; $p2_table[$p2serialtemp][2] = $p0_table[$i][2]; $p2_table[$p2serialtemp][3] = $p0_table[$i][3]; $p2_table[$p2serialtemp][4] = $p0_table[$i][4]; $p2_table[$p2serialtemp][5] = $p0_table[$i][5]; $p2_table[$p2serialtemp][6] = $p0_table[$i][6]; $p2_table[$p2serialtemp][7] = $p0_table[$i][7]; $p2_table[$p2serialtemp][8] = $p0_table[$i][8]; $p2_table[$p2serialtemp][9] = $p0_table[$i][9]; $p2_table[$p2serialtemp][10] = $p0_table[$i][10]; $p2_table[$p2serialtemp][11] = $p0_table[$i][11];
        $p2serialtemp++;
    }
}

#select CP from P2
$p2a = 0; $p2b = 0;
while ($p2a < $p2tablesize){
    $p2b = 0; @p2_height = ""; %p2height2p0serial ="";
    while ($p2_table[$p2a][4] eq $p2_table[$p2a + $p2b][4]){
        push (@p2_height, $p2_table[$p2a + $p2b][3]);
        $p2height2p0serial{$p2_table[$p2a + $p2b][3]} = $p2_table[$p2a + $p2b][5]; #height >> p0serial
        $p2b++;
    }
    #$p2b = p2# in a cluster
    @p2_height = sort {$b <=> $a} @p2_height;
    $p0serial_p2highest = $p2height2p0serial{$p2_height[0]};
    $p0_table[$p0serial_p2highest][11] = "CP";
    $p2a = $p2a + $p2b;
}




#output file = @tablesplit
for ($i = 0; $i < $p0tablesize; $i++){
    for ($j = 0; $j < 12; $j++){
        print OUT "$p0_table[$i][$j]\t";
    }
    print OUT "\n";
}
close (OUT);

#output file = @tablesplit

