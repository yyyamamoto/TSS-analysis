TSS-analysis.  
clustering TSS tags in genomics.  
TSS-seq data are mapped into a genome, simlified form tag read to position of the tags, and grouped into TSS clusters with a width around tens of bps. From clustered data, peak position, positions of both edges, and tag number of the cluster are to be identified.


TSSClusteringFromMappedData.pl

  2017.01.31
  Yoshiharu Y. Yamamoto, Faculty of Applied Biological Sciences, Gifu University
  yyy@gifu-u.ac.jp >> yamamoto.yoshiharu.h5@f.gifu-u.ac.jp (2024.12.13)
  Solexa TSS analysis
  Preparation of TSS clusers from mapped data
  rule 1: cut if distance is eq or more than 20 bp (=$gapdistance1)
  rule 2: cut if distance between secondary peaks is eq or more than 100 bp (=$gapdistance2)
  cutting site: bottom of P1 peaks, longest space among TSS in the bottom
  ref: Tokizawa M, Kusunoki K, Koyama H, Kurotani A, Sakurai T, Suzuki Y, Kurata T, Yamamoto YY. "Identification of Arabidopsis genic and non-genic promoters by pair-end sequencing of TSS tags". Plant J, 89: 671-681, 2017.

  argument 1: input file name
  input file: mapped TSS data (tab delimited table)
  1st column: chromosome number (Chr01~)
  2nd column: chromosomal position
  3rd column: direction (F or R)
  4th column: number of TSS tag
 Table should be sorted according to 1)chromosome number, 2)direction, and 3)chromosomal position in this order.

  argument 2: output file name
  output file:
  1st column: chromosome number
  2nd column: chromosomal position
  3rd column: direction (F or R)
  4th column: number of TSS tag
  5th column: cluster ID (temporal, rename recommended)
  6th column: serial number of TSS
  7th column: u/s/e/d mark (ignore this)
  8th column: P1 peak mark
  9th column: P1 serial number
  10th column: P2 peak mark
  11th column: P2 serial number
  12th column: CP peak mark (=cluster peak)

2024.12.13
Yoshiharu Y. Yamamoto, yamamoto.yoshiharu.h5@f.gifu-u.ac.jp
new script uploaded: TSSClusteringFromMappedData_1_2024.pl
A bug of no marking of new Cluster Peak after final cutting of a long cluster, is fixed.
Used in our recent report: "Genetic stability of genic and non-genic promoters identified by comparative TSS-seq of Arabidopsis natural variations", Sugekawa K, Ezeh OS, Yamamoto YY (manuscript in preparation)
new script uploaded: TSSClusteringFromMappedData_2_prepClusterTable2018.pl
Apply this to the output file from TSSClusteringFromMappedData_1_2024.pl

How to use  25.8.16
1
Prep TSS data like this file: TSSclusteringfromMappedData_inputtest.tab
memo: Consider Cap-Singature (first G or GG in the TSS tag seq) when mapping TSS tags into genome seq. If the first G of a tag matches with the genome seq, there are two possibilities: the G is artificially added, or it comes from the corresponding genome seq. I recommend to assume the latter situation first, and the assumption would not cause any critical problems for many analyses. If you are trying to analyze minor promoter sections or events, you may recall this matter. See Fig S3 & 4 in Tokizawa et al. Plant J 90: 587, 2017 more detail.

2
Use TSSClusteringFromMappedData.pl
There are two arguments as described in the top of the script.
The second argument is the name of the output file like this:
TSSclusteringfromMappedData_outputtest.tab

3
Use TSSClusteringFromMappedData_1_2024


