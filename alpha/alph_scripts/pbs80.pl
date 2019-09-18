#!/usr/bin/perl

$wall="20:00:00";
$cpuspernode = 1;

if ($#ARGV < 0 || $#ARGV > 1) {die("usage: energies80.pl inputfile.80 [outputfile]\n")};

if (-e "zeros.dat") { open(ZEROS,"<zeros.dat"); foreach $line (<ZEROS>) {chomp($line); $zero[++$k]=$line; print "z $k $zero[$k]\n";}}

$cwd = qx ! pwd !;
chomp($cwd);
$autocmi = 219474.63067;

if ($#ARGV == 0) {$output = ""};
if ($#ARGV == 1) {$output = $ARGV[1]};
print " Reading $ARGV[0] and writing $output \n";

$" = "";
open(QC,"<qc.mol"); @tmp = (<QC>); $qc = "@tmp" ; close(QC);
open(IN,"<$ARGV[0]"); @dump = (<IN>); close(IN);
open(AI,">ai.dat");
open(OUTXX,">$output"); 

#read geo file, count geometries, and store input files
$k=1;
$natom = $dump[0];
while ($i < $#dump) {
if ($natom != $dump[$i]) {die("number of atoms changed from $natom to $dump[$i]\n");}
$geo[$k] = "";
for ($j=0;$j<$natom;$j++) { $geo[$k].=$dump[$i+$j+2] ; }
chomp($geo[$k]);
@dat = split/ +/,$dump[$i+1]; $rr[$k] = $dat[2];
$qc0[$k] = $qc;
$qc0[$k] =~ s/GEOMETRY/$geo[$k]/;
$k++;
$i+=($natom+2);
}
$geoms=$k-1;

#look for input and output files
for ($k=1;$k<=$geoms;$k++) {
if (!(-e "qc$k.in")) { 
print " Writing input file for $k th geometry\n";
open(OUT,">qc$k.in"); print OUT $qc0[$k]; close(OUT);
}
if (!(-e "qc$k.out")) {
@todolist = (@todolist,"qc$k.in");
} else {
open(MOLOUT,"<qc$k.out");
foreach $line (<MOLOUT>) {
$line =~ s/D/E/;
if ($line =~ /SETTING ENERGY([0-9]+) += +([E\+0-9\.\-]+)/) {
$energy[$1][$k] = ($2-$zero[$1])*$autocmi;
if ($1 > $maxEind) { $maxEind = $1 }
}
}
print "$k $rr[$k]";
for ($z=1;$z<=$maxEind;$z++) { printf("%15.5f",$energy[$z][$k]); }
print "\n";

print AI "$k $rr[$k]";
for ($z=1;$z<=$maxEind;$z++) { printf AI ("%15.5f",$energy[$z][$k]); }
print AI "\n";

print OUTXX "$natom";
print OUTXX " $k $rr[$k]";
for ($z=1;$z<=$maxEind;$z++) { printf OUTXX ("%15.5f",$energy[$z][$k]); }
print OUTXX "\n";
print OUTXX "$geo[$k]\n";
}
}
$needed = $#todolist+1;

print " $geoms geometries found in $ARGV[0]\n";
print " $needed output files missing\n";

$node=int($needed/$cpuspernode)+1;

print " Need $node $cpuspernode-core nodes\n";

if ($needed > 0) {
for($i=1;$i<=$node;$i++) {
open(OUT,">job$i.x");
print OUT <<EOF;
#!/bin/sh
#
#PBS -N job$i
#PBS -l nodes=1:ppn=16
#PBS -l walltime=$wall
#PBS -j oe

mkdir -p /scratch/ajasper
cd \$PBS_O_WORKDIR

EOF

for($j=0;$j<$cpuspernode;$j++) {
if ($todolist[$l] ne "") {print OUT "molpro -n 1 -d /home/ajasper/scratch $todolist[$l] & \n";}
$l++;
}
print OUT "wait\n";
close(OUT);
qx!chmod u+x job$i.x!;
qx!qsub job$i.x!;
}
}



