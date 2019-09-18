#!/usr/bin/perl

$temperature = $ARGV[0];
$energy = $ARGV[1];
$target = $ARGV[2];
$bath = $ARGV[3];
$com = $ARGV[4];
$bmax = $ARGV[5];
$ntraj = $ARGV[6];

open(OPT,"<dint.geo");
$natom=-2;
foreach $line (<OPT>) {
$natom++;
if ($natom == -1) {$zero = $line; chomp($zero)}
if ($natom ==  0) {chomp($line); @b = split/ +/,$line; shift @b; @b = sort @b; if (($b[2]-$b[1]) < ($b[1]-$b[0])) { $b2 = $b[0]; $b1 = ($b[1]+$b[2])/2.} else {$b2 = $b[2]; $b1 = ($b[1]+$b[0])/2.}; print "DEBUG BROT: @b symmetrized to $b1 $b2 \n"};
if ($natom > 0) {$geo = $geo.$line}
}
chomp($geo);
close(OPT);

if (-e "dint.brot") {
open(BROT,"<dint.brot");
$line=<BROT>;
@dat=split/ +/,$line;
$b1=$dat[1];
$b2=$dat[2];
print "DEBUG BROT: using $b1 and $b2 instead \n";
}

$exe = "dint.x.opt";
$nsub=20;

open(IN,"<r0/input.template");
open(OUT,">r0/input");
foreach $line (<IN>) {
        $line =~ s/ZERO/$zero/;
        $line =~ s/NATOM/$natom/;
        $line =~ s/GEOMETRY/$geo/;
        $line =~ s/TEMP/$temperature/g;
        $line =~ s/ENERGY/$energy/;
        $line =~ s/BROT1/$b1/;
        $line =~ s/BROT2/$b2/;
        $line =~ s/COM/$com/;
        $line =~ s/BMAX/$bmax/;
        $line =~ s/NTRAJ/$ntraj/;
        if ($line =~ /BATH/) {
                if ($bath eq "He") {print OUT "1 0 1 0 0.\n He 4.002603 0. 0. 0.\n";}
                if ($bath eq "Ar") {print OUT "1 0 1 0 0.\n Ar 39.962383 0. 0. 0.\n";}
                if ($bath eq "Kr") {print OUT "1 0 1 0 0.\n Kr 83.911507 0. 0. 0.\n";}
                if ($bath eq "N2") {print OUT "2 6 -1 0 0.\n 20000 F n2x.$temperature n2p.$temperature F\nN2 14.00307\nN2 14.00307\n";}
        } else {
        print OUT $line;
        }
}

@tmp = qx!ls!;
foreach $i (@tmp) {if ($i =~ /DINT/) {$j++} }
if ($j == 0) {
for ($i=1;$i<=$nsub;$i++) { qx ! ./qqdint !; }
} elsif ($j != $nsub) {
die("Found $j DINT#### directories instead of $nsub\n");
} else {
$tmp = $nsub*32*$ntraj;
qx !echo "$tmp 1 $bmax $target $bath $temperature" > dint.traj!;
qx !cat D*/*/fort.31 >> dint.traj!;
}
