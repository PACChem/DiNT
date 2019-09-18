#!/usr/bin/perl

$energy = $ARGV[0];
if ($energy eq "") {die("DINTSAMP.PL DIED -- USAGE: dintsamp.pl energy-in-eV \n")}

open(IN,"<r0/input.template");
open(OPT,"../dint.geo");

$natom=-2;
foreach $line (<OPT>) {
$natom++;
if ($natom == -1) {$zero = $line; chomp($zero)}
if ($natom ==  0) { ; }
if ($natom > 0) {$geo = $geo.$line}
}
chomp($geo);

open(OUT,">r0/input");
foreach $line (<IN>) {
$line =~ s/ZERO/$zero/;
$line =~ s/ENERGY/$energy/;
$line =~ s/NATOM/$natom/;
$line =~ s/GEOMETRY/$geo/;
#$ran = int(rand(1000000000));
#$line =~ s/RANSEED/$ran/;

print OUT $line;
}

qx ! cp ../r0/dint.x.opt r0/. !;
if (!(-e "fort.80")) { 
print "Cant find samp/fort.80. Launching job.\n"; 
qx ! ./qqdint !; 
qx ! cp fort.80 ../r0/xyz.dat ; cp fort.81 ../r0/ppp.dat ! 
} else { 
print "Found samp/fort.80\n"; 
qx ! cp fort.80 ../r0/xyz.dat ; cp fort.81 ../r0/ppp.dat ! 
}
