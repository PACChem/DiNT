#!/usr/bin/perl

$geofile = $ARGV[0];
if ($geofile eq "") {die("DINTOPT.PL DIED -- USAGE: dintopt.pl GEO.FILE \n")}

open(IN,"<input.template");
open(GEO,"../$geofile");

foreach $line (<GEO>) {
if ($line =~ /[A-Za-z0-9]/) {$natom++};
$line =~ s/C /C 12./;
$line =~ s/H /H 1.007825/;
$line =~ s/O /O 15.9949/;
$line =~ s/Ho /Ho 1.007825/;
$line =~ s/N /N 14.003074/;
$geo = $geo.$line;
}
chomp($geo);

open(OUT,">input");
foreach $line (<IN>) {
$line =~ s/NATOM/$natom/;
$line =~ s/GEOMETRY/$geo/;
print OUT $line;
}

qx ! ../r0/dint.x.opt < input > output !;
