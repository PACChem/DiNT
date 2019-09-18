#!/usr/bin/perl

$mm=1; # MM PES
#$mm=0; # TB PES

@do1 = ("c1");
@heavy = (4);
@do2 = ("N2");
@do3 = (300);

$cwd = qx ! pwd !;
chomp($cwd);

foreach $i1 (@do1) { $nheavy = $heavy[$hindex++]; foreach $i2 (@do2) { foreach $i3 (@do3) {
$dir = $i1.".".$i2.".".$i3;
$dir1 = $i1.".".$do2[0].".".$do3[0];
$dir2 = $i1.".".$i2.".".$do3[0];
print "Doing $dir \n";
if (-e "$dir") { print "Directory $dir already exists.\n";} else {
qx ! cp -R script $dir !;
chdir "$cwd/$dir";
open(IN,"<../script/dintmom.script");  # ALWAYS USE TEMPLATE IN SCRIPT/ ?
open(OUT,">dintmom.pl");
foreach $line (<IN>) {
$line =~ s/SCRIPT1/c1/;
$line =~ s/SCRIPT2/$i2/;
$line =~ s/SCRIPT3/$i3/;
$line =~ s/SCRIPT4/$nheavy/;
print OUT $line;
}
qx ! cp ../c1.xyz .!;
if ($mm == 1) {
qx ! cp ../c1.tinker opt/tinker.xyz!;
qx ! cp ../c1.tinker samp/r0/tinker.xyz!;
qx ! cp ../c1.tinker r0/tinker.xyz!;
qx ! cp ../c1.key opt/tinker.key!;
qx ! cp ../c1.key samp/r0/tinker.key!;
qx ! cp ../c1.key r0/tinker.key!;
}
if (-e "../amfit.dat") { qx ! cp ../amfit.dat r0/. !};
if (-e "../$i1.basis") { qx ! cp ../$i1.basis r0/basis.dat !};
if (-e "../$i1.coef") { qx ! cp ../$i1.coef r0/coef.dat !};
if ($i3 ne @do3[0] ) { qx !cp $cwd/$dir2/dint.lj $cwd/$dir/.! ; }
if ($i3 ne @do3[0] || $i2 ne $do2[0]) {qx !cp $cwd/$dir1/dint.brot $cwd/$dir/.! ;
		     qx !cp $cwd/$dir1/dint.geo $cwd/$dir/.! ;
		     qx !cp $cwd/$dir1/r0/*.dat $cwd/$dir/r0/.!	}; #reuse target-specific info
qx ! chmod u+x dintmom.pl !;
chdir "$cwd";
}
chdir "$cwd/$dir";
qx ! perl ./dintmom.pl !;
qx ! cat dint.summary >> ..\/all.summary !;
chdir "$cwd";
}}}
