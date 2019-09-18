#!/usr/bin/perl

qx! cat job.batch*/fort.80 job.?/fort.80 job.??/fort.80 job.???/fort.80 > fort.80 !;
qx! cat job*/minimum.xyz > tmp !;

open(IN,"<tmp");
@dump = (<IN>);
close(IN);

$natom = $dump[0];

while ($i < $#dump) {
if ($natom != $dump[$i]) {die("number of atoms changed from $natom to $dump[$i]\n");}
$natom = $dump[$i];
@dat = split/ +/,$dump[$i+1];
$energy = $dat[3];
if ($energy < $minimum) {
print "Finding minimum... $energy\n";
$minimum=$energy;
$save = "";
for ($j=0;$j<$natom+2;$j++) { $save.=$dump[$i+$j] ; }
}
$i+=($natom+2);
}

#if ($#tmp > 0 && !(-e minimum.xyz)) {
open(OUT,">minimum.xyz");
print OUT $save;
#}

print "\n";
print "Rcom(A)  ECC  EnoCC (cm-1) \n";
@dump = qx ! cat */output | grep ' 1       1' !;
foreach $line (@dump) {
	chomp($line);
	@dat = split/ +/,$line;
	if ($printcount++ < 100) {print "$dat[1] $dat[3] $dat[6] $dat[7]\n";}
#bin
$binhi=7.25;
$binlo=2.25;
$binss=.1;
$ind=int(($dat[3]-$binlo)/$binss);
$count[$ind]++;
$weight=(200./($dat[6]-$minimum+200.))**2;
if ($dat[6] <= 3000.) {$countin[$ind]++}
if ($dat[6] <= 3000.) {$wcountin[$ind]+=$weight}
$val[$ind]+=$dat[6];
$val2[$ind]+=$dat[6]**2;
if ($dat[6] > $max[$ind] || $max[$ind] eq "") {$max[$ind]=$dat[6]}
if ($dat[6] < $min[$ind] || $min[$ind] eq "") {$min[$ind]=$dat[6]}
if ($weight > $wmax[$ind] || $wmax[$ind] eq "") {$wmax[$ind]=$weight}
if ($weight < $wmin[$ind] || $wmin[$ind] eq "") {$wmin[$ind]=$weight}
#bin
}

print "\n";
for($i=0;$i<$#count;$i++) {
$x=$binlo+($i+.5)*$binss;
$a=$count[$i]+0.;
if($a == 0) {$a = 1;}
$aa=$countin[$i]+0.;
$ab=($wcountin[$i]+0.)/$a;
$ac=$wmax[$i]+0.;
$ad=$wmin[$i]+0.;
$b=$val[$i]/$a+0.;
$c=$val2[$i]/$a+0.;
$d=sqrt($c-$b**2);
$e=$min[$i];
$f=$max[$i];
print "$x $a $aa $ab $ac $ad $b $e $f $d\n";
}

print "\nFound ".($#dump+1)." energies \n\n";

for ($j=81;$j<=86;$j++) {
if (-e "fort.$j") {
	$i=0;
	$k=0;
	@dump = qx ! cat fort.$j !;
	while ($i < $#dump) {
		chomp($dump[$i]);
		$natom = $dump[$i];
		chomp($dump[$i+1]);
		@dat = split/ +/,$dump[$i+1];
		$r[$k]=$dat[2];
		$e[$k][$j-81]=$dat[3];
		$k++;
		$i+=$natom+2;
	}
}
}

for ($l=0;$l<$k;$l++) {
#for ($l=0;$l<10;$l++) {
print " $r[$l] ";
for ($j=0;$j<=5;$j++) {
	print " $e[$l][$j] ";
}
print "\n";
}

print "\nMinimum energy = $minimum\n";
