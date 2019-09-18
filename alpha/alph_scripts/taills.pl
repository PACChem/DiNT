#!/usr/bin/perl

$target = $ARGV[0];
$tmax = $target*1.3;

open(IN,"<summary");
@dump = (<IN>);
#print "@dump";

foreach $line (@dump) {
chomp($line);
@dat = split/ +/,$line;

#$emax[$i]=$dat[9];
$emax[$i]=$dat[6]+2.*$dat[7];
$dirs[$i]=$dat[1];
if ($emax[$i] < $tmax && $ithis == 0) {$ithis = $i};
print " $line     $emax[$i] \n";
$i++;
}

$a = -1;
$b = 1;

for ($i=$ithis+$a;$i<=$ithis+$b;$i++) {
$n++;
$avgx+=($emax[$i]-$tmax);
$avgy+=log($dirs[$i]);
print " $i $dirs[$i] $emax[$i] \n";
}
$avgx/=$n;
$avgy/=$n;

for ($i=$ithis+$a;$i<=$ithis+$b;$i++) {
$slopen += ($emax[$i]-$tmax-$avgx)*(log($dirs[$i])-$avgy);
$sloped += ($emax[$i]-$tmax-$avgx)**2;
}
$b = $avgy - $slopen/$sloped*$avgx;

print "\n";
print " target= $tmax \n";
print " m= ".(exp($b))." \n";
