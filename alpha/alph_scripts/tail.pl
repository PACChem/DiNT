#!/usr/bin/perl

$target = $ARGV[0];
$tmax = $target*1.3;

@dirs = qx ! ls -d j??; ls -d j???; ls -d j????; ls -d j????? !;
@output = qx !  tail -q -n1 j??/summary; tail -q -n1 j???/summary; tail -q -n1 j????/summary; tail -q -n1 j?????/summary !;

for ($i=0;$i<=$#dirs;$i++) {
chomp($dirs[$i]);
chomp($output[$i]);
$dirs[$i] =~ s/j//;
print "$dirs[$i] $output[$i] \n";

@dat = split/ +/,$output[$i];
$emax[$i]=$dat[8];
$emax[$i]=$dat[5]+2.*$dat[6];
if ($emax[$i] < $tmax && $ithis == 0) {$ithis = $i};

}
print "\n";

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
