#!/usr/bin/perl

$cpupernode=36;

$eee = 4.5;
$ttt = 2000;
$xxx = "zz";

if ($ARGV[0] ne "") { $xxx = $ARGV[0] };

if ($xxx eq "za") { $eee = 5.35; $ttt = 1000; }

$ete = $eee/8.61692E-05/3;

$r1 = 1;
$r2 = $cpupernode;
$"="";
$pre = "sub$xxx";
$exe= "dint-h3o.x.opt";

for ($i=$r1;$i<=$r2;$i++) { @dirs=(@dirs,"$xxx$i"); }

open(SUB,">$pre-$r1.x");
print SUB "cd $dirs[0]\n";
foreach $dd (@dirs) {
	$count++;
	if ($count == 1) {
		$job++;
		open(SUB,">$pre-$job.x");
		print SUB "cd $dirs[0]\n";
	}
	if (-e $dd) {die("directory $dd already exists!\n")};
	mkdir $dd;
	qx ! cp r0/$exe $dd/. !;
	qx ! cp r0/fort.70 $dd/. !;
	qx ! cp r0/basis.dat $dd/. !;
	qx ! cp r0/coef.dat $dd/. !;
	open(IN,"<r0/input");
	open(OUT,">$dd/input");
	$ran = int(rand(1000000000));
	foreach $line (<IN>) {
		$line =~ s/RANSEED/$ran/;
		$line =~ s/TTT/$ttt/g;
		$line =~ s/ETE/$ete/g;
		$line =~ s/EEE/$eee/g;
		print OUT $line;
	}
	close(IN);
	close(OUT);
	print SUB "cd ../$dd\n";
	print SUB "time ./$exe < input > output &\n";
	if ($count == $cpupernode) {
		$count=0;
		print SUB "wait\n";
		close(SUB);
		qx ! chmod u+x $pre-$job.x !;
		qx ! qqbb q ./$pre-$job.x !;
	}
}

