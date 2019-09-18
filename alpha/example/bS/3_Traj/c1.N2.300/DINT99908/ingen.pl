#!/usr/bin/perl

$suf="a";

$cpupernode=32;

$r1 = $ARGV[0];
$r2 = $ARGV[1];
$"="";
$pre = "sub";
$exe = "dint.x.opt";

for ($i=$r1;$i<=$r2;$i++) { @dirs=(@dirs,"$suf$i"); }

foreach $dd (@dirs) {
	$count++;
	if ($count == 1) {
		$job++;
		open(SUB,">$pre-$job.x");
print SUB "set OMP_THREAD_LIMIT=1\nexport OMP_THREAD_LIMIT=1\n";
		print SUB "cd $dirs[0]\n";
	}
	if (-e $dd) {die("directory $dd already exists!\n")};
	mkdir $dd;
	qx ! cp r0/$exe $dd/. !;
	qx ! cp r0/tinker.xyz $dd/. !;
	qx ! cp r0/tinker.key $dd/. !;
	qx ! cp r0/amfit.dat $dd/. !;
	qx ! cp r0/coef.dat $dd/. !;
	qx ! cp r0/basis.dat $dd/. !;
	open(IN,"<r0/input");
	open(OUT,">$dd/input");
	$ran = int(rand(1000000000));
	foreach $line (<IN>) {
		$line =~ s/RANSEED/$ran/;
		print OUT $line;
	}
	close(IN);
	close(OUT);
	print SUB "cd ../$dd\n";
	print SUB "time ./$exe < input >  output &\n";
	if ($count == $cpupernode) {
		$count=0;
		print SUB "wait\n";
		close(SUB);
		qx ! chmod u+x $pre-$job.x !;
		qx ! qq4 ./$pre-$job.x $pre-$job !;
	}
}
qx ! chmod u+x $pre-$job.x !;
