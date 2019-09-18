#!/usr/bin/perl


$ndo = 16;

open(IN,"<$ARGV[0]");

@dump = (<IN>);
$tmp = "@dump";

$tmp =~ /(        aa=.*\n)  -/s;

@p = split/\n/,$1;

for ($i=0;$i<$ndo;$i++){ print "$p[$i]\n"; }

if ($ndo == 12) {
for ($i=0;$i<4;$i++){ print "$p[$i]\n"; }
}


