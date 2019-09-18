#!/usr/bin/perl

open(IN,"<$ARGV[0]");

$i=1;
$j=1;
foreach $line (<IN>) {
if ($j<0) {$j++; next}  # skip N first lines
$line = " $line";
if ($line =~ /\-\-\-/) {last}
# get R and Vfit
@dat = split/ +/,$line; $rr = $dat[2]; $vv = $dat[4];
# save unique R
if ($rr < $rrlast) {@rlist = @empty;$i++}; 
@rlist = (@rlist,$rr);
$rrlast=$rr;
$key=$rr."-".$i; $v{$key}=$vv;
}

foreach $rr (@rlist) {
print " $rr ";
for ($i=1;$i<7;$i++) {
$key=$rr."-".$i; $vv=$v{$key};
print " $vv ";
}
print "\n";
}

