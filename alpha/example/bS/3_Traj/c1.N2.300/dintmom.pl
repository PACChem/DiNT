#!/usr/bin/perl

# TO DO:
# 1. Improve error catching
# 2. Improve generality (ignore input case, add atom types, etc)
# 3. Parallelize xyz sampling?
# 4. Add diatomic baths (sampling over their coordinates and extra lines to dinttraj.pl)


####################
# SET THESE
####################
$exe         = "rohMP";
#$exe         = "tbplusexp6-2";
$ntraj       =  100;            # Per core; total = 32*NTRAJ = 32*NTRAJ
$target      = "c1" ;     # Name of the target
$bath        = "N2" ;     # Name of the bath gas, case sensitive
$temperature = 300 ;       # Temperature
$geofile = $target.".xyz";  # ESTOK-formatted geometry, in the launch directory
$energykcal = 40.; 	    # Initial energy (kcal/mol)
# Guess BMAX and COM 
# From PROCI paper: bmax = 9-12 A for C1-C8(+He), +1 for (+!He)
# Override below and set bmax and com, if necesseary
$nheavy = 4;	    # Number of heavy atoms, used only for BMAX and COM
$bmax = 6. + $nheavy/2.;    # Maximum impact parameter
if ($bath ne "He") {$bmax++};
$com = $bmax + 1.;          # Initial and final center of mass distance
####################

####################
# RESTART
# on = always run, write new dint.XXXX
# off = never run, read existing dint.XXXX
# check = run if dint.XXXX cannot be found, read dint.XXXX otherwise
####################
#$optonoff = "on";    # optimize geometry, write dint.geo
$optonoff = "off";
$optonoff = "check";
#$ljonoff = "on";     # calculate LJ parameters, write dint.lj
$ljonoff = "off";
$ljonoff = "check";
#$samponoff = "on";   # sample unimolecular target, write dint.samp
$samponoff = "off";
$samponoff = "check";
#$trajonoff = "on";   # run collision trajectories, write dint.traj
$trajonoff = "off";
$trajonoff = "check";
#$momonoff = "on";    # calculate moments, write dint.mom
$momonoff = "off";
$momonoff = "check";
$cleanonoff = "on";   # clean up unnecessary files
$cleanonoff = "off";
####################

# CONSTANTS
$autoev=27.2113961;
$autocmi=219474.63067;
$autokcal=627.5095;
$energy = $energykcal/$autokcal*$autoev;

# GET CURRENT DIRECTORY
$cwd = qx ! pwd !;
chomp($cwd);

# SET EXE
qx ! cp -f r0/dint-$exe.x.opt r0/dint.x.opt !;
qx ! cp -f lj/auto1dmin-$exe.x lj/auto1dmin.x !;

# OPTIMIZE GEOMETRY IN ./opt
if ($optonoff eq "on" || ($optonoff eq "check" && (!(-e "dint.geo")))) {
chdir "$cwd/opt" ;
if (-e "../r0/amfit.dat") { qx ! cp ../r0/amfit.dat . !}
if (-e "../r0/coef.dat") { qx ! cp ../r0/coef.dat . !}
if (-e "../r0/basis.dat") { qx ! cp ../r0/basis.dat . !}
qx ! ./dintopt.pl $geofile !; # If successful, this will generate opt/dint.geo with the optimized geometry and energy
qx ! cp dint.geo ..\/. !;
chdir "$cwd" ;
if (-e "dint.geo") { print " Geometry optimized! \n";} else {die" Geometry optimization failed -- no dint.geo file was found!\n";}
} else {
if (-e "dint.geo") { print " Geometry optimization skipped! \n";} else {die" Geometry optimization skipped but no dint.geo file was found!\n";}
}

# SAMPLE INITIAL COORDINATES OF THE UNIMOLECUALR SPECIES
if ($samponoff eq "on" || ($samponoff eq "check" && (!(-e "r0/xyz.dat")))) {
chdir "$cwd/samp" ;
if (-e "../r0/amfit.dat") { qx ! cp ../r0/amfit.dat r0/. !}
if (-e "../r0/coef.dat") { qx ! cp ../r0/coef.dat r0/. !}
if (-e "../r0/basis.dat") { qx ! cp ../r0/basis.dat r0/. !}
qx ! ./dintsamp.pl $energy !; # If successful, this will generate sampled coordinates for the unimolecular species at the target energy
chdir "$cwd" ;
} else {
if (-e "r0/xyz.dat") { print " Unimolecular species coordinates sampling skipped! \n";} else {die" Coordinate sampling skipped but no xyz.dat file was found! \n";}
}
if (-e "dint.brot") { print "Found dint.brot\n"; } else {
chdir "$cwd/samp" ;
if (-e "fort.80") {qx ! ./brot.x !;} else { die("Couldn't find sampled fort.80\n") }
qx ! cp dint.brot ..\/. !;
chdir "$cwd" ;
}

# GET LENNARD-JONES PARAMETERS
if ($ljonoff eq "on" || ($ljonoff eq "check" && (!(-e "dint.lj")))) {
open(LJ,">lj/input");
print LJ "$geofile\n$bath\n";
close(LJ);
if (-e "r0/amfit.dat") { qx ! cp r0/amfit.dat . !}
if (-e "r0/coef.dat") { qx ! cp r0/coef.dat . !}
if (-e "r0/basis.dat") { qx ! cp r0/basis.dat . !}
qx ! ./lj/auto1dmin.x < lj/input !;
qx ! mv lj.dat dint.lj !;
if (-e "dint.lj") { print " Lennard-Jones parameters calculated! \n";} else {die" Lennard-Jones calculation failed -- no dint.lj file was found! \n";}
} else {
if (-e "dint.lj") { print " Lennard-Jones parameters calculation skipped! \n";} else {die" Lennard-Jones calculation skipped but no dint.lj file was found! \n";}
}

# RUN COLLISION TRAJECTORIES
if ($trajonoff eq "on" || ($trajonoff eq "check" && (!(-e "dint.traj")))) {
qx ! cp r0/dinttraj.pl .!;
$tmp = qx ! ./dinttraj.pl $temperature $energy $target $bath $com $bmax $ntraj !;
print $tmp;
if (-e "dint.traj") { print " Found dint.traj ! \n";} else {die" No dint.traj file was found! \n";}
} else {
if (-e "dint.traj") { print " Trajectory calculation skipped! \n";} else {die" Trajectory calculation skipped but no dint.traj file was found! \n";}
}

# CALCULATE MOMENTS
open(IN,"<dint.lj");
@dump = (<IN>);
close(IN);
$tmp = "@dump";
@dump = split/ +/,$tmp;
$ep = $dump[6];
$si = $dump[9];
$m1 = $dump[12];
$m2 = $dump[13];
chomp($m2);
if ($momonoff eq "on" || ($momonoff eq "check" && (!(-e "dint.mom")))) {
open(INCLUDE,">mom/param.inc");
print INCLUDE "       ep=$ep\n";
print INCLUDE "       si=$si\n";
print INCLUDE "       mu1=$m1\n";
print INCLUDE "       mu2=$m2\n";
close(INCLUDE);
qx ! gfortran -O3 mom/moments2.f mom/zonz.f -o mom/mom2.x !;
$tmp = qx ! mom/mom2.x < dint.traj !;
open(OUT,">mom/mom.out");
print OUT $tmp;
close(OUT);
} else {
if (-e "dint.mom") { print " Moments calculation skipped! \n";} else {die" Moments calculation skipped but no dint.mom file was found! \n";}
}

open(IN,"<dint.mom");
foreach $line (<IN>) {
if ($line =~ /alpha =/) {chomp($line);@dat=split/ +/,$line; $alpha = $dat[4]." ".$dat[5]};
if ($line =~ /LJ/) {chomp($line);@dat=split/ +/,$line; $zlj = $dat[4]};
}

open(SUM,">dint.summary");
print SUM "$target $bath $temperature $si $ep $zlj $alpha\n";
print "$target $bath $temperature $si $ep $zlj $alpha\n";

# CLEAN UP
if ($cleanonoff eq "on") { qx !./clean.x! }

