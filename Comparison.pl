use strict;

#Script to control changes in the TOPLATS code

#INPUT VARIABLES
my $OLD_DATA_ROOT = "/home/ice/nchaney/PROJECTS/TOPLATS_DEVELOPMENT/DATA/LittleRiver/OLD";
my $NEW_DATA_ROOT = "/home/ice/nchaney/PROJECTS/TOPLATS_DEVELOPMENT/DATA/LittleRiver/NEW";
#my @FLAG_VARS = ("EF","CWB","PRI","ET","STZB","WTB","SSR","RMEC","SWE");
my @FLAG_VARS = ();

#LOCAL VARIABLES
my $FLAG;
my $i;

#EXECUTABLES
#my $TOPLATS_OLD = "../TOPLATS_3.33/TOPLATS_3.0";
my $TOPLATS_NEW = "./TOPLATS_3.0";

#Recompile the newest TOPLATS version
#print "Recompiling the newest version of TOPLATS\n";
#chdir("../TOPLATS_2.01");
#system("make > log.txt");

#1. Run both old and newer TOPLATS
print "Running the old version\n";
#system("$TOPLATS_OLD ../DATA/LittleRiver/OLD/toplats.io.files > log1.txt");
print "Running the new version\n";
system("$TOPLATS_NEW /home/ice/nchaney/PROJECTS/TOPLATS_DEVELOPMENT/DATA/LittleRiver/NEW/toplats.io.files");
#print "$TOPLATS_NEW ../DATA/LittleRiver/NEW/toplats.io.files\n";
#system("valgrind $TOPLATS_NEW ../DATA/LittleRiver/NEW/toplats.io.files");

#2. Compare output files to ensure we are obtaining the same results
for ($i=0;$i<@FLAG_VARS;$i++){
	$FLAG = `diff $OLD_DATA_ROOT/$FLAG_VARS[$i].txt $NEW_DATA_ROOT/$FLAG_VARS[$i].txt`; 
	if ($FLAG){print "DIFFERENCES IN $FLAG_VARS[$i].txt\n"};
	}
