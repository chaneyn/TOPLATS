use strict;

#Script to control changes in the TOPLATS code

#LOCAL VARIABLES
my $FLAG;
my $i;

#EXECUTABLES
my $TOPLATS_NEW = "./TOPLATS_3.0";
my $UNIT_TESTS = "./TESTS_DRIVER";
my $CATCHMENT_TESTS = "./CATCHMENT_TESTS_DRIVER";

print "Running the model\n";
system("$TOPLATS_NEW /home/ice/nchaney/PROJECTS/TOPLATS_DEVELOPMENT/DATA/LittleRiver/NEW/toplats.io.files");

print "Running the catchment tests\n";
system("$CATCHMENT_TESTS");
