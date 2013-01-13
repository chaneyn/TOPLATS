use strict;

#Script to control changes in the TOPLATS code

#LOCAL VARIABLES
my $FLAG;
my $i;
my $GENERAL_FILE = "/home/ice/nchaney/PROJECTS/TOPLATS_DEVELOPMENT/DATA/LittleRiver/GLOBAL_PARAMETER_TEST.txt";

#EXECUTABLES
my $TOPLATS_NEW = "./TOPLATS_3.0";
my $UNIT_TESTS = "./TESTS_DRIVER";
my $CATCHMENT_TESTS = "./CATCHMENT_TESTS_DRIVER";

print "Running the model\n";
system("$TOPLATS_NEW $GENERAL_FILE");

print "Running the catchment tests\n";
system("$CATCHMENT_TESTS");
