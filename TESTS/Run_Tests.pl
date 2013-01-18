use strict;

#Script to control changes in the TOPLATS code

#LOCAL VARIABLES
my $FLAG;
my $i;
my $GENERAL_FILE = "INPUT/GENERAL_PARAMETER.txt";
my $NEW_REGIONAL_FILE = "OUTPUT/Regional_Variables.txt";
my $OLD_REGIONAL_FILE = "OUTPUT_REFERENCE/Regional_Variables.txt";

#EXECUTABLES
my $TOPLATS_NEW = "../TOPLATS_3.0";
my $UNIT_TESTS = "./TESTS_DRIVER";
my $CATCHMENT_TESTS = "./CATCHMENT_TESTS_DRIVER";

print "Running the model\n";
system("$TOPLATS_NEW $GENERAL_FILE");

print "Running the catchment tests\n";
system("$CATCHMENT_TESTS $NEW_REGIONAL_FILE $OLD_REGIONAL_FILE");

print "Running the unit tests\n";
system("$UNIT_TESTS");
