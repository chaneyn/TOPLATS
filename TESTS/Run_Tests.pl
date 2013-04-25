use strict;

#Script to control changes in the TOPLATS code

#LOCAL VARIABLES
my $FLAG;
my $i;
my $GENERAL_FILE = "TESTS/INPUT/GENERAL_PARAMETER.txt";
my $NEW_REGIONAL_FILE = "TESTS/OUTPUT/Regional_Variables.txt";
my $OLD_REGIONAL_FILE = "TESTS/OUTPUT_REFERENCE/Regional_Variables.txt";

#EXECUTABLES
my $TOPLATS_NEW = "./TOPLATS_3.0";
my $UNIT_TESTS = "./TESTS/TESTS_DRIVER";
my $CATCHMENT_TESTS = "./TESTS/CATCHMENT_TESTS_DRIVER";

chdir("../");

print "Running the unit tests\n";
#print $UNIT_TESTS."\n";
system("$UNIT_TESTS");

#print "Running the model\n";
system("$TOPLATS_NEW $GENERAL_FILE");

#print "Running the catchment tests\n";
system("$CATCHMENT_TESTS $NEW_REGIONAL_FILE $OLD_REGIONAL_FILE");
