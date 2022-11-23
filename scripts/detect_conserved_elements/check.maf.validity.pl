use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

##############################################################
##############################################################

sub printHelp{
    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script checks if the MAF file is valid.\n";
    print "\n";
    
    print "Options:\n";
    
    foreach my $par (@{$parnames}){
	print "--".$par."  [  default value: ".$parvalues->{$par}."  ]\n";
    }
    print "\n";
}

##############################################################
##############################################################

my %parameters;

$parameters{"pathMAF"}="NA";

my @defaultpars=("pathMAF");
my %defaultvalues;

foreach my $par (keys %parameters){
    $defaultvalues{$par}=$parameters{$par};
}

## update arguments

my $nbargs=@ARGV;

for(my $i=0;$i<$nbargs; $i++){
    my $arg=$ARGV[$i];
    $arg=substr $arg,2;
    
    my @s=split("=",$arg);
    my $parname=$s[0];
    my $parval=$s[1];
    
    if(exists $parameters{$parname}){
	$parameters{$parname}=$parval;
    }
    else{
	print "Error: parameter ".$parname." was not recognized!!!\n";
	printHelp(\@defaultpars, \%defaultvalues);
	exit(1);
    }
}

## show parameters

print "\n";

print "Running program with the following parameters:\n";

foreach my $par (@defaultpars){
    print "--".$par."=".$parameters{$par}."\n";
}

print "\n";

##############################################################
##############################################################

print "Reading MAF and checking sequences...\n";

my $nbdone=0;

my $input;
my $pathmaf=$parameters{"pathMAF"};
my @s=split("\\.");
my $ext=$s[-1];

if($ext eq "gz"){
    open($input, "zcat $pathmaf |");
} else{
    open($input, $pathmaf);
}

my $line=<$input>;

while($line){
    my $prefix=substr $line, 0, 1;

    if($prefix eq "s"){
	chomp $line;
	$line =~ tr/ //s;

	my @s=split("\t", $line);
	my $nbfields=@s;

	if($nbfields!=7){
	    print "Found ".$nbfields." fields in alignment line ".$line."\n";
	    print "MAF file invalid\n";
	    exit(1);
	}
	
	my $size=$s[3]+0;
	my $sequence=$s[6];
	$sequence =~ s/-//gi;
	my $len=length $sequence;

	if($size != $len){
	    print "Incompatible declared and real size in alignment line ".$line.".\n";
	    print "MAF file invalid\n";
	    exit(1);
	}
    }

    if($prefix eq "a"){
	$nbdone++;

	if($nbdone%10000==0){
	    print $nbdone." alignments done.\n";
	}
    }
    
    $line=<$input>;
}

close($input);

print "Everything ok.\n";
print "Done.\n";

##############################################################
