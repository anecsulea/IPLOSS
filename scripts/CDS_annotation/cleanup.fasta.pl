#!/usr/bin/perl
use strict;

##############################################################

sub printHelp{

    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script modifies sequence names to remove spaces and replaces ambiguous nucleotides with Ns.\n";
    print "\n";
    
    print "Options:\n";
    
    foreach my $par (@{$parnames}){
	print "--".$par."  [  default value: ".$parvalues->{$par}."  ]\n";
    }
    print "\n";
}

##############################################################
##############################################################

## parameters 

my %parameters;
$parameters{"pathInput"}="NA";
$parameters{"pathOutput"}="NA";


my @defaultpars=("pathInput","pathOutput");
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

print "Reading input and writing output...\n";

my $pathin=$parameters{"pathInput"};
my @s=split("\\.",$pathin);
my $ext=$s[-1];

my $input;

if($ext eq "gz"){
    open($input, "zcat $pathin |");
} else{
    open($input, $pathin);
}

open(my $output, ">".$parameters{"pathOutput"});
my $line=<$input>;

while($line){
    my $firstchar=substr $line, 0, 1;

    if($firstchar eq ">"){
	chomp $line;
	my @s=split(" ", $line);
	print $output $s[0]."\n";
    } else{
	## replace ambiguous nucleotides
	$line =~ tr/K|R|M|W|Y|k|r|m|w|y|/N/;
	
	print $output $line;
    }
    
    $line=<$input>;
}

close($output);
close($input);

print "Done.\n";

##############################################################

