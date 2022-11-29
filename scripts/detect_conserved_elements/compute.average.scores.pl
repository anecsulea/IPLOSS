use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

##############################################################

sub readScores{
    ## Wiggle (variableStep and fixedStep) is the only format defined by UCSC that uses a 1-based format, for historical reasons.
    
    my $pathin=$_[0];
    my $refphast=$_[1];
    
    my @s=split("\\.",$pathin);
    my $nbs=@s;
    my $ext=$s[$nbs-1];

    my $input;
    
    if($ext eq "gz"){
	open($input,"zcat $pathin |");
    }
    else{
	open($input,$pathin);
    }
    
    my $line=<$input>;
    my $currentpos="NA";
    my $step="NA";
    
    while($line){
	my $first=substr $line,0,1;
	
	if($first eq "f"){
	    my @s=split(" ",$line);
	    my $start=$s[2];
	    my @t=split("=",$start);
	    $currentpos=$t[1]+0; ## 1-based

	    my $stepinfo=$s[3];
	    my @u=split("=",$stepinfo);
	    $step=$u[1]+0;
	    
	    $line=<$input>;
	    next;
	}
	else{
	    chomp $line;
	    my $val=$line+0.0;

	    $refphast->{$currentpos}=$val;
	    
	    $currentpos+=$step; 
	    $line=<$input>;
	}
    }
    
    close($input);
}

##############################################################

sub computeScore{
    my $coords=$_[0];
    my $phast=$_[1];
    
    my $nbelements=@{$coords->{"start"}};

    for(my $i=0; $i<$nbelements; $i++){
	my $start=${$coords->{"start"}}[$i];
	my $end=${$coords->{"end"}}[$i];

	my $sumscore=0;
	my $nbpos=0;

	for(my $j=$start; $j<=$end; $j++){
	    
	    if(exists $phast->{$j}){
		$nbpos++;
		$sumscore+=$phast->{$j};
	    }
	}
	
	if($nbpos>0){
	    ${$coords->{"score"}}[$i]=($sumscore+0.0)/($nbpos+0.0);
	} 
	
	${$coords->{"coveredbases"}}[$i]=$nbpos+0.0;
	
    }
}

##############################################################

sub readElementCoordinates{
    my $pathin=$_[0];
    my $refblocks=$_[1];
    
    $refblocks->{"id"}=[];
    $refblocks->{"start"}=[];
    $refblocks->{"end"}=[];
    $refblocks->{"score"}=[];
    $refblocks->{"coveredbases"}=[];
    
    open(my $input, $pathin);
    my $line=<$input>;
    
    while($line){
	chomp $line;
	my @s=split("\t",$line);
	
	my $id=$s[3];
	my $start=$s[1]+1; ## 1-based now
	my $end=$s[2]+0;

	push(@{$refblocks->{"id"}}, $id);
	push(@{$refblocks->{"start"}}, $start);
	push(@{$refblocks->{"end"}}, $end);
	push(@{$refblocks->{"score"}}, 0);
	push(@{$refblocks->{"coveredbases"}}, 0);
	    
	$line=<$input>;
    }
    close($input);
}


##########################################################################

sub printHelp{

    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script computes phastCons or phyloP scores for most conserved elements. \n";
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

$parameters{"chr"}="NA";
$parameters{"pathCoords"}="NA";
$parameters{"pathScores"}="NA";
$parameters{"pathOutput"}="NA";


my %defaultvalues;
my @defaultpars=("chr", "pathCoords",  "pathScores",  "pathOutput");

my %numericpars;

foreach my $par (keys %parameters){
    $defaultvalues{$par}=$parameters{$par};
}

## check if help was asked 

foreach my $arg (@ARGV){
    if($arg eq "--help"){
	printHelp(\@defaultpars, \%defaultvalues);
	exit(0);
    }
}

## check new parameters

my $nbargs=@ARGV;

for(my $i=0; $i<$nbargs; $i++){
    my $arg=$ARGV[$i];
    $arg=substr $arg,2;
    my @s=split("=",$arg);
    my $parname=$s[0];
    my $parval=$s[1];
    
    if(exists $parameters{$parname}){
	$parameters{$parname}=$parval;
	
	if(exists $numericpars{$parname}){
	    $parameters{$parname}=$parval+0.0;
	}
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

print "Reading genomic coordinates...\n";
my %elements;

print "coordinate convention: 0-based (bed format)\n";

readElementCoordinates($parameters{"pathCoords"}, \%elements);

print "Done.\n";

##############################################################

print "Computing score and writing output...\n";

my $chr=$parameters{"chr"};

open(my $output,">".$parameters{"pathOutput"});
print $output "ID\tChr\tStart\tEnd\tScore\tCoveredLength\n";

my $path=$parameters{"pathScores"};

if(-e $path){
 
    my %phastCons;
    readScores($path,\%phastCons);
    
    computeScore(\%elements,\%phastCons);
    
    my $nbel=@{$elements{"start"}};
		
    for(my $i=0;$i<$nbel;$i++){
	    
	my $id=${$elements{"id"}}[$i];
	my $start=${$elements{"start"}}[$i];
	my $end=${$elements{"end"}}[$i];
	my $score=${$elements{"score"}}[$i];
	my $cov=${$elements{"coveredbases"}}[$i];
	
	print $output $id."\t".$chr."\t".$start."\t".$end."\t".$score."\t".$cov."\n";
	
    }
}

close($output);

print "Done.\n";


##############################################################
