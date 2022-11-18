use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

##############################################################
##############################################################

sub readChromosomesGFF{
    my $pathin=$_[0];
    my $chrsizes=$_[1];

    my @s=split("\\.", $pathin);
    my $ext=$s[-1];

    my $input;

    if($ext eq "gz"){
	open($input, "zcat $pathin | ");
    } else{
	open($input, $pathin);
    }
    
    my $line=<$input>;

    my $nbchr=0;
    
    while($line){
	chomp $line;
	my $prefix=substr $line,0,17;

	if($prefix eq "##sequence-region"){
	    my @s=split(" ",$line);
	    my $chr=$s[1];
	    my $size=$s[3]+0;

	    if(exists $chrsizes->{$size}){
		$chrsizes->{$size}{$chr}=1;
	    } else{
		$chrsizes->{$size}={$chr=>1};
	    }

	    $nbchr++;
	}
		
	$line=<$input>;
    }

    close $input;

    print "Saw ".$nbchr." sequence regions.\n";
}

##############################################################

sub splitGFF{
    my $pathin=$_[0];
    my $chrparts=$_[1];
    my $dirOutput=$_[2];
    my $prefixOutput=$_[3];

    my @s=split("\\.", $pathin);
    my $ext=$s[-1];

    my $input;

    if($ext eq "gz"){
	open($input, "zcat $pathin | ");
    } else{
	open($input, $pathin);
    }
    
    my $line=<$input>;

    my $nbchr=0;

    my $output;
    my $currentpart="NA";
    
    while($line){
	my $prefix=substr $line,0,1;

	if($prefix ne "#"){
	    my @s=split("\t", $line);
	    my $chr=$s[0];

	    if(!exists $chrparts->{$chr}){
		print "cannot find assigned part for ".$chr."\n";
		exit(1);
	    } else{
		my $part=$chrparts->{$chr};

		if($currentpart ne $part){
		    ## we close previous file

		    if($currentpart ne "NA"){
			close $output;
		    }

		    ## we open new file
		    open($output, ">>".$dirOutput."/".$prefixOutput.".part".${part}.".gff");
		    print $output $line;

		    $currentpart=$part;
		} else{
		    print $output $line;
		}
	    }
	    
	}
	
	$line=<$input>;
    }

    close $input;

    if($currentpart ne "NA"){
	close $output;
    }
}

##########################################################################

sub printHelp{

    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script splits GFF annotations into parts.\n";
    print "\n";
    
    print "Options:\n";
    
    foreach my $par (@{$parnames}){
	print "--".$par."  [  default value: ".$parvalues->{$par}."  ]\n";
    }
    print "\n";
}

##########################################################################
##########################################################################

## parameters 

my %parameters;
$parameters{"pathGFF"}="NA";
$parameters{"maxChrSize"}="NA";
$parameters{"dirOutput"}="NA";
$parameters{"prefixOutput"}="NA";

my @defaultpars=("pathGFF", "maxChrSize", "dirOutput", "prefixOutput");


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

print "Reading chromosome sizes...\n";

my %chrsizes;
readChromosomesGFF($parameters{"pathGFF"}, \%chrsizes);

print "Done.\n";

##############################################################

print "Splitting chromosomes by size...\n";

my $limit=$parameters{"maxChrSize"}+0;

print "limit ".$limit."\n";

my @sizes=keys %chrsizes;
my @sortedsizes=sort{$a<=>$b} @sizes;
my @decreasingsizes=reverse @sortedsizes;

my $currentpart=0;
my $currentsize=0;

my %parts;

my $nbparts=0;

foreach my $size (@decreasingsizes){
    foreach my $chr (keys %{$chrsizes{$size}}){
	my $thissize=$currentsize+$size;
	if($thissize<=$limit){
	    $currentsize=$thissize;
	    $parts{$chr}=$currentpart;
	} else{
	    $currentpart++;
	    $currentsize=$size;
	    $parts{$chr}=$currentpart;
	    
	    $nbparts++;
	}
    }
}

print "There will be ".$nbparts." parts.\n";

print "Done.\n";

##############################################################

print "Reading GFF and writing output...\n";

splitGFF($parameters{"pathGFF"}, \%parts, $parameters{"dirOutput"}, $parameters{"prefixOutput"});
 
print "Done.\n";

##############################################################
