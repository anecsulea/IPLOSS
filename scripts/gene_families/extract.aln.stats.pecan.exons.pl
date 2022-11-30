use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

#########################################################################

sub readAlignments{
    my $pathin=$_[0];
    my $minlen=$_[1]; ## minimum alignment length (excluding gaps, for each species)
    my $aln=$_[2];

    my @s=split("\\.",$pathin);
    my $ext=$s[-1];

    my $input;

    if($ext eq "gz"){
	open($input,"zcat $pathin |");
    }
    else{
	open($input, $pathin);
    }
    
    my $line=<$input>;

    my %fasta;
    
    while($line){
	my $b=substr $line,0,1;
	
	if($b eq ">"){
	    chomp $line;
	    my $id=substr $line,1;

	    my @s=split(" ",$id);
	    $id=$s[0];
	    
	    $fasta{$id}="";

	    $line=<$input>;
	    $b=substr $line,0,1;
	    
	    while($line && !($b eq ">")){
		chomp $line;
		$fasta{$id}.=$line;
		$line=<$input>;
		$b=substr $line,0,1;
	    }
	}
    }

    close($input);

    ## now format alignments like for MAF

    my $alllong=1;

    my $index=0;

    $aln->{$index}={};
    
    foreach my $sp (keys %fasta){
	my $sequence=$fasta{$sp};
	my $nbgaps = $sequence =~ tr/-//;
	my $totlen = length $fasta{$sp};
	my $ungappedlen=$totlen-$nbgaps;
	
	## these are stitched multiple alignments! always on the + strand, always starting at 0, covering the entire length of the sequence

	$aln->{$index}{$sp}={"start"=>0,  "strand"=>"+", "sequence"=>$sequence, "chrlen"=>$ungappedlen, "ungappedlen"=>$ungappedlen};
	
	if($ungappedlen<$minlen){
	    $alllong=0;
	    last;
	}
    }
    	
    if($alllong==0){
	delete $aln->{$index};
    }
   
}

################################################################################

sub extractAlignmentStats{
    my $aln=$_[0];
    my $alnstats=$_[1];
  
    my @allindexes=keys %{$aln};
    my $firstindex=$allindexes[0];
    my @splist=keys %{$aln->{$firstindex}};
    my $firstsp=$splist[0];

    my $nbungapped=0;
    my $nbidentical=0;
    
    foreach my $index (keys %{$aln}){
	## now go over the alignment base by base

	my $alnlength=length $aln->{$index}{$firstsp}{"sequence"};

	for(my $i=0; $i<$alnlength; $i++){
	    my %bases;

	    my $allungap=1;
	    
	    foreach my $sp (keys %{$aln->{$index}}){
		my $base=uc (substr $aln->{$index}{$sp}{"sequence"}, $i, 1);

		if($base eq "-"){
		    $allungap=0;
		}
		else{
		    $bases{$base}=1;
		}
	    }
	    
	    if($allungap==1){
		$nbungapped++;

		my $nbdiffbases=keys %bases;

		if($nbdiffbases==1){
		    $nbidentical++;
		}
	    }
	}
    }

    $alnstats->{"ungapped"}=$nbungapped;
    $alnstats->{"identical"}=$nbidentical;
}

################################################################################

sub readClusters{
    my $pathin=$_[0];
    my $sp1=$_[1];
    my $sp2=$_[2];
    my $clusters=$_[3];
    
    open(my $input, $pathin);
    my $line=<$input>;
    chomp $line;
    my @s=split("\t", $line);
    my %header;
    for(my $i=0; $i<@s; $i++){
	$header{$s[$i]}=$i;
    }

    $line=<$input>;
    
    my $index=0;

    while($line){
	chomp $line;
	my @s=split("\t", $line);
	
	my $g1=$s[$header{"Genes.".$sp1}];
	my $g2=$s[$header{"Genes.".$sp2}];
	
	if($g1 ne "NA" && $g2 ne "NA"){
	    my @genes1=split(";", $g1);
	    my @genes2=split(";", $g2);
	    
	    $clusters->{$index}={$sp1=>[], $sp2=>[]};

	    push(@{$clusters->{$index}{$sp1}}, @genes1);
	    push(@{$clusters->{$index}{$sp2}}, @genes2);
		
	    $index++;
	}
	
	$line=<$input>;
    }
        
    close($input);
}
################################################################################

sub printHelp{
    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script extracts exon alignment stats from TBA alignments. \n";
    print "\n";
    print "Options:\n";
    
    foreach my $par (@{$parnames}){
	print "--".$par."  [  default value: ".$parvalues->{$par}."  ]\n";
    }
    print "\n";
}

################################################################################
################################################################################

my %parameters;
$parameters{"species1"}="NA";
$parameters{"species2"}="NA";
$parameters{"pathClusters"}="NA";
$parameters{"dirPecan"}="NA";
$parameters{"minAlignmentLength"}="NA";
$parameters{"pathOutput"}="NA";

my %defaultvalues;
my @defaultpars=("species1", "species2", "pathClusters", "dirPecan", "minAlignmentLength", "pathOutput");

my @numericpars=();


my %numericpars;

foreach my $par (@numericpars){
    $numericpars{$par}=1;
}

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

print "Running program with the following parameters:\n";

foreach my $par (@defaultpars){
    print "--".$par."=".$parameters{$par}."\n";
}

print "\n";

#####################################################################
#####################################################################

my $sp1=$parameters{"species1"};
my $sp2=$parameters{"species2"};
my @species=($sp1, $sp2);

#####################################################################

print "Reading clusters...\n";

my %clusters;
readClusters($parameters{"pathClusters"}, $sp1, $sp2, \%clusters);

my $nbclust=keys %clusters;

print "Found ".$nbclust." non-empty clusters.\n";
print "Done.\n";

#####################################################################

print "Reading alignments and writing output...\n";

my $minlength=$parameters{"minAlignmentLength"}+0;

print "minimum length: ".$minlength."\n";

open(my $output, ">".$parameters{"pathOutput"});
print $output "ID.".$sp1."\tID.".$sp2."\tLengthUngapped\tLengthIdentical\n";

foreach my $idclust (keys %clusters){
    foreach my $gene1 (@{$clusters{$idclust}{$sp1}}){
	foreach my $gene2 (@{$clusters{$idclust}{$sp2}}){
	    my $pathPecan=$parameters{"dirPecan"}."/".$gene1.".".$gene2.".mfa.gz";
	    
	    if(-e $pathPecan){
		my %aln;
		readAlignments($pathPecan, $minlength, \%aln);
		
		my %alnstats;
		extractAlignmentStats(\%aln, \%alnstats);
		
		print $output $gene1."\t".$gene2."\t".$alnstats{"ungapped"}."\t".$alnstats{"identical"}."\n";
		
	    } else{
		print "Weird! cannot find ".$pathPecan." for ".$gene1." and ".$gene2."\n";
	    }
	}
    }
}

close($output);

print "Done.\n";

#####################################################################
#####################################################################
