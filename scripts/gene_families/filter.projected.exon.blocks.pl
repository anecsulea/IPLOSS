use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

##############################################################
##############################################################

sub readExonBlocks{
    my $pathin=$_[0];
    my $genes=$_[1];

    open(my $input, $pathin);

    my $line=<$input>;

    while($line){
	chomp $line;
	my @s=split("\t", $line);
	my $geneid=$s[0];

	my $chr=$s[2];
	my $start=$s[3]+0; ## 1-based
	my $end=$s[4]+0;
	my $strand=$s[5];

	my $exonid=$chr.",".$start.",".$end.",".$strand;

	if(exists $genes->{$geneid}){
	    $genes->{$geneid}{$exonid}=1;
	}
	else{
	    $genes->{$geneid}={$exonid=>1};
	}

	$line=<$input>;
    }

    close($input);
}

##############################################################

sub readProjectedExons{
    my $pathin=$_[0];
    my $exons=$_[1];

    open(my $input, $pathin);

    my $line=<$input>;

    my %multiple;

    while($line){
	chomp $line;

	my @s=split("\t", $line);
	my $chr=$s[0];
	my $start=$s[1]+1; ## 1-based coordinates, included
	my $end=$s[2]+0; ## 1-based coordinate, included
	my $id=$s[3];
	my $strand=$s[5];

	if(exists $exons->{$id}){
	    push(@{$exons->{$id}{"chr"}}, $chr);
	    push(@{$exons->{$id}{"start"}}, $start);
	    push(@{$exons->{$id}{"end"}}, $end);
	    push(@{$exons->{$id}{"strand"}}, $strand);
	}
	else{
	    $exons->{$id}={"chr"=>[$chr], "start"=>[$start], "end"=>[$end], "strand"=>[$strand]};
	}

	$line=<$input>;
    }

    close($input);

}
##############################################################

sub makeBlocks{
    my $coords=$_[0];
    my $margin=$_[1];
    my $blocks=$_[2];

    $blocks->{"chr"}=[];
    $blocks->{"start"}=[];
    $blocks->{"end"}=[];
    $blocks->{"strand"}=[];

    my $nb=@{$coords->{"chr"}};
    my %hashpos;

    for(my $i=0; $i<$nb; $i++){
	my $chr=${$coords->{"chr"}}[$i];
	my $start=${$coords->{"start"}}[$i];
	my $end=${$coords->{"end"}}[$i];
	my $strand=${$coords->{"strand"}}[$i];

	if(exists $hashpos{$chr}){
	    if(exists $hashpos{$chr}{$strand}){
		if(exists $hashpos{$chr}{$strand}{$start}){
		    if($end > $hashpos{$chr}{$strand}{$start}){
			$hashpos{$chr}{$strand}{$start}=$end;
		    }
		} else{
		    $hashpos{$chr}{$strand}{$start}=$end;
		}
	    } else{
		$hashpos{$chr}{$strand}={$start=>$end};
	    }
	} else{
	    $hashpos{$chr}={$strand=>{$start=>$end}};
	}
    }

    foreach my $chr (keys %hashpos){
	foreach my $strand (keys %{$hashpos{$chr}}){
	    my @uniquestart = keys %{$hashpos{$chr}{$strand}};
	    my @sortedstart = sort {$a <=> $b} @uniquestart;

	    my $nbex=@sortedstart;

	    my $currentstart=$sortedstart[0];
	    my $currentend=$hashpos{$chr}{$strand}{$currentstart};

	    for(my $u=1; $u<$nbex; $u++){
		my $thisstart=$sortedstart[$u];
		my $thisend=$hashpos{$chr}{$strand}{$thisstart};

		## cluster blocks if they overlap

		if($thisstart>=$currentstart && $thisstart<=($currentend+$margin)){
		    ## we only change the end if it's larger than the current position
		    if($thisend>$currentend){
			$currentend=$thisend;
		    }
		}
		else{
		    push(@{$blocks->{"chr"}},$chr);
		    push(@{$blocks->{"strand"}},$strand);
		    push(@{$blocks->{"start"}},$currentstart);
		    push(@{$blocks->{"end"}},$currentend);

		    $currentstart=$thisstart;
		    $currentend=$thisend;
		}
	    }

	    ## don't forget the last block
	    push(@{$blocks->{"chr"}},$chr);
	    push(@{$blocks->{"strand"}},$strand);
	    push(@{$blocks->{"start"}},$currentstart);
	    push(@{$blocks->{"end"}},$currentend);
	}
    }
}

##############################################################

sub combineProjectedCoordinates{
    my $projections=$_[0];
    my $margin=$_[1];
    my $combined=$_[2];
    my $outputrejected=$_[3];

    my %multiple;

    foreach my $id (keys %{$projections}){
	my %blocks;
	makeBlocks($projections->{$id}, $margin, \%blocks);

	my $nb=@{$blocks{"chr"}};

	if($nb==1){
	    my $chr=${$blocks{"chr"}}[0];
	    my $start=${$blocks{"start"}}[0];
	    my $end=${$blocks{"end"}}[0];
	    my $strand=${$blocks{"strand"}}[0];

	    $combined->{$id}={"chr"=>$chr, "start"=>$start, "end"=>$end, "strand"=>$strand};
	} else{
	    $multiple{$id}=1;
	    print $outputrejected "MultipleHits\t".$id."\t".join(",", @{$blocks{"chr"}})."\t".join(",", @{$blocks{"start"}})."\t".join(",", @{$blocks{"end"}})."\t".join(",", @{$blocks{"strand"}})."\n";
	}
    }

    my $nbmulti=keys %multiple;

    print "Found ".$nbmulti." duplicated or split exons.\n";
}

##############################################################

sub filterProjectedExons{
    my $projexons=$_[0];
    my $minsizeratio=$_[1];
    my $maxsizeratio=$_[2];
    my $filteredexons=$_[3];
    my $outLog=$_[4];

     my $nbbadsize=0;

    foreach my $id (keys %{$projexons}){
	my @s=split(",", $id);
	my $chr=$s[0];
	my $start=$s[1];
	my $end=$s[2];
	my $strand=$s[3];

	my $newchr=$projexons->{$id}{"chr"};
	my $newstart=$projexons->{$id}{"start"};
	my $newend=$projexons->{$id}{"end"};
	my $newstrand=$projexons->{$id}{"strand"};


	## we initialize projections at the expected positions
	## if transcript boundaries, we don't check further

	my $oldsize=$end-$start+1.0;
	my $newsize=$newend-$newstart+1.0;

	my $sizeratio=$newsize/$oldsize;

	if($sizeratio>=$minsizeratio && $sizeratio<=$maxsizeratio){
	    $filteredexons->{$id}={"chr"=>$newchr, "start"=>$newstart, "end"=>$newend, "strand"=>$newstrand};
	}
	else{
	    $nbbadsize++;
	    print $outLog $id." bad_size ".$oldsize." ".$newsize." ".$sizeratio."\n";
	}

    }

    print "Rejected ".$nbbadsize." exons with bad size ratios.\n";

}

##############################################################

sub printHelp{
    my $parnames=$_[0];
    my $parvalues=$_[1];

    print "\n";
    print "This script filters projected annotations.\n";
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

$parameters{"pathExonBlocks"}="NA";
$parameters{"pathProjectedExons"}="NA";
$parameters{"collapse"}="NA";
$parameters{"minSizeRatio"}="NA";
$parameters{"maxSizeRatio"}="NA";
$parameters{"pathOutputFilteredExons"}="NA";
$parameters{"pathOutputRejectedExons"}="NA";
$parameters{"pathOutputLog"}="NA";

my @defaultpars=("pathExonBlocks", "pathProjectedExons", "collapse",  "minSizeRatio", "maxSizeRatio", "pathOutputFilteredExons", "pathOutputRejectedExons", "pathOutputLog");

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

open(my $outputrejected, ">".$parameters{"pathOutputRejectedExons"});

##############################################################

print "Reading exon blocks...\n";

my %genes;

readExonBlocks($parameters{"pathExonBlocks"},  \%genes);

print "Done.\n";

##############################################################

print "Reading projected exons...\n";

my %initialprojectedexons;

readProjectedExons($parameters{"pathProjectedExons"}, \%initialprojectedexons);

my $nbproj=keys %initialprojectedexons;

print "Found ".$nbproj." projected exons.\n";

print "Done.\n";

##############################################################

print "Combining liftOver coordinates...\n";

my %projectedexons;

my $collapse=$parameters{"collapse"}+0;

print "Joining liftOver hits if closer than ".$collapse."bp.\n";

combineProjectedCoordinates(\%initialprojectedexons, $collapse, \%projectedexons, $outputrejected);

print "Done.\n";

##############################################################

print "Filtering projected exons...\n";

my $minratio=$parameters{"minSizeRatio"}+0.0;
my $maxratio=$parameters{"maxSizeRatio"}+0.0;

print "Size ratio: ".$minratio." ".$maxratio."\n";

open(my $outLog, ">".$parameters{"pathOutputLog"});

my %filteredexons;

filterProjectedExons(\%projectedexons,  $minratio, $maxratio, \%filteredexons, $outLog);

my $nbkept=keys %filteredexons;

print "Kept ".$nbkept." exons.\n";

close($outLog);

print "Done.\n";

##############################################################

print "Writing output for filtered exons...\n";

open(my $output, ">".$parameters{"pathOutputFilteredExons"});

print $output "GeneID\tExonID\tChr\tStart\tEnd\tStrand\n";

foreach my $gene (keys %genes){
    foreach my $id (keys %{$genes{$gene}}){
	if(exists $filteredexons{$id}){
	    my $outchr=$filteredexons{$id}{"chr"};

	    my $newstrand="NA";

	    if($filteredexons{$id}{"strand"} eq "+"){
		$newstrand="1";
	    }
	    else{
		if($filteredexons{$id}{"strand"} eq "-"){
		    $newstrand="-1";
		} else{
		    if($filteredexons{$id}{"strand"} eq "-1" || $filteredexons{$id}{"strand"} eq "1"){
			$newstrand=$filteredexons{$id}{"strand"};
		    } else{
			print "Weird strand for ".$id." ".$filteredexons{$id}{"strand"} ."!!!\n";
			exit(1);
		    }
		}
	    }

	    print $output $gene."\t".$id."\t".$outchr."\t".$filteredexons{$id}{"start"}."\t".$filteredexons{$id}{"end"}."\t".$newstrand."\n";
	}
    }
}
close($output);

print "Done.\n";

##############################################################
print "Writing output for rejected exons...\n";

foreach my $gene (keys %genes){
    foreach my $id (keys %{$genes{$gene}}){
	if((exists $projectedexons{$id}) && (!exists $filteredexons{$id})){
	    print $outputrejected "SizeRatio\t".$gene."\t".$id."\t".$projectedexons{$id}{"chr"}."\t".$projectedexons{$id}{"start"}."\t".$projectedexons{$id}{"end"}."\t".$projectedexons{$id}{"strand"}."\n";
	}
    }
}
close($outputrejected);

print "Done.\n";

##############################################################
