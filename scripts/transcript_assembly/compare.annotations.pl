#!/usr/bin/perl

use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

##############################################################

sub readGTF{
    my $pathin=$_[0];
    my $chromo=$_[1];
    my $exoncoords=$_[2];
    my $exontx=$_[3];
    my $genetx=$_[4];
    my $txgene=$_[5];
    my $txex=$_[6];

    my $nbaccchr=keys %{$chromo};

    if($nbaccchr == 0){
	print "no filtering, we keep all chromosomes.\n";
    }

    open(my $input, $pathin);

    my $line=<$input>;
    my $prefix=substr $line, 0,1;

    while($prefix eq "#"){
	$line=<$input>;
	$prefix=substr $line, 0,1;
    }

    my $nbunstranded=0; ## we discard genes with "." strand

    my %rttx;

    while($line){
	chomp $line;
	my @s=split("\t", $line);

	my $type=$s[2];

	if($type eq "exon"){

	    my $chr=$s[0];

	    if($nbaccchr == 0 || exists $chromo->{$chr}){
		my $start=$s[3]+0; ## 1-based
		my $end=$s[4]+0;
		my $strand=$s[6];

		if($strand eq "+"){
		    $strand="1";
		} else{
		    if($strand eq "-"){
			$strand="-1";
		    } else{
			if($strand ne "."){
			    print "Unknown strand ".$strand." at line ".$line."\n";
			    exit(1);
			}
		    }
		}

		if($strand ne "."){
		    my $info=$s[8];
		    my @t=split(";", $info);
		    my $geneid=findInfo("gene_id", \@t);
		    my $txid=findInfo("transcript_id", \@t);
		    my $classcode=findInfo("class_code", \@t);
		    my $oldid=findInfo("oId", \@t);

		    if($txid eq "NA"){
			print "could not find transcript in ".$line."\n";
			exit(1);
		    }

		    if($geneid eq "NA"){
			print "could not find gene in ".$line."\n";
			exit(1);
		    }

		    my $exonid=$chr.",".$start.",".$end.",".$strand;

		    ## fill in exon coords

		    $exoncoords->{$exonid}={"chr"=>$chr,"start"=>$start, "end"=>$end, "strand"=>$strand};

		    ## exon - tx correspondence

		    if(exists $exontx->{$exonid}){
			$exontx->{$exonid}{$txid}=1;
		    }
		    else{
			$exontx->{$exonid}={$txid=>1};
		    }

		    if(exists $genetx->{$geneid}){
			$genetx->{$geneid}{$txid}=1;
		    }
		    else{
			$genetx->{$geneid}={$txid=>1};
		    }

		    if(exists $txgene->{$txid}){
			if($txgene->{$txid} ne $geneid){
			    print "Weird! ".$txid." is associated to more than one gene: ".$geneid. " ".$txgene->{$txid}."\n";
			}
		    } else{
			$txgene->{$txid}=$geneid;
		    }

		    if(exists $txex->{$txid}){
			$txex->{$txid}{$exonid}=1;
		    }
		    else{
			$txex->{$txid}={$exonid=>1};
		    }
		}
		else{
		    $nbunstranded++;
		}
	    }
	}

	$line=<$input>;
    }

    close($input);

    print "We discarded ".$nbunstranded." exons with undefined strand.\n";

}

##############################################################

sub findInfo{
    my $pattern=$_[0];
    my $array=$_[1];

    my $res="NA";

    my @grepres=grep(/${pattern}/,@{$array});

    my $nbg=@grepres;
    my $nbreal=0;

    if(@grepres==1){
	my @t=split("\"",$grepres[0]);
	$res=$t[1];
    } else{
	my $nbreal=0;

	foreach my $g (@grepres){
	    $g =~ s/^\s+//; ## remove whitespace
	    my @u=split(" ",$g);

	    if($u[0] eq $pattern){
		$nbreal++;
		my @t=split("\"",$g);
		$res=$t[1];
	    }
	}
    }

    if($nbreal>1){
	return "NA";
    }
    return $res;
}

##############################################################

sub computeTranscriptLength{
    my $txexons=$_[0];
    my $exoncoords=$_[1];
    my $txlen=$_[2];

    foreach my $tx (keys %{$txexons}){
	my $chr="NA";
	my $strand="NA";
	my $len=0;

	foreach my $exid (keys %{$txexons->{$tx}}){
	    my $thischr=$exoncoords->{$exid}{"chr"};
	    my $thisstrand=$exoncoords->{$exid}{"strand"};
	    my $thisstart=$exoncoords->{$exid}{"start"};
	    my $thisend=$exoncoords->{$exid}{"end"};

	    if($chr eq "NA"){
		$chr=$thischr;
	    } else{
		if($chr ne $thischr){
		    print "Weird! two different chromosomes for ".$tx."\n";
		    exit(1);
		}
	    }

	    if($strand eq "NA"){
		 $strand=$thisstrand;
	    } else{
		if($strand ne $thisstrand){
		    print "Weird! two different strands for ".$tx."\n";
		    exit(1);
		}
	    }

	    if($thisstart>$thisend){
		print "Weird coordinates for ".$exid."\n";
		exit(1);
	    }

	    $len+=($thisend-$thisstart+1);
	}

	$txlen->{$tx}=$len;
    }
}

##############################################################

sub orderExons{
    my $exons=$_[0];
    my $ordered=$_[1];

    my %hashpos;

    foreach my $exid (keys %{$exons}){
	my $chr=$exons->{$exid}{"chr"};
	my $strand=$exons->{$exid}{"strand"};
	my $start=$exons->{$exid}{"start"};
	my $end=$exons->{$exid}{"end"};

	if(exists $hashpos{$chr}){
	    if(exists $hashpos{$chr}{$start}){
		push(@{$hashpos{$chr}{$start}{"id"}},$exid);
		push(@{$hashpos{$chr}{$start}{"end"}},$end);
		push(@{$hashpos{$chr}{$start}{"strand"}},$strand);
	    }
	    else{
		$hashpos{$chr}{$start}={"id"=>[$exid],"end"=>[$end],"strand"=>[$strand]};
	    }
	}
	else{
	    $hashpos{$chr}={$start=>{"id"=>[$exid],"end"=>[$end],"strand"=>[$strand]}};
	}
    }

    foreach my $chr (keys %hashpos){
	$ordered->{$chr}={"start"=>[],"end"=>[],"strand"=>[],"id"=>[]};

	my @uniquestart=keys %{$hashpos{$chr}};
	my @sortedstart=sort {$a<=>$b} @uniquestart;

	foreach my $start (@sortedstart){
	    my $nb=@{$hashpos{$chr}{$start}{"end"}};

	    for(my $i=0; $i<$nb; $i++){
		my $end=${$hashpos{$chr}{$start}{"end"}}[$i];
		my $strand=${$hashpos{$chr}{$start}{"strand"}}[$i];
		my $id=${$hashpos{$chr}{$start}{"id"}}[$i];

		push(@{$ordered->{$chr}{"start"}}, $start);
		push(@{$ordered->{$chr}{"end"}}, $end);
		push(@{$ordered->{$chr}{"strand"}}, $strand);
		push(@{$ordered->{$chr}{"id"}}, $id);
	    }
	}
    }
}

##############################################################

sub extractOverlap{
    my $coords1=$_[0]; ## ordered coordinates
    my $coords2=$_[1]; ## ordered coordinates
    my $margin=$_[2];
    my $type=$_[3];
    my $overlap=$_[4];

    foreach my $chr (keys %{$coords1}){
	if(exists $coords2->{$chr}){
	    my $nbex1=@{$coords1->{$chr}{"start"}};
	    my $nbex2=@{$coords2->{$chr}{"start"}};

	    my $firstj=0;

	    for(my $i=0; $i<$nbex1; $i++){

		my $start1=${$coords1->{$chr}{"start"}}[$i]-$margin;
		my $end1=${$coords1->{$chr}{"end"}}[$i]+$margin;
		my $strand1=${$coords1->{$chr}{"strand"}}[$i];
		my $id1=${$coords1->{$chr}{"id"}}[$i];

		my $j=$firstj;

		while($j<$nbex2 && ${$coords2->{$chr}{"end"}}[$j]<$start1){
		    $j++;
		}

		$firstj=$j;

		while($j<$nbex2 && ${$coords2->{$chr}{"start"}}[$j]<=$end1){

		    my $start2=${$coords2->{$chr}{"start"}}[$j];
		    my $end2=${$coords2->{$chr}{"end"}}[$j];
		    my $strand2=${$coords2->{$chr}{"strand"}}[$j];
		    my $id2=${$coords2->{$chr}{"id"}}[$j];

		    if(($strand1 eq $strand2 && $type eq "sense") || ($strand1 ne $strand2 && $type eq "antisense")){

			my $M=max($start1,$start2);
			my $m=min($end1,$end2);

			if($M<=$m){
			    my $fr1=($m-$M+1)/($end1-$start1+1);
			    my $fr2=($m-$M+1)/($end2-$start2+1);

			    if(exists $overlap->{$id1}){
				$overlap->{$id1}{$id2}={"start"=>$M,"end"=>$m};
			    }
			    else{
				$overlap->{$id1}={$id2=>{"start"=>$M,"end"=>$m}};
			    }
			}
		    }

		    $j++;
		}
	    }
	}
    }
}

##########################################################################

sub analyzeTranscriptOverlap{
    my $exonoverlap=$_[0];
    my $exontx1=$_[1];
    my $exontx2=$_[2];
    my $txoverlap1=$_[3];
    my $txoverlap2=$_[4];

    foreach my $ex1 (keys %{$exonoverlap}){
	foreach my $ex2 (keys %{$exonoverlap->{$ex1}}){
	    my $startov=$exonoverlap->{$ex1}{$ex2}{"start"};
	    my $endov=$exonoverlap->{$ex1}{$ex2}{"end"};
	    my $lenov=$endov-$startov+1;

	    foreach my $tx1 (keys %{$exontx1->{$ex1}}){
		if(!exists $txoverlap1->{$tx1}){
		    $txoverlap1->{$tx1}={};
		}

		foreach my $tx2 (keys %{$exontx2->{$ex2}}){
		    if(!exists $txoverlap1->{$tx1}{$tx2}){
			$txoverlap1->{$tx1}{$tx2}=$lenov;
		    } else{
			$txoverlap1->{$tx1}{$tx2}+=$lenov;
		    }

		    if(!exists $txoverlap2->{$tx2}){
			$txoverlap2->{$tx2}={};
		    }

		    if(!exists $txoverlap2->{$tx2}{$tx1}){
			$txoverlap2->{$tx2}{$tx1}=$lenov;
		    } else{
			$txoverlap2->{$tx2}{$tx1}+=$lenov;
		    }
		}
	    }
	}
    }
}

##########################################################################

sub printHelp{

    my $parnames=$_[0];
    my $parvalues=$_[1];

    print "\n";
    print "This script compares two sets of annotations.\n";
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
$parameters{"pathAnnot1"}="NA";
$parameters{"pathAnnot2"}="NA";
$parameters{"pathOutput"}="NA";

my @defaultpars=("pathAnnot1", "pathAnnot2", "pathOutput");


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

my %chromo;

print "We keep all chromosomes.\n";

##############################################################

print "Reading first set of annotations...\n";

my %exoncoords1;
my %exontx1;
my %genetx1;
my %txgene1;
my %txexon1;

readGTF($parameters{"pathAnnot1"}, \%chromo, \%exoncoords1, \%exontx1, \%genetx1, \%txgene1, \%txexon1);

my $nbexons=keys %exoncoords1;
my $nbtx=keys %txgene1;
my $nbg=keys %genetx1;

print "Found ".$nbexons." exons, ".$nbtx." transcripts and ".$nbg." genes in first set of annotations.\n";

##############################################################

print "Reading second set of annotations...\n";

my %exoncoords2;
my %exontx2;
my %genetx2;
my %txgene2;
my %txexon2;

readGTF($parameters{"pathAnnot2"},  \%chromo, \%exoncoords2, \%exontx2, \%genetx2, \%txgene2, \%txexon2);

my $nbexons=keys %exoncoords2;
my $nbtx=keys %txgene2;
my $nbg=keys %genetx2;

print "Found ".$nbexons." exons, ".$nbtx." transcripts and ".$nbg." genes in second set of annotations.\n";

##############################################################

print "Computing transcript length...\n";

my %txlen1;
computeTranscriptLength(\%txexon1, \%exoncoords1, \%txlen1);

my %txlen2;
computeTranscriptLength(\%txexon2, \%exoncoords2, \%txlen2);

print "Done.\n";

##############################################################

print "Ordering exons...\n";

my %orderedexons1;
orderExons(\%exoncoords1, \%orderedexons1);

my %orderedexons2;
orderExons(\%exoncoords2, \%orderedexons2);

print "Done.\n";

##############################################################

print "Extracting overlap between exons...\n";

my %overlap21;
extractOverlap(\%orderedexons2, \%orderedexons1, 0, "sense", \%overlap21);

print "Done.\n";

##############################################################

print "Analyzing transcript overlap...\n";

my %txoverlap12;
my %txoverlap21;

analyzeTranscriptOverlap(\%overlap21, \%exontx2, \%exontx1, \%txoverlap21, \%txoverlap12);

print "Done.\n";

##############################################################

print "Writing output...\n";

open(my $output, ">".$parameters{"pathOutput"});

print $output "TranscriptID1\tGeneID1\tTranscriptID2\tGeneID2\tLengthTranscript1\tLengthTranscript2\tLengthOverlap\n";

foreach my $tx1 (keys %txgene1){
    my $gene1=$txgene1{$tx1};
    my $len1=$txlen1{$tx1};

    if(exists $txoverlap12{$tx1}){
	foreach my $tx2 (keys %{$txoverlap12{$tx1}}){
	    my $gene2=$txgene2{$tx2};
	    my $len2=$txlen2{$tx2};
	    
	    my $lenov=$txoverlap12{$tx1}{$tx2};
	    
	    print $output $tx1."\t".$gene1."\t".$tx2."\t".$gene2."\t".$len1."\t".$len2."\t".$lenov."\n";
	}
    } else{
	print $output $tx1."\t".$gene1."\tNA\tNA\t".$len1."\tNA\t0\n";
    }
}

close($output);

print "Done.\n";

##############################################################
