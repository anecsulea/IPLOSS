#!/usr/bin/perl

use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

##############################################################

sub readGTF{
    my $pathin=$_[0];
    my $acctype=$_[1];
    my $exoncoords=$_[2];
    my $exontx=$_[3];
    my $genetx=$_[4];
    my $txgene=$_[5];
    my $txex=$_[6];

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

	if($type eq $acctype){
	    
	    my $chr=$s[0];
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

sub extractIntrons{
    my $txexons=$_[0];
    my $exoncoords=$_[1];
    my $introncoords=$_[2];
    
    foreach my $tx (keys %{$txexons}){
	my %hashexons;
	
	foreach my $exid (keys %{$txexons->{$tx}}){
	    my $thisstart=$exoncoords->{$exid}{"start"};
	    my $thisend=$exoncoords->{$exid}{"end"};

	    $hashexons{$thisstart}=$thisend;
	}

	my @sortedstart=sort {$a<=>$b} (keys %hashexons);
	my @sortedend;

	foreach my $start (@sortedstart){
	    push(@sortedend, $hashexons{$start});
	}

	my $nbex=@sortedstart;

	if($nbex>=2){
	    $introncoords->{$tx}={};
	    
	    for(my $i=0; $i<($nbex-1); $i++){
		my $end1=$sortedend[$i];
		my $start1=$sortedstart[$i+1];

		if($end1>$start1){
		    print "Weird! exons are not ordered for ".$tx."\n";
		    exit(1);
		}

		my $idint=($end1+1)."-".($start1-1);

		$introncoords->{$tx}{$idint}=1;
	    }
	}
    }
}

##############################################################

sub computeTranscriptInfo{
    my $txexons=$_[0];
    my $exoncoords=$_[1];
    my $txlen=$_[2];
    my $txcoords=$_[3];

    foreach my $tx (keys %{$txexons}){
	my $chr="NA";
	my $strand="NA";
	my $len=0;

	my $minstart="NA";
	my $maxend="NA";

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

	    if($minstart eq "NA"){
		$minstart=$thisstart;
	    } else{
		if($thisstart<$minstart){
		    $minstart=$thisstart;
		}
	    }

	     if($maxend eq "NA"){
		$maxend=$thisend;
	    } else{
		if($thisend>$maxend){
		    $maxend=$thisend;
		}
	    }
	}

	$txlen->{$tx}=$len;
	$txcoords->{$tx}={"start"=>$minstart, "end"=>$maxend};
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


#######################################################################

sub readFasta{
    my $path=$_[0];
    my $reffasta=$_[1];
    my $txgenes=$_[2];

    my @s=split("\\.",$path);
    my $ext=$s[-1];

    my $input;

    if($ext eq "gz"){
	open($input,"zcat $path |");
    }
    else{
	open($input, $path);
    }

    my $line=<$input>;

    while($line){
	my $b=substr $line,0,1;

	if($b eq ">"){
	    chomp $line;
	    my $info=substr $line,1;
	    
	    my @s=split(" ",$info);
	    my $id=$s[0];

	    my $geneid="NA";

	    foreach my $item (@s){
                my $prefix1=substr $item, 0, 5;

                if($prefix1 eq "gene:"){
                    my @t=split(":", $item);
                    $geneid=$t[1];
                    last;
                }
            }

	    if($geneid ne "NA"){
		$txgenes->{$id}=$geneid;
	    }

	    $reffasta->{$id}="";

	    $line=<$input>;
	    $b=substr $line,0,1;

	    while($line && !($b eq ">")){
		chomp $line;
		$reffasta->{$id}.=$line;
		$line=<$input>;
		$b=substr $line,0,1;
	    }
	}
    }

    close($input);
}

##########################################################################

sub writeSequence{
    my $sequence=$_[0];
    my $name=$_[1];
    my $output=$_[2];

    my $n=length $sequence;

    print $output ">".$name."\n";

    my $i=0;

    while($i<($n-60)){

        my $subseq=substr $sequence,$i,60;

        print $output $subseq ."\n";

        $i+=60;
    }

    if($i<$n){
        my $subseq=substr $sequence,$i;
        print $output $subseq ."\n";
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
$parameters{"pathNewCDSAnnotation"}="NA";
$parameters{"pathTranscriptAnnotation"}="NA";
$parameters{"minFractionOverlap"}="NA";
$parameters{"pathOutputOverlap"}="NA";
$parameters{"pathNewProteins"}="NA";
$parameters{"pathKnownProteins"}="NA";
$parameters{"pathOutputProteins"}="NA";

my @defaultpars=("pathNewCDSAnnotation",  "pathTranscriptAnnotation", "minFractionOverlap", "pathOutputOverlap", "pathNewProteins", "pathKnownProteins", "pathOutputProteins");


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

print "Reading CDS annotations...\n";

my %exoncoords1;
my %exontx1;
my %genetx1;
my %txgene1;
my %txexon1;

readGTF($parameters{"pathNewCDSAnnotation"}, "CDS", \%exoncoords1, \%exontx1, \%genetx1, \%txgene1, \%txexon1);

my $nbexons=keys %exoncoords1;
my $nbtx=keys %txgene1;
my $nbg=keys %genetx1;

my %txlen1;
my %txcoords1;
computeTranscriptInfo(\%txexon1, \%exoncoords1, \%txlen1, \%txcoords1);

my %introns1;
extractIntrons(\%txexon1,  \%exoncoords1, \%introns1);
  
print "Found ".$nbexons." exons, ".$nbtx." transcripts and ".$nbg." genes in first set of annotations (new CDS).\n";

##############################################################

print "Reading transcript annotations...\n";

my %exoncoords2;
my %exontx2;
my %genetx2;
my %txgene2;
my %txexon2;

readGTF($parameters{"pathTranscriptAnnotation"}, "exon", \%exoncoords2, \%exontx2, \%genetx2, \%txgene2, \%txexon2);

my $nbexons=keys %exoncoords2;
my $nbtx=keys %txgene2;
my $nbg=keys %genetx2;

my %txlen2;
my %txcoords2;
computeTranscriptInfo(\%txexon2, \%exoncoords2, \%txlen2, \%txcoords2);

my %introns2;
extractIntrons(\%txexon2,  \%exoncoords2, \%introns2);

print "Found ".$nbexons." exons, ".$nbtx." transcripts and ".$nbg." genes in second set of annotations (transcript annotations).\n";

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

open(my $output, ">".$parameters{"pathOutputOverlap"});

print $output "CDS\tGeneID1\tCDSLength\tNbIntronsCDS\tTranscriptID\tGeneID2\tTranscriptLength\tOverlapLength\tNbIntronsTranscript\tNbCommonIntrons\tAbsentIntronsOverlapCDS\n";

my %acceptablecds;

my $minfrov=$parameters{"minFractionOverlap"}+0.0;

print "We accept CDS if they overlap on at least ".$minfrov." of their length with annotated transcripts.\n";

foreach my $tx1 (keys %txgene1){
    my $gene1=$txgene1{$tx1};
    my $len1=$txlen1{$tx1};
    my $start1=$txcoords1{$tx1}{"start"};
    my $end1=$txcoords1{$tx1}{"end"};

    if(exists $txoverlap12{$tx1}){

	foreach my $tx2 (keys %{$txoverlap12{$tx1}}){
	    my $gene2=$txgene2{$tx2};
	    my $len2=$txlen2{$tx2};
	    
	    my $lenov=$txoverlap12{$tx1}{$tx2};

	    my $frov=($lenov+0.0)/($len1+0.0);

	    if($frov>=$minfrov){
		if(exists $acceptablecds{$tx1}){
		    $acceptablecds{$tx1}{$gene2}=1;
		} else{
		    $acceptablecds{$tx1}={$gene2=>1};
		}
	    }
	    
	    my $nbint1=keys %{$introns1{$tx1}};
	    my $nbint2=keys %{$introns2{$tx2}};
	    
	    my $nbcommonintrons=0;
	    
	    if(exists $introns1{$tx1}){
		if(exists $introns2{$tx2}){
		    foreach my $in1 (keys %{$introns1{$tx1}}){
			if(exists $introns2{$tx2}{$in1}){
			    $nbcommonintrons++;
			}
		    }
		}
		
		## check which introns from tx2 are not in tx1
		
		my %absentintrons;
		
		if(exists $introns2{$tx2}){
		    foreach my $in2 (keys %{$introns2{$tx2}}){
			if(!exists $introns1{$tx1}{$in2}){
			    $absentintrons{$in2}=1;
			}
		    }
		}
		
		my %absentoverlap;
		
		foreach my $in2 (keys %absentintrons){
		    my @s=split("-", $in2);
		    my $startint2=$s[0]+0;
		    my $endint2=$s[1]+0;
		    
		    my $M=max($start1, $startint2);
		    my $m=min($end1, $endint2);
		    
		    ## if intron from tx2 is not in CDS, then coordinates should not overlap with CDS
		    if($M<=$m){
			$absentoverlap{$in2}=1;
		    }
		}
		
		my $nbabsov=keys %absentoverlap;
		
		print $output $tx1."\t".$gene1."\t".$len1."\t".$nbint1."\t".$tx2."\t".$gene2."\t".$len2."\t".$lenov."\t".$nbint2."\t".$nbcommonintrons."\t".$nbabsov."\n";
		
	    }
	}
    } 
}

close($output);

print "Done.\n";

my $nbacc=keys %acceptablecds;

print "Found ".$nbacc." CDS that overlapped on at least ".$minfrov." of their length with an annotated transcript.\n";

##############################################################

print "Checking if the CDS could be attributed to a single gene.\n";

my %filtered;

foreach my $tx (keys %acceptablecds){
    my @genes=keys %{$acceptablecds{$tx}};

    if(@genes==1){
	my $gene=$genes[0];
	
	$filtered{$tx}=$gene;
    }
}

my $nbfiltered=keys %filtered;

print "Kept ".$nbfiltered." CDS.\n";

print "Done.\n";

##############################################################

print "Reading protein sequences...\n";

my %knownproteins;
my %knownprotgene;

readFasta($parameters{"pathKnownProteins"}, \%knownproteins, \%knownprotgene);

my %newproteins;
my %newprotgene;

readFasta($parameters{"pathNewProteins"}, \%newproteins, \%newprotgene);

print "Done.\n";

##############################################################

print "Writing combined protein sequences...\n";

open(my $output, ">".$parameters{"pathOutputProteins"});

my $censoredX=0;

foreach my $id (keys %knownproteins){
    if(exists $knownprotgene{$id}){
	my @s=split("\\.", $knownprotgene{$id});
	my $gene=$s[0];

	my $seq=$knownproteins{$id};
	my $name=$id." gene:".$gene;

	my $nbX = ($seq =~ tr/X//);

	if($nbX==0){
	    writeSequence($seq, $name, $output);
	} else{
	    $censoredX++;
	}
    } else{
	print "Weird! cannot find gene id for ".$id."\n";
	exit(1);
    }
}

print "Discarded ".$censoredX." sequences with Xs from known proteins.\n";

$censoredX=0;

foreach my $id (keys %newproteins){
    if(exists $filtered{$id}){
	my $gene=$filtered{$id};

	my $seq=$newproteins{$id};
	my $name=$id." gene:".$gene;
	
	my $nbX = ($seq =~ tr/X//);

	if($nbX==0){
	    writeSequence($seq, $name, $output);
	} else{
	    $censoredX++;
	}
    }
}

print "Discarded ".$censoredX." sequences with Xs from new proteins.\n";

close($output);

print "Done.\n";

##############################################################

