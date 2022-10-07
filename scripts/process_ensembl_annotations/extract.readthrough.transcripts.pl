use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

##################################################################################################
##################################################################################################

sub readExonsTranscripts{
     my $path=$_[0];
     my $exonstx=$_[1];
     my $txexons=$_[2];

     open(my $input,$path);

     my $line=<$input>; ## header
     $line=<$input>;

     while($line){
	 chomp $line;
	 my @s=split("\t",$line);

	 my $exonid=$s[0];
	 my $txid=$s[1];

	 if(exists $exonstx->{$exonid}){
	     $exonstx->{$exonid}{$txid}=1;
	 }
	 else{
	     $exonstx->{$exonid}={$txid=>1};
	 }

	 if(exists $txexons->{$txid}){
	     $txexons->{$txid}{$exonid}=1;
	 }
	 else{
	     $txexons->{$txid}={$exonid=>1};
	 }


	 $line=<$input>;
     }

     close($input);
}

##################################################################################################

sub readTranscriptInfo{
    my $path=$_[0];
    my $txinfo=$_[1];

    open(my $input,$path);

     my $line=<$input>; ## header
     $line=<$input>;

     while($line){
	 chomp $line;
	 my @s=split("\t",$line);

	 my $geneid=$s[0];
	 my $txid=$s[1];
	 my $bio=$s[2];

	 $txinfo->{$txid}={"GeneID"=>$geneid, "Biotype"=>$bio};

	 $line=<$input>;
     }

     close($input);
}

##################################################################################################

sub readExonCoords{
    my $pathExons=$_[0];
    my $exons=$_[1];

    open(my $input,$pathExons);

    my $line=<$input>; ## header
    $line=<$input>;

    while($line){
	chomp $line;
	my @s=split("\t",$line);

	my $id=$s[0];
	my $chr=$s[1];
	my $start=$s[2]+0;
	my $end=$s[3]+0;
	my $strand=$s[4]+0;

	$exons->{$id}={"chr"=>$chr,"start"=>$start,"end"=>$end,"strand"=>$strand};

	$line=<$input>;
    }

    close($input);
}

###################################################################################################

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

###################################################################################################

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
				$overlap->{$id1}{$id2}={"start"=>$M,"end"=>$m,"fractionoverlap1"=>$fr1,"fractionoverlap2"=>$fr2};
			    }
			    else{
				$overlap->{$id1}={$id2=>{"start"=>$M,"end"=>$m,"fractionoverlap1"=>$fr1,"fractionoverlap2"=>$fr2}};
			    }
			}
		    }

		    $j++;
		}
	    }
	}
    }
}



############################################################################################################

sub addToCluster{
    my $refconnected=$_[0];
    my $refclusters=$_[1];
    my $indexclusters=$_[2];
    my $key=$_[3];

    if(exists $refconnected->{$key}){

        my $indexcluster="NA";

        if(exists $indexclusters->{$key}){
            $indexcluster=$indexclusters->{$key}
        }

        if($indexcluster eq "NA"){
            my $nbclusters=keys %{$refclusters};
            $indexcluster=$nbclusters+1;
            $refclusters->{$indexcluster}=[$key];
            $indexclusters->{$key}=$indexcluster;
        }

        foreach my $connection (keys %{$refconnected->{$key}}){
            if(!(exists $indexclusters->{$connection})){
                push(@{$refclusters->{$indexcluster}},$connection);
                $indexclusters->{$connection}=$indexcluster;
                addToCluster($refconnected,$refclusters,$indexclusters,$connection);
            }
        }

        delete $refconnected->{$key};
    }
}

#############################################################################################################

sub extractClusters{
    my $refconnected=$_[0];
    my $refclusters=$_[1];

    my %indexclusters;

    my $nbconnected=keys %{$refconnected};

    my $round=0;

    my %alreadyprinted;

    while($nbconnected>0){

        foreach my $key (keys %{$refconnected}){
            addToCluster($refconnected,$refclusters,\%indexclusters,$key);
            $nbconnected=keys %{$refconnected};
        }

        $round++;

        $nbconnected=keys %{$refconnected};
    }
}


####################################################################################################

sub readGeneInfo{
    my $pathInfo=$_[0];
    my $info=$_[1];

    open(my $input,$pathInfo);

    my $line=<$input>; ## header
    my %header;
    chomp $line;
    my @s=split("\t",$line);

    for(my $i=0;$i<@s;$i++){
	$header{$s[$i]}=$i;
    }

    $line=<$input>;

    while($line){
	chomp $line;
	my @s=split("\t",$line);

	my $id=$s[$header{"stable_id"}];
	my $thisbio=$s[$header{"biotype"}];
	my $descr=$s[$header{"description"}];

	$info->{$id}={"biotype"=>$thisbio,"description"=>$descr};

	$line=<$input>;
    }

    close($input);
}

###############################################################################################################

sub printHelp{

    my $parnames=$_[0];
    my $parvalues=$_[1];

    print "\n";
    print "This script identifies read-through transcripts. \n";
    print "\n";
    print "Options:\n";

    foreach my $par (@{$parnames}){
	print "--".$par."  [  default value: ".$parvalues->{$par}."  ]\n";
    }
    print "\n";
}

###############################################################################################################
###############################################################################################################

my %parameters;

$parameters{"pathExonCoords"}="NA";
$parameters{"pathExonsTranscripts"}="NA";
$parameters{"pathTranscriptInfo"}="NA";
$parameters{"pathGeneInfo"}="NA";
$parameters{"monoexonicBiotypes"}="NA";
$parameters{"pathOutput"}="NA";

my %defaultvalues;
my @defaultpars=("pathExonCoords", "pathExonsTranscripts", "pathTranscriptInfo", "pathGeneInfo", "monoexonicBiotypes","pathOutput");

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

#####################################################################
#####################################################################


print "Reading Ensembl gene info...\n";
my %geneinfo;

readGeneInfo($parameters{"pathGeneInfo"}, \%geneinfo);

print "Done.\n";

##############################################################

print "Extracting genes and exons with monoexonic biotypes...\n";

my @m=split(",", $parameters{"monoexonicBiotypes"});

print join(", ", @m)." are the monoexonic biotypes.\n";

my %monobio;

foreach my $mm (@m){
    $monobio{$mm}=1;
}

my %monogenes;

foreach my $gene (keys %geneinfo){
    if(exists $geneinfo{$gene}){

	my $bio=$geneinfo{$gene}{"biotype"};

	if(exists $monobio{$bio}){
	    $monogenes{$gene}=1;
	}
    }
}

my $nbmonogenes=keys %monogenes;

print "Found ".$nbmonogenes." monoexonic genes.\n";

print "Done.\n";

##############################################################

print "Reading exon coordinates...\n";

my %exoncoords;
readExonCoords($parameters{"pathExonCoords"}, \%exoncoords);
my $nbex=keys %exoncoords;

print "Found ".$nbex." exons.\n";

print "Done.\n";

print "Ordering exons...\n";

my %orderedexons;

orderExons(\%exoncoords, \%orderedexons);

print "Done.\n";

#####################################################################

print "Extracting overlaps...\n";

my %exonoverlap;

extractOverlap(\%orderedexons, \%orderedexons, 0, "sense", \%exonoverlap);

my $nbov=keys %exonoverlap;

print "Found ".$nbov." overlapping exons.\n";

print "Done.\n";

#####################################################################

print "Reading transcript info...\n";

my %transcriptinfo;

readTranscriptInfo($parameters{"pathTranscriptInfo"}, \%transcriptinfo);

my $nbtx=keys %transcriptinfo;

print "Found ".$nbtx." transcripts.\n";

print "Done.\n";

#####################################################################

print "Reading transcript-exon correspondence...\n";

my %exonstx;
my %txexons;

readExonsTranscripts($parameters{"pathExonsTranscripts"}, \%exonstx, \%txexons);

print "Done.\n";

#####################################################################

print "Extracting overlapping transcripts...\n";

my %txoverlap;

foreach my $tx1 (keys %txexons){
    my $gene1=$transcriptinfo{$tx1}{"GeneID"};

    foreach my $exon1 (keys %{$txexons{$tx1}}){
	if(exists $exonoverlap{$exon1}){
	    foreach my $exon2 (keys %{$exonoverlap{$exon1}}){
		if($exon1 ne $exon2){
		    foreach my $tx2 (keys %{$exonstx{$exon2}}){

			my $gene2=$transcriptinfo{$tx2}{"GeneID"};

			if($tx1 ne $tx2 && $gene1 ne $gene2){
			    if(exists $txoverlap{$tx1}){
				$txoverlap{$tx1}{$tx2}=1;
			    }
			    else{
				$txoverlap{$tx1}={$tx2=>1};
			    }
			}
		    }
		}
	    }
	}
    }
}

my $nbtxov=keys %txoverlap;

print "Found ".$nbtxov." overlapping transcripts.\n";

print "Done.\n";

#####################################################################

print "Extracting read-through transcripts...\n";

my %ovgenes;
my %ovtx;

foreach my $tx1 (keys %txoverlap){
    my $gene1=$transcriptinfo{$tx1}{"GeneID"};

    if(!exists($monogenes{$gene1})){
	foreach my $tx2 (keys %{$txoverlap{$tx1}}){
	    my $nbexons2=keys %{$txexons{$tx2}};
	    my $gene2=$transcriptinfo{$tx2}{"GeneID"};

	    if($gene1 ne $gene2 && (!exists($monogenes{$gene2}))){
		if(exists $ovgenes{$tx1}){
		    $ovgenes{$tx1}{$gene2}=1;
		}
		else{
		    $ovgenes{$tx1}={$gene2=>1};
		}

		if(exists $ovtx{$tx1}){
		    $ovtx{$tx1}{$tx2}=1;
		}
		else{
		    $ovtx{$tx1}={$tx2=>1};
		}
	    }
	}
    }
}

print "Done.\n";

#####################################################################

print "Filtering overlaps and writing output...\n";

open(my $output, ">".$parameters{"pathOutput"});

print $output "TranscriptID\tGeneID\tOtherTranscripts\tOtherGenes\n";

foreach my $tx1 (keys %txoverlap){
    my $gene1=$transcriptinfo{$tx1}{"GeneID"};

    if(exists $ovgenes{$tx1} && exists $ovtx{$tx1}){
	my %realovgenes;
	my %realovtx;

	foreach my $tx2 (keys %{$ovtx{$tx1}}){
	    my $gene2=$transcriptinfo{$tx2}{"GeneID"};
	    my $nbthisgenes=0;

	    if(exists $ovgenes{$tx2}){
		$nbthisgenes=keys %{$ovgenes{$tx2}};

		if($nbthisgenes<2){
		    $realovgenes{$gene2}=1;
		    $realovtx{$tx2}=1;
		}
	    }
	}

	my $nbreal=keys %realovgenes;

	## this transcript overlaps with at least 2 genes which cannot be classified as read-through themselves (they do not overlap with two genes)
	if($nbreal>=2){
	    print $output $tx1."\t".$gene1."\t".join(",", keys %realovtx)."\t".join(",", keys %realovgenes)."\n";
	}
    }
}

close($output);

print "Done.\n";

#####################################################################
#####################################################################
