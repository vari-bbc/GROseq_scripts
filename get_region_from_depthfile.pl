#!/usr/bin/perl -w
use strict;
=head1
Author: Ian Beddows
Date: 3/9/18 

Short description:	
		
		
		./get.pl -region_file *.txt -depth deliverables/depth_of_coverage.filtered.txt
		
=cut
use Getopt::Long;

my $usage = <<EOF;
OPTIONS:
-region_file
-depth_file
-outfile
OPTIONAL:
-h|help = print usage
-b = uppercase arguments
EOF
my($help,$bold,$region_file,$depth_file,$outfile);
#=======================================================================
GetOptions(
	'region_file=s' => \$region_file,	# string
	'depth_file=s' => \$depth_file,	# string
	'outfile=s' => \$outfile,	# string
    #~ '=i' => \$,	# integer
	#~ '=f' => \$,	# floating point number
    #~ 'bold' => \$bold,	# flag
	'h|help' => \$help	# flag for help
);

if (defined($help)){die print "HELP:\n",$usage;}


# First load all depth of coverage data:
open(my $in,'<',$depth_file)||die;
open(my $out,'>',$outfile)||die;
my %samples=();
my %d=(); # depth
my $c=0;
while(<$in>){
	chomp;
	my @data = split('\t',$_);
	my $chr = shift @data;
	#~ next if $chr ne 'X';
	my $pos = shift @data;
	if(1..1){ # get header
		for(my $i=0; $i<@data; $i++){
			$samples{$i} = $data[$i];
			#~ print "Sample $i: $samples{$i}\n";
		}
		next;
	}
	
	#~ for(my $i=0; $i<@data; $i++){
		#~ my $sample = $samples{$i};
		#~ $d{$chr}{$pos}{$sample}=$data[$i];
		@{$d{$chr}{$pos}}=@data;
		
		#~ print "$sample\n\t$chr $pos\n\t\t$d{$chr}{$pos}{$sample}\n";
		#~ print "$sample\n\t@data\n";
	#~ }
	$c++;
	if($c % 1000000 == 0){
		print STDOUT "$c\t$chr\t$pos\n";
	}
	#~ if($chr>1){last;}
	#~ if($c>100){last;}
}
close($in);


foreach my $c (keys %d){
	#~ foreach my $p (keys %{$d{$c}}){
		#~ print "$c\t$p\n";
		#~ print join("\t",@{$d{$c}{$p}}),"\n";
	#~ }
	print "Found ",scalar keys %{$d{$c}}," positions on chromosome |$c|\n";
}

#~ exit;

# Then load the regions and foreach subset the depth data and save that as the output
open($in,'<',$region_file)||die;

while(<$in>){
	chomp;
	next if 1..1;
	my($ext_gene,$ensembl_gene_id,$exon_start,$exon_end,$distance,$strand,$chromosome_name)=split('\t',$_);
	print "$ext_gene\n";
	#~ $genes{$ext_gene}{'chromosome_name'}=$chromosome_name;
	#~ $genes{$ext_gene}{'exon_start'}=$exon_start;
	#~ $genes{$ext_gene}{'exon_end'}=$exon_end;
	#~ $genes{$ext_gene}{'strand'}=$strand;
	
	# go through all positions that are in the gene body +/- 50kb, keeping track of how far they are
	my $tss_distance;
	my $found_pos=0;
	if($strand==1){
		$tss_distance = -50000; # always starts -50kb from first exon
	}else{ # else if reverse strand then flip
		$tss_distance = $distance - 50000; # because distance includes +/- 50k
	}
	print "\tnow searching for postions from $exon_start to $exon_end\n";
	for(my $p=$exon_start;$p<=$exon_end;$p++){ # go from -50k to +50k of gene body 
		if(exists $d{$chromosome_name}{$p}){
			print $out "$ext_gene\t$p\t$tss_distance\t";
			print $out join("\t",@{$d{$chromosome_name}{$p}}),"\n";
			$found_pos++;
			#~ print STDOUT "Found position |$p| on |$chromosome_name| for $ext_gene\n";
		}else{
			#~ print STDOUT "No position |$p| on |$chromosome_name| for $ext_gene\n";
			#~ print $out "$ext_gene\t$p\t$tss_distance","\n";
		}
		
		if($strand==1){$tss_distance++;}else{$tss_distance--} # add or subtract distance as needed.
	}
	
	print "\tFound $found_pos positions for $ext_gene\n";
	
	#~ last; # test with ACTN2
	
	
}
close($in);

print "Done\n";
exit;



#=======================================================================
#( Subroutines                  )
# ------------------------------------ 
#  o
#   o   \_\_    _/_/
#    o      \__/
#           (oo)\_______
#           (__)\       )\/\
#               ||----w |
#               ||     ||











