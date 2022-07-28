#!/usr/bin/perl
use strict;
use List::Util qw( min max sum reduce);
use Data::Dumper;

our $AUTHOR =   '$Author: Thomas Nalpathamkalam <thomas.nalpathamkalam@sickkids.ca> $';

## Script to extract gene level information from annovar output
## Loads variant and exonic.varaints outputs and extract gene_symbols, assign gene-id, calculate splie 
## distance etc.
## called by main script 


my $var_gene       = $ARGV[0];
my $vex_gene       = $ARGV[1];
my $rlx_gene       = $ARGV[2];
my $gene_refFlat   = $ARGV[3];
my $aprs_gene      = $ARGV[4];
my $aprs_genex     = $ARGV[5];
my $var_chr        = $ARGV[6];
my $ent_gene_map   = $ARGV[7];
my $morbidmap_file = $ARGV[8];
my $ext_file       = $ARGV[9];
my $stat_file = $ARGV[10];
my $varsel_gene = $ARGV[11];
my $varxsel_gene = $ARGV[12];
my $wr_chr = $var_chr;

warn "Processing chrom: $wr_chr\n";

open( VGN, $ext_file ) or die "$! $ext_file";
open( VXN, $vex_gene ) or die "$! $vex_gene";
open( RLX, $rlx_gene ) or die "$! $rlx_gene";
open(STA, ">$stat_file") or die "$! $stat_file\n";


open( SGN, $varsel_gene ) or die "$! $varsel_gene ";
open( SXN, $varxsel_gene ) or die "$! $varxsel_gene";

 
my $fh_gen = IO::File->new( $ent_gene_map, q{<} )
  or usage("gene info file  , $ent_gene_map , is not readable ( $! )");
my $fh_omi = IO::File->new( $morbidmap_file, q{<} )
  or usage("Morbidmap file , $morbidmap_file , is not readable ( $! )");

#open(VXO, ">$out_gene") or die "$! $out_gene";

my $fh_aprsg;
my $fh_aprsx;

if ( $aprs_gene ne "0" ) {
	$fh_aprsg = IO::File->new( $aprs_gene, q{<} )
	  or die("gene appris gene file  , $aprs_gene , is not readable ( $! )");
	$fh_aprsx = IO::File->new( $aprs_genex, q{<} )
	  or
	  die("gene appris gene ex file  , $aprs_genex , is not readable ( $! )");
}

my $fh_rfl = IO::File->new( $gene_refFlat, q{<} )
  or die("gene refFlat file  , $gene_refFlat , is not readable ( $! )");
my %vxn     = ();
my %cd_lens = ();

my $stat_count = 0;
my %gene_exons = ();
my %tr_counts  = ();
my %refflats   = ();

while ( my $line = <$fh_rfl> ) {
	chomp($line);
	my $key     = ( split( "\t", $line ) )[1];
	my $rfl_chr = ( split( "\t", $line ) )[2];
	$rfl_chr =~ s/chr//;
	my $key_reff = $key . ":" . $rfl_chr;
	$refflats{$key_reff} = $line;

	# for transcript count
	my $key_gene = ( split( "\t", $line ) )[0];
	my $key_gene_rf = $rfl_chr . ":" . ( split( "\t", $line ) )[0];
	my $exon_starts = ( split( "\t", $line ) )[9];
	my $exon_ends   = ( split( "\t", $line ) )[10];
	if ( $key !~ m/^NR_/ ) {

		#print "$key -- not NR , $key_gene_rf\n";
		$gene_exons{$key_gene_rf}{'exstart'} =
		  ( $gene_exons{$key_gene_rf}{'exstart'} )
		  ? $gene_exons{$key_gene_rf}{'exstart'} . $exon_starts
		  : $exon_starts;
		$gene_exons{$key_gene_rf}{'exend'} =
		  ( $gene_exons{$key_gene_rf}{'exend'} )
		  ? $gene_exons{$key_gene_rf}{'exend'} . $exon_ends
		  : $exon_ends;
	}

	if ( exists $tr_counts{$key_gene_rf} ) {
		$tr_counts{$key_gene_rf} =
		  $tr_counts{$key_gene_rf} . ":" . ( split( "\t", $line ) )[1];
	}
	else {
		$tr_counts{$key_gene_rf} = ( split( "\t", $line ) )[1];
	}
	$stat_count++;

}

close $fh_rfl;

#load genex1
my %gx = ();

while ( my $line = <VXN> ) {
	chomp($line);
	my $st = "x";
	my @fields = split( "\t", $line );
	my $key =
	    $fields[3] . ":"
	  . $fields[4] . ":"
	  . $fields[5] . ":"
	  . $fields[6] . ":"
	  . $fields[7];
	if($wr_chr ne $fields[3]) { next; }else{  
	foreach my $each_gex ( split( ",", $fields[2] ) ) {
		
		
		
		my $this_gex_gene = ( split( ":", $each_gex ) )[0];
		#print "VXN: $this_gex_gene , $line\n";
		$gx{$key}->{$this_gex_gene}->{rfq} =
		  ( exists $gx{$key}->{$this_gex_gene}->{rfq} )
		  ? $gx{$key}->{$this_gex_gene}->{rfq} . "," . $each_gex
		  : $each_gex;
		$gx{$key}->{$this_gex_gene}->{eff} =
		  ( exists $gx{$key}->{$this_gex_gene}->{eff} )
		  ? $gx{$key}->{$this_gex_gene}->{eff} . "," . $fields[1]
		  : $fields[1];
		$gx{ $key . ":" . $st }->{$this_gex_gene}->{ $fields[1] } =
		  ( exists $gx{ $key . ":" . $st }->{$this_gex_gene}->{ $fields[1] } )
		  ? $gx{ $key . ":" . $st }->{$this_gex_gene}->{ $fields[1] } . ","
		  . $each_gex
		  : $each_gex;
	}
	}
}

my %gxr = ();

while ( my $line = <RLX> ) {
	chomp($line);
	my $st = "ux";
	my @fields = split( "\t", $line );
	my $key =
	    $fields[3] . ":"
	  . $fields[4] . ":"
	  . $fields[5] . ":"
	  . $fields[6] . ":"
	  . $fields[7];
	if($wr_chr ne $fields[3]) { next; }else{  
	foreach my $each_gex ( split( ",", $fields[2] ) ) {
		
		my $this_gex_gene = ( split( ":", $each_gex ) )[0];
		$gxr{$key}->{$this_gex_gene}->{rfq} =
		  ( exists $gxr{$key}->{$this_gex_gene}->{rfq} )
		  ? $gxr{$key}->{$this_gex_gene}->{rfq} . "," . $each_gex
		  : $each_gex;
		$gxr{$key}->{$this_gex_gene}->{eff} =
		  ( exists $gxr{$key}->{$this_gex_gene}->{eff} )
		  ? $gxr{$key}->{$this_gex_gene}->{eff} . "," . $fields[1]
		  : $fields[1];
		$gxr{ $key . ":" . $st }->{$this_gex_gene}->{ $fields[1] } =
		  ( exists $gxr{ $key . ":" . $st }->{$this_gex_gene}->{ $fields[1] } )
		  ? $gxr{ $key . ":" . $st }->{$this_gex_gene}->{ $fields[1] } . ","
		  . $each_gex
		  : $each_gex;
	}
	}
}

#intronic        NOC2L   1       889159  889159  A       C
#intronic        NOC2L   1       889638  889638  G       C
#intronic        NOC2L   1       889713  889713  C       A
#intronic        NOC2L   1       891021  891021  G       A
#intronic        NOC2L   1       891059  891059  C       T

my %apr_gn = ();
my %apr_gx = ();

while ( my $line = <$fh_aprsg> ) {
	chomp($line);
	my @fields = split( "\t", $line );
	if($fields[2] ne $wr_chr){ next ; }else{
	my $this_gn = get_genename( $fields[0], $fields[1] );
	my $key =
	    $fields[2] . ":"
	  . $fields[3] . ":"
	  . $fields[4] . ":"
	  . $fields[5] . ":"
	  . $fields[6];
	my $seq_overlap_type_apr = $fields[0];
	foreach my $g ( split( ",", $this_gn ) ) {
		if ( $seq_overlap_type_apr =~ /splicing/ ) {
			$apr_gn{$key}->{$g}->{eff} = $seq_overlap_type_apr;
		}
	}
	}
}

#close $fh_aprsg;

while ( my $line = <$fh_aprsx> ) {
	chomp($line);
	my @fields = split( "\t", $line );
	my $key =
	    $fields[3] . ":"
	  . $fields[4] . ":"
	  . $fields[5] . ":"
	  . $fields[6] . ":"
	  . $fields[7];
	if($wr_chr ne $fields[3]){ next; }else{  
	foreach my $each_gex ( split( ",", $fields[2] ) ) {
		if($wr_chr ne $fields[3]){ next; }else{
		my $g = ( split( ":", $each_gex ) )[0];
		$apr_gn{$key}->{$g}->{eff} =
		  ( exists $apr_gn{$key}->{$g}->{eff} )
		  ? $apr_gn{$key}->{$g}->{eff} . "," . $fields[1]
		  : $fields[1];
	}
	}
}
}

#refseq select ************************************************************
my %sel_gn = ();
my %sel_gx = ();

while ( my $line = <SGN> ) {
	chomp($line);
	my @fields = split( "\t", $line );
	my $this_gn = get_genename( $fields[0], $fields[1] );
	my $key =
	    $fields[2] . ":"
	  . $fields[3] . ":"
	  . $fields[4] . ":"
	  . $fields[5] . ":"
	  . $fields[6];
	
	if($wr_chr ne $fields[2]) { next; }else{  
	my $seq_overlap_type_apr = $fields[0];
	for my $g ( split( ",", $this_gn ) ) {
		$sel_gn{$key}->{$g}->{typeseq} = (exists $sel_gn{$key}->{$g}->{typeseq} ) ? $sel_gn{$key}->{$g}->{typeseq} . "," .   $seq_overlap_type_apr : $seq_overlap_type_apr;
	}
	}
}

#close $fh_aprsg;

while ( my $line = <SXN> ) {
	chomp($line);
	my @fields = split( "\t", $line );
	my $key =
	    $fields[3] . ":"
	  . $fields[4] . ":"
	  . $fields[5] . ":"
	  . $fields[6] . ":"
	  . $fields[7];
	if($wr_chr ne $fields[3]){ next; }else{  
	foreach my $each_gex ( split( ",", $fields[2] ) ) {
		
		my $g = ( split( ":", $each_gex ) )[0];
		my $tr = ( split( ":", $each_gex ) )[1];
		
		$sel_gn{$key}->{$g}->{eff} =
		  ( exists $sel_gn{$key}->{$g}->{eff} )
		  ? $sel_gn{$key}->{$g}->{eff} . "," . $fields[1]
		  : $fields[1];
		$sel_gn{$key}->{$g}->{trid} =
		  ( exists $sel_gn{$key}->{$g}->{trid} )
		  ? $sel_gn{$key}->{$g}->{trid} . "," . $tr
		  : $tr;
	}
	}
}

#***************************************************

#close $fh_aprsg;

my %entz_id           = ();
my %transcripts_entrz = ();
my %gene_desc         = ();
my %gene_mim          = ();
my %gene_types        = ();
while ( my $line = <$fh_gen> ) {
	chomp($line);
	my @fields    = split( "\t", $line );
	my $gen_tr_id = $fields[0];
	my $gsymbol   = $fields[1];
	my $ent_id    = $fields[2];

	#	print "Gene = $gsymbol\n";
	$entz_id{$gsymbol}             = $ent_id;
	$transcripts_entrz{$gen_tr_id} = $ent_id;
	my $gdisc = $fields[4];
	my $gtype = $fields[5];
	my $gmim  = $fields[6];
	my $cdl   = $fields[7];
	$gene_desc{$ent_id}  = $gdisc;
	$gene_types{$ent_id} = $gtype;
	$gene_mim{$ent_id}   = $gmim;

	#print "$gsymbol $ent_id $gmim\n";
	$cd_lens{$gen_tr_id} = $cdl;
	$stat_count++;
}
close $fh_gen;
warn "Done populating Gene-id hash\n";

my %db_giov = ();
while ( my $line = <$fh_omi> ) {
	chomp($line);
	my @fields = split( "\t", $line );
	my $egID = $fields[0];
	$db_giov{$egID}{'OMIM.ID'}          = $fields[3];
	$db_giov{$egID}{'OMIM.name'}        = $fields[4];
	$db_giov{$egID}{'gene.description'} = $fields[6];
	$db_giov{$egID}{'gene.type'}        = $fields[7];
	$stat_count++;
}
close $fh_omi;

my $rt_gene_h = load_genev_hash( $var_gene, $var_chr );

my %no_entrez_id_genes_list = ();
my $st                      = "x";
my %full_summary            = ();
my %stats_final_merg        = ();
my $multi_gene              = 0;

<VGN>;
while ( my $line = <VGN> ) {
	chomp($line);
	my @fields = split( "\t", $line );
	if ( $fields[0] eq $var_chr ) {
		my $key =
		    $fields[0] . ":"
		  . $fields[1] . ":"
		  . $fields[2] . ":"
		  . $fields[3] . ":"
		  . $fields[4];
		my $gene_symbol             = "N_A";
		my $seq_overlap_type        = "NA";
		my $refSeq_id               = "NA";
		my $gene_iter               = 0;
		my @dis_apr_effs            = ();
		my %var_outs                = ();
		my %var_seqOverlaps         = ();
		my %var_effects             = ();
		my %primary_var_seqOverlaps = ();
		my $primary_seqOverlap      = "NA";
		my $pri_refseq              = "NA";
		my $pri_refSeq_id_sorted    = "NA";
		my $cds_prop                = "NA";
		my $per_transcript_affected = "NA";
		my $eff_category            = "NA";
		my $eff_category_pri        = "NA";
		my $refSeq_id_sorted        = "NA";
		my $aa_flag                 = "NA";
		my ( $exon_dist_valueL, $exon_dist_valueR ) = ( "NA", "NA" );
		my $apr_eff      = "NA";
		my $typeseq_iter = 0;
		my %for_gene_iter = ();

		#my $aaflag = "NA";

		my %refseqid_for_effect = ();
		if ( exists $rt_gene_h->{$key} ) {
			my @gvar_kys = keys %{ $rt_gene_h->{$key} };
			foreach ( sort { $b cmp $a } @gvar_kys ) {
				$typeseq_iter++ unless ( $_ eq "refseqid" );

				if ( $_ eq "refseqid" ) { next; }
				$gene_symbol      = $rt_gene_h->{$key}->{$_};
				$seq_overlap_type = $_;
				$refSeq_id        = $rt_gene_h->{$key}->{'refseqid'};
				$gene_iter        = 0;
				$full_summary{$_}{'count'} =
				    ( exists $full_summary{$_}{'count'} )
				  ? ( $full_summary{$_}{'count'} + 1 )
				  : 1;
				$stats_final_merg{seq_overlap}{$seq_overlap_type} =
				    ( exists $stats_final_merg{seq_overlap}{$seq_overlap_type} )
				  ? ( $stats_final_merg{seq_overlap}{$seq_overlap_type} + 1 )
				  : 1;

				foreach my $each_gene ( split( ",", $gene_symbol ) ) {
					my $ucsc_anno = 0;
					$gene_iter++;
					$for_gene_iter{$each_gene} = (exists $for_gene_iter{$each_gene}) ? $for_gene_iter{$each_gene} + 1 : 1;
					my $exb_distance = getExb_distanceLow( $key, $each_gene );
					my $each_entz_id = "NA";
					$each_entz_id            = "NA";
					$cds_prop                = "NA";
					$eff_category            = "NA";
					$eff_category_pri        = "NA";
					$pri_refseq              = "NA";
					$refSeq_id_sorted        = "NA";
					$per_transcript_affected = "NA";
					$aa_flag                 = "NA";
					( $exon_dist_valueL, $exon_dist_valueR ) = ( "NA", "NA" );

					$apr_eff = "NA";
					my ( $exon_dist_valueL, $exon_dist_valueR ) =
					  ( "NA", "NA" );

					if ( exists $entz_id{$each_gene} ) {
						$each_entz_id = $entz_id{$each_gene};
					}
					else {
						$no_entrez_id_genes_list{$each_gene} = 1
						  unless ( $gene_symbol eq "N_A"
							|| $each_gene eq "NA" );
					}

					my $current_gene_desc =
					  (      $each_entz_id eq ""
						  || $db_giov{$each_entz_id}{'gene.description'} eq "" )
					  ? "NA"
					  : $db_giov{$each_entz_id}{'gene.description'};
					my $current_gene_type =
					  (      $each_entz_id eq ""
						  || $db_giov{$each_entz_id}{'gene.type'} eq "" )
					  ? "NA"
					  : $db_giov{$each_entz_id}{'gene.type'};
					my $current_gene_mim =
					  (      $each_entz_id eq ""
						  || $db_giov{$each_entz_id}{'OMIM.ID'} eq "" )
					  ? "NA"
					  : $db_giov{$each_entz_id}{'OMIM.ID'};
					my $current_mim_desc =
					  (      $each_entz_id eq ""
						  || $db_giov{$each_entz_id}{'OMIM.name'} eq "" )
					  ? "NA"
					  : $db_giov{$each_entz_id}{'OMIM.name'};

					my $mim_entry = "NA";
					if ( $current_gene_mim != "" || $current_mim_desc != "" ) {
						$mim_entry =
						  $current_gene_mim . "|" . $current_mim_desc;
					}

#			        my $current_mpo = ($each_entz_id eq "" || $gene_mpo{$each_entz_id} eq "") ? "NA" : $gene_mpo{$each_entz_id};
#    	    		my $current_hpo = ($each_entz_id eq "" || $gene_hpo{$each_entz_id} eq "" ) ? "NA" : $gene_hpo{$each_entz_id};
#        			my $current_cgd = ($each_entz_id eq "" || $gene_cgd{$each_entz_id} eq "") ? "NA\@NA" : $gene_cgd{$each_entz_id};
#        			my ($current_cgd_dis, $current_cgd_inh) = split('@', $current_cgd);
#					my $current_acmg = ($each_entz_id eq "" || $acmg{$each_entz_id} eq "" ) ? "NA" : $acmg{$each_entz_id};
					my $var_chr0 = $var_chr;
					$var_chr0 =~ s/chr//;
					my $each_gene_w_chr = $var_chr0 . ":" . $each_gene;
					my $all_isoforms    = $tr_counts{$each_gene_w_chr};
					my $no_of_all_isoforms =
					  scalar( split( ":", $all_isoforms ) );
					if ( $seq_overlap_type eq "exonic" ) {
						$refSeq_id =
						  ( exists $gx{$key}->{$each_gene}->{rfq} )
						  ? $gx{$key}->{$each_gene}->{rfq}
						  : "NA";
						my $eff_category0 =
						  ( exists $gx{$key}->{$each_gene}->{eff} )
						  ? $gx{$key}->{$each_gene}->{eff}
						  : "NA";
						my @dis_eff_cate = do {
							my %seen;
							grep { !$seen{$_}++ } split( ",", $eff_category0 );
						};
						$eff_category = join( ",", @dis_eff_cate );
						$eff_category_pri =
						  get_priortized_effect($eff_category);
						$pri_refseq =
						  exists( $gx{ $key . ":" . $st }->{$each_gene}
							  ->{$eff_category_pri} )
						  ? $gx{ $key . ":" . $st }->{$each_gene}
						  ->{$eff_category_pri}
						  : "NA";

#$apr_eff =   ( exists $apr_gn{$key}->{$each_gene}->{eff} )   ? $apr_gn{$key}->{$each_gene}->{eff}  : "NA";
						if ( $eff_category eq "NA" ) {

							$st        = "ux";
							$ucsc_anno = 1;
							$stats_final_merg{counts}{'refseq_UKNOWN_effect'} = (exists $stats_final_merg{counts}{'refseq_UKNOWN_effect'}) ? $stats_final_merg{counts}{'refseq_UKNOWN_effect'} + 1 : 1 ;
							$refSeq_id =
							  ( exists $gxr{$key}->{$each_gene}->{rfq} )
							  ? $gxr{$key}->{$each_gene}->{rfq}
							  : "unknown";
							my $eff_categoryU0 =
							  ( exists $gxr{$key}->{$each_gene}->{eff} )
							  ? $gxr{$key}->{$each_gene}->{eff}
							  : "UNKNOWN";

							my @dis_eff_cateU = do {
								my %seen;
								grep { !$seen{$_}++ }
								  split( ",", $eff_categoryU0 );
							};
							push( @dis_eff_cateU, 'unknown' );
							$eff_category = join( ",", @dis_eff_cateU );
							$eff_category_pri =
							  get_priortized_effect($eff_category);
							$pri_refseq =
							  exists( $gxr{ $key . ":" . $st }->{$each_gene}
								  ->{$eff_category_pri} )
							  ? $gxr{ $key . ":" . $st }->{$each_gene}
							  ->{$eff_category_pri}
							  : "NA";

						}
						$refSeq_id_sorted = sort_refseq_id($refSeq_id);
						$aa_flag          = get_aaflag($refSeq_id_sorted);
						$stats_final_merg{'exonic'}{$eff_category_pri} = (exists $stats_final_merg{'exonic'}{$eff_category_pri}) ? ($stats_final_merg{'exonic'}{$eff_category_pri} + 1) : 1; 

						#$pri_refSeq_id_sorted = sort_refseq_id($pri_refseq);

						#cds per
						$cds_prop =
						  sprintf( "%.3f", compute_cdsProp($refSeq_id) );
#						$cds_prop =  compute_cdsProp($pri_refseq);
						my $refSeq_id_tagged = get_tagged($refSeq_id , $sel_gn{$key}->{$each_gene}->{trid});
						$refSeq_id =   $refSeq_id_tagged;

						#transcript percentage
						my $u_tr_count = get_uniq_tr_counts($refSeq_id);
						$per_transcript_affected =
						  ( $u_tr_count / $no_of_all_isoforms ) * 100
						  unless ( $no_of_all_isoforms == 0 );

					}    #exonic

					#aprs
					$apr_eff =
					  ( exists $apr_gn{$key}->{$each_gene}->{eff} )
					  ? $apr_gn{$key}->{$each_gene}->{eff}
					  : "NA";
					@dis_apr_effs = do {
						my %seen;
						grep { !$seen{$_}++ } split( ",", $apr_eff );
					};
					my @dis_eff_cate = do {
						my %seen;
						grep { !$seen{$_}++ } split( ",", $eff_category );
					};
					my $eff_category_uniq = join( ",", @dis_eff_cate );

					if (
						check_precedence_ext(
							$seq_overlap_type,
							$var_seqOverlaps{$each_gene},
							$primary_var_seqOverlaps{$each_gene}
						) == 1
					  )
					{
						
						#my $refSeq_id_tagged = get_tagged($refSeq_id , $sel_gn{$key}->{$each_gene}->{trid});
						$var_outs{$each_gene}{'refseq_id'} =
						    ( $ucsc_anno == 1 && $fields[4] ne '0' )
						  ? ( "(Warning! wrong ORF) " . $refSeq_id )
						  : $refSeq_id;
						$var_outs{$each_gene}{'aa_flag'} = $aa_flag;
						$var_outs{$each_gene}{'effect'}  = $eff_category;
						$var_outs{$each_gene}{'effect_priority'} =
						  $eff_category_pri;

				   #$var_outs{$each_gene}{'effect_appris'} = $this_apprs_effect;
						$var_outs{$each_gene}{'gene_symbol'} = $each_gene;
						$var_outs{$each_gene}{'entrez_id'}   = $each_entz_id;
						$var_outs{$each_gene}{'gene_desc'} = $current_gene_desc;
						$var_outs{$each_gene}{'gene_type'} = $current_gene_type;

#$stats_final_merge{gene_types}{$current_gene_type} =  (exists $stats_final_merge{gene_types}{$current_gene_type})  ?  ($stats_final_merge{gene_types}{$current_gene_type} + 1) : 1;
						$var_outs{$each_gene}{'omim_id'} = $current_gene_mim;
						$var_outs{$each_gene}{'omim_phenotype'} =
						  $current_mim_desc;

			 #						$var_outs{$each_gene}{'MPO'} = $current_mpo;
			 #						$var_outs{$each_gene}{'HPO'} = $current_hpo;
			 #						$var_outs{$each_gene}{'CGD_disease'} = $current_cgd_dis;
			 #						$var_outs{$each_gene}{'CGD_inheritance'} = $current_cgd_inh;
			 #						$var_outs{$each_gene}{'ACMG_disease'} = $current_acmg;
						$var_outs{$each_gene}{'per_cds_affected'} = $cds_prop;
						$var_outs{$each_gene}{'per_transcripts_affected'} =
						  $per_transcript_affected;
						$var_outs{$each_gene}{'distance_spliceJunction'} =
						  $exb_distance;
						my @this_appris_effs = (exists $apr_gn{$key}->{$each_gene}->{eff}) ? split(",", $apr_gn{$key}->{$each_gene}->{eff}) : ("NA");
						my %seen =() ;
						my @this_appris_effs_U = grep { ! $seen{$_}++ } @this_appris_effs ;						
						$var_outs{$each_gene}{'effect_appris'} = join("," , @this_appris_effs_U);
						
						my @this_sel_effs = (exists $sel_gn{$key}->{$each_gene}->{eff}) ? split(",", $sel_gn{$key}->{$each_gene}->{eff}) : ("NA");
						%seen =() ;
						my @this_sel_effs_U = grep { ! $seen{$_}++ } @this_sel_effs ;						
						$var_outs{$each_gene}{'effect_select'} = join("," , @this_sel_effs_U);


						my @this_sel_typs = (exists $sel_gn{$key}->{$each_gene}->{typeseq}) ? split(",", $sel_gn{$key}->{$each_gene}->{typeseq}) : ("NA");
#						%seen =() ;
#						my @this_sel_effs_U = grep { ! $seen{$_}++ } @this_sel_effs ;						
						$var_outs{$each_gene}{'typeseq_select'} = join("," , @this_sel_typs);


#$sel_gn{$key}->{$g}->{typeseq} 

#($each_entz_id eq "" || $gene_exmat{$each_entz_id}{'mis_z'} eq "" ) ? "NA" : $gene_exmat{$each_entz_id}{'mis_z'};
#				foreach my $xm_hd (@xma__required){
#					my $this_xmt_key = "gnomAD_" . $xm_hd;
#
#					$var_outs{$each_gene}{$this_xmt_key} = ($each_entz_id eq "" || $gene_exmat{$each_entz_id}{$xm_hd} eq "" ) ? "NA" : $gene_exmat{$each_entz_id}{$xm_hd};
#				}

						$primary_seqOverlap = $seq_overlap_type;
						$primary_var_seqOverlaps{$each_gene} =
						  $primary_seqOverlap;

					}
					$var_seqOverlaps{$each_gene} =
					  ( exists $var_seqOverlaps{$each_gene} )
					  ? (
						$var_seqOverlaps{$each_gene} . ";" . $seq_overlap_type )
					  : $seq_overlap_type;

				}

			}    #@gvar_kys

		}    #exists $rt_gene_h->{$key}

		foreach my $this_gn ( keys %var_seqOverlaps ) {
			my $sorted_seqOverlaps =
			  join( ";", sort ( split( ";", $var_seqOverlaps{$this_gn} ) ) );

#my $this_ent_id = (exists $var_outs{$this_gn}{'entrez_id'}) ? $var_outs{$this_gn}{'entrez_id'} : die "Error, error with parsing entrez-id for $this_gn\n";
#$multi_typeseq_per_gene = ($sorted_seqOverlaps =~ m/;/) ? ($multi_typeseq_per_gene + 1) : $multi_typeseq_per_gene;
			my $primary_typseq_for_this_gene =
			  get_priortized_typeseq($sorted_seqOverlaps);
			my $sorted_primary_typseq_for_this_gene =
			  join( ";", sort ( split( ";", $primary_typseq_for_this_gene ) ) );
			  
			print $key, "\t", $var_outs{$this_gn}{'gene_symbol'}, "\t",
			  $var_outs{$this_gn}{'entrez_id'}, "\t",
			  $var_outs{$this_gn}{'refseq_id'}, "\t", $sorted_seqOverlaps, "\t",
			  $sorted_primary_typseq_for_this_gene, "\t",
			  $var_outs{$this_gn}{'effect'},          "\t",
			  $var_outs{$this_gn}{'effect_priority'}, "\t",
			  $var_outs{$this_gn}{'per_cds_affected'},
			  "\t", $var_outs{$this_gn}{'per_transcripts_affected'}, "\t",
			  $var_outs{$this_gn}{'distance_spliceJunction'}, "\t",
			  $var_outs{$this_gn}{'aa_flag'},                 "\t",
			  $var_outs{$this_gn}{'gene_desc'},               "\t",
			  $var_outs{$this_gn}{'gene_type'},               "\t",
			  $var_outs{$this_gn}{'omim_id'},                 "\t",
			  $var_outs{$this_gn}{'omim_phenotype'},          "\t",
			  $var_outs{$this_gn}{'effect_appris'} ,  "\t" , $var_outs{$this_gn}{'typeseq_select'} , "\t" , $var_outs{$this_gn}{'effect_select'} ,   "\n";
			  

		}    #$this_gn

#		if($gene_iter > 1){
#			$stats_final_merg{counts}{'multi-gene-overlap'} = (exists $stats_final_merg{counts}{'multi-gene-overlap'}) ? $stats_final_merg{counts}{'multi-gene-overlap'} + 1 : 1 ;
#			print "Multi =  $key\n";
#		}

		if(scalar(keys %for_gene_iter)  > 1){
			$multi_gene++;
		}
		if($typeseq_iter > 1){
			$stats_final_merg{counts}{'multi-sequenceOverlap'} = (exists $stats_final_merg{counts}{'multi-sequenceOverlap'}) ? $stats_final_merg{counts}{'multi-sequenceOverlap'} + 1 : 1 ;
		}

	}    #if chr
}    #end

close VGN;
$stats_final_merg{counts}{'multi-gene-overlap'} = $multi_gene ;



#print STA Dumper( \%stats_final_merg ), "\n";


foreach my $ck (keys %{$stats_final_merg{'counts'}}  ){
	print STA "General->" . $ck , " = " , $stats_final_merg{'counts'}{$ck} , "\n";
}

foreach my $ck (keys %{$stats_final_merg{'seq_overlap'}}  ){
	print STA "Typeseq->" . $ck , " = " , $stats_final_merg{'seq_overlap'}{$ck} , "\n";
}

foreach my $ck (keys %{$stats_final_merg{'exonic'}}  ){
	print STA "Effect->" . $ck , " = " , $stats_final_merg{'exonic'}{$ck} , "\n";
}

#close VXO;
warn "Done!\n";

sub getExb_distanceLow {
	my ( $sub_key, $sub_gene ) = @_;
	my $sub_dist   = "NA";
	my @key_fields = split( ":", $sub_key );
	my $sub_chr    = $key_fields[0];
	$sub_chr =~ s/chr//;
	my $gene_key = $sub_chr . ":" . $sub_gene;

	#print "Gene-key = $gene_key\n";
	my @dists = ();
	if (   exists $gene_exons{$gene_key}{'exstart'}
		&& exists $gene_exons{$gene_key}{'exend'} )
	{

		#my @exsts = split("," , $gene_exons{$gene_key}{'exstart'});
		#my @exeds = split("," , $gene_exons{$gene_key}{'exend'});
		foreach ( split( ",", $gene_exons{$gene_key}{'exstart'} ) ) {
			push( @dists, abs( $key_fields[1] - ( $_ + 1 ) ) );
		}
		foreach ( split( ",", $gene_exons{$gene_key}{'exend'} ) ) {
			push( @dists, abs( $key_fields[1] - $_ ) );
		}
		$sub_dist = min(@dists);
	}

	return $sub_dist;
}

sub get_genename {
	my ( $seq_overlap, $gn ) = @_;

	#print "$seq_overlap --- $gn\n";
	my $gene_name = "NA";
	my @parsed    = ();
	if ( $gn =~ m/(.+)\((.+)\)/ ) {
		my @spl = split( /\(.+?\)/, $gn );
		foreach my $g (@spl) {

			#$g =~ s/,//;
			#	print "g = $g\n";
			if ( $g =~ /(^,).+/ ) {
				$g =~ s/(^,)(.+)/$2/;
			}
			push( @parsed, $g );
		}
		$gene_name = join( ",", @parsed );
	}
	else {
		$gene_name = $gn;
	}

	if ( $seq_overlap eq "intergenic" ) {
		$gene_name = "NA";
	}
	return $gene_name;

}

sub get_priortized_effect {
	my ($ovrs) = @_;

	#	warn "Ins sub: $ovrs\n";
	my %prec = (
		'frameshift insertion'             => 1,
		'frameshift deletion'              => 2,
		'frameshift block substitution'    => 3,
		'stopgain'                         => 4,
		'stoploss'                         => 5,
		'nonframeshift insertion'          => 6,
		'nonframeshift deletion'           => 7,
		'nonframeshift block substitution' => 8,
		'nonsynonymous SNV'                => 9,
		'synonymous SNV'                   => 10,
		'UNKNOWN'                          => 11,
		'unknown'                          => 11,
		'unknown*'                         => 11,
	);

	my @ov       = split( ",", $ovrs );
	my $first    = pop @ov;
	my $pr_first = $prec{$first};
	my @out      = ($first);
	foreach (@ov) {
		my $this_pr = $prec{$_};
		if ( $this_pr == $pr_first ) {
			push( @out, $_ );
		}
		if ( $this_pr < $pr_first ) {
			@out      = ($_);
			$pr_first = $this_pr;
		}
	}
	return join( ";", @out );

}

sub sort_refseq_id {

#rid = PDE4DIP:NM_001002811:exon8:c.T2047C:p.L683L,PDE4DIP:NM_001002812:exon12:c.T1558C:p.L520L,PDE4DIP:NM_014644:exon12:c.T1558C:p.L520L,PDE4DIP:NM_001198832:exon15:c.T1756C:p.L586L,PDE4DIP:NM_001198834:exon12:c.T1558C:p.L520L,
	my ($rid)    = @_;
	my %tmp_rids = ();
	my @out_rid  = ();
	foreach my $this_rid ( split( ",", $rid ) ) {
		my $rid_tr = ( split( ":", $this_rid ) )[1];

		#$tmp_rids{$this_rid} = $cd_lens{$rid_tr};
	}

	foreach my $id ( sort { $tmp_rids{$b} <=> $tmp_rids{$a} } keys %tmp_rids ) {
		push( @out_rid, $id );
	}

	return join( ",", @out_rid );

}

sub get_uniq_tr_counts {
	my ($rfqid) = @_;
	my %trs = ();
	foreach my $r ( split( ",", $rfqid ) ) {
		my $this_trs = ( split( ":", $r ) )[1];
		$trs{$this_trs} = 1;
	}
	return scalar( keys %trs );
}

sub load_genev_hash {
	my ( $mm, $gchr ) = @_;
	warn "Loading $mm for $gchr\n";
	my %gva = ();

	#open(GV, $mm) or warn "$mm $!";
	my $fh_gv = IO::File->new( $mm, q{<} ) or warn "$mm  $!\n";
	while ( my $line = <$fh_gv> ) {
		chomp($line);
		my @fields = split( "\t", $line );
		my $key =
		    $fields[2] . ":"
		  . $fields[3] . ":"
		  . $fields[4] . ":"
		  . $fields[5] . ":"
		  . $fields[6];
		if ( $fields[2] eq $gchr ) {
			my $gn = get_genename( $fields[0], $fields[1] );

			#			warn "Gene = $gn\n";
			$gva{$key}->{ $fields[0] } = $gn;
			$gva{$key}->{'refseqid'} =
			    ( exists $gva{$key}->{'refseqid'} )
			  ? ( $gva{$key}->{'refseqid'} . ";" . $fields[1] )
			  : $fields[1];
		}

		#$gva{$key}->{'refseqid'} = $fields[1] ;

	}
	close $fh_gv;
	return \%gva;
}

sub check_precedence_ext {
	my ( $cur_ovr, $ovrs, $povr ) = @_;
	my %primary_ovrls = ();
	my %prec          = (
		'splicing'       => 1,
		'exonic'         => 1,
		'ncRNA_splicing' => 2,
		'ncRNA_exonic'   => 2,
		'UTR3'           => 5,
		'UTR5'           => 5,
		'intronic'       => 6,
		'ncRNA_intronic' => 7,
		'upstream'       => 8,
		'downstream'     => 8,
		'intergenic'     => 9,
	);
	my $ret = 1;
	my $this_prec = ( exists $prec{$cur_ovr} ) ? ( $prec{$cur_ovr} ) : 1100;
	foreach my $ov ( split( ";", $ovrs ) ) {
		my $pr = exists( $prec{$ov} ) ? ( $prec{$ov} ) : 1100;
		if ( $pr < $this_prec ) {
			$ret = 0;
		}
	}
	return $ret;
}

sub get_priortized_typeseq {
	my ($ovrs) = @_;

	#	warn "Ins sub: $ovrs\n";
	my %prec = (
		'exonic'         => 1,
		'splicing'       => 1,
		'ncRNA_splicing' => 2,
		'ncRNA_exonic'   => 2,
		'UTR3'           => 3,
		'UTR5'           => 3,
		'intronic'       => 4,
		'ncRNA_intronic' => 5,
		'upstream'       => 6,
		'downstream'     => 6,
		'intergenic'     => 7,
	);

	my @ov       = split( ";", $ovrs );
	my $first    = pop @ov;
	my $pr_first = $prec{$first};
	my @out      = ($first);
	foreach (@ov) {
		my $this_pr = $prec{$_};
		if ( $this_pr == $pr_first ) {
			push( @out, $_ );
		}
		if ( $this_pr < $pr_first ) {
			@out      = ($_);
			$pr_first = $this_pr;
		}
	}
	return join( ";", @out );

}

sub compute_cdsProp {
	my ($rid) = @_;
	my @cdsprps = ();
	foreach my $this_rid ( split( ",", $rid ) ) {
		my $rid_tr = ( split( ":", $this_rid ) )[1];
		my $this_tr_len = ( exists $cd_lens{$rid_tr} ) ? $cd_lens{$rid_tr} : 0;
		if ( $this_tr_len == 0 or $this_tr_len eq "" ) {
			push( @cdsprps, 0 );

			#	warn "$this_rid , " , 0 , "\n";
		}
		else {
			my $cdPos = 0;
			my $cdPos_rid = ( split( ":", $this_rid ) )[3];
			if ( $cdPos_rid =~ m/c\..*?(\d+).+/ ) {
				$cdPos = $1;
			}
			my $this_prop = $cdPos / $this_tr_len;
			my $cds_prop = ( $this_prop <= 1 ) ? $this_prop : 1;
			push( @cdsprps, ( 1 - $cds_prop ) );

		  #	warn "$cdPos , $this_tr_len ,  $this_rid , " , (1-$cds_prop) , "\n";
		}
	}
	return max(@cdsprps);
}

sub get_aaflag {
	my ($rfqid) = @_;
	my @trs            = split( ",", $rfqid );
	my $this_aa_change = "";
	my @am_changes     = ();
	my $flag           = 0;

	foreach (@trs) {
		my $aa = ( split( ":", $_ ) )[4];
		if ( $aa =~ m/p?(.)\d+(.)/ ) {
			if ( $this_aa_change ne "$1$2" ) {
				push( @am_changes, $aa );
			}
			$this_aa_change = "$1$2";

		}

	}
	if ( ( scalar @am_changes ) > 1 ) {
		$flag = 1;
	}
	return $flag;
}

sub get_tagged{
	my ($rq, $sel_tr0) = @_;
	my $sel_tr = (split('\.' , $sel_tr0))[0];
	my @rq_spl = split(',' , $rq);
	my $idx_spl = -1;
	my $this_idx = 0;
	foreach (@rq_spl){
		my $this_tr = (split(':' , $_))[1];
		if($this_tr eq $sel_tr){
			$idx_spl = $this_idx;
		}
		$this_idx++;
	}
	
	if($idx_spl == -1){
		warn "Cannot match trascript from Select : $rq vs $sel_tr\n";
	}else{
		my $tag_tr = $rq_spl[$idx_spl] . ':Select';
		$rq_spl[$idx_spl] = $tag_tr;
	}	
	
	return join("," , @rq_spl);
	
}
