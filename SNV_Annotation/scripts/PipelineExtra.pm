#!/usr/bin/perl
package Parse::Annovar;
use strict;
use IO::File;
use List::Util qw( min max sum );
use File::Basename;

our $AUTHOR =   '$Author: Thomas Nalpathamkalam <thomas.nalpathamkalam@sickkids.ca> $';

## Module to format annotation outputs, called by the main script 

use constant VERSION    => "OCT-25-2019";

sub get_genename{
	my ($seq_overlap, $gn) = @_;
	#print "$seq_overlap --- $gn\n";
	my $gene_name = "NA";
	my @parsed = ();
	if($gn  =~ m/(.+)\((.+)\)/){
		my @spl =  split(/\(.+?\)/, $gn);
		foreach my $g (@spl){
			#$g =~ s/,//;
		#	print "g = $g\n";
			if($g =~ /(^,).+/){
				$g =~ s/(^,)(.+)/$2/;
			}
			push(@parsed , $g);
		}
		$gene_name = join(",", @parsed);
	}else{
		$gene_name = $gn
	}
	
	if($seq_overlap eq "intergenic"){
		$gene_name = "NA";
	}
	return $gene_name;
 	
}	

sub get_genename_old{
	my ($seq_overlap, $gn) = @_;
	my $gene_name = "NA";
	my @parsed = ();
	if($gn  =~ m/(.+)\((.+)\)/){
		my @spl =  split(/\(.+?\)/, $gn);
		foreach my $g (@spl){
			#$g =~ s/,//;
			push(@parsed , $g);
		}
		$gene_name = join(",", @parsed);
	}else{
		$gene_name = $gn
	}
	
	if($seq_overlap eq "intergenic"){
		$gene_name = "NA";
	}
	return $gene_name;
 	
}




sub get_uniq_tr_counts{
	my ($rfqid) = @_;
	my %trs = ();
	foreach my $r (split("," , $rfqid)){
		my $this_trs = (split(":", $r))[1];
		$trs{$this_trs} = 1;
	}
	return scalar (keys %trs);
}

sub get_aaflag{
	my ($rfqid) = @_;
	my @trs            = split( ",", $rfqid );
	my $this_aa_change = "";
	my @am_changes     = ();
	my $flag = 0;
	
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

sub get_annoVar_type{
	my ($rf, $at) = @_;
	my $vtype ="unknown";
	if($rf =~ /^[AGCT-]/ && $at =~ /^[AGCT-]/){
	if($rf eq "-"){
		$vtype = "ins";
	}elsif($at eq "-"){
		$vtype = "del";
	}elsif(length($rf) == length($at)){
		$vtype = (length($rf) == 1) ? "snp" : "mnp";
	}else{
		$vtype = "complex";
	}
	
	if(($rf eq $at) ){
		$vtype = ($rf eq "-") ?  "unknown" : "ref";
	}
	}
	return $vtype;
}


sub format_annotation{
	my ($ann_file , $smp_file, $type_p, $out_final,$out_filter, $add_cols, $ex_target_dist_out, $no_bgzip, $dng_smps) = @_;
	my $spl_thr = 100;
	my $rt_val = -1;
	#open(ANN, $ann_file) or die "$ann_file missing;  $!";
	my $fh_ann = IO::File->new($ann_file , q{<}) or die "$! $ann_file\n";
	#open(SMP, $smp_file) or die "$smp_file missing ;  $!";
	my $fh_smp = IO::File->new($smp_file , q{<}) or die "$! $smp_file\n";
	#open(OUT, ">$out_final") or die "error making $out_final ;  $!";
	my $fh_out = undef;
	if($no_bgzip != 0 ){
		$fh_out = IO::File->new($out_final , q{>}) or die "$! $out_final\n";
	}else{
	  open ($fh_out, "| bgzip -c > $out_final") or die "error writing $out_final $!\n";		
	}
		
	#open(OUTF, ">$out_filter") or die "error making $out_filter ;  $!";
	my $fh_outf = IO::File->new($out_filter , q{>}) or die "$! out_filter\n";
	
	if($type_p > 4 || $type_p < 0){
		die "Error: $type_p not allowed\n";
	}
	my @dng_cols = ("NULL_CONFIG","ML_NULL","PP_NULL","DNM_CONFIG","ML_DNM","PP_DNM","MQ","RD", "DNG_confidence");
#	#coverage file , modified Aug 17,2016
#	my %cov_info = ();
#	my @coverage_header = ();
	
	#added Jan 24,2019
	#For extra columns
	my @extra_cols = split("," , $add_cols);
	
	#target distance for exome 
	my %target_info = ();
	my $fh_ex_target;
	if(($type_p == 3 || $type_p == 3.7) && defined $ex_target_dist_out ){
		warn "exome file, loading target distances..\n";
		$fh_ex_target = IO::File->new($ex_target_dist_out , q{<}) or warn "$! $ex_target_dist_out\n";
		if(defined $fh_ex_target){
			while(my $tr_line = <$fh_ex_target>){
				chomp($tr_line);
				my ($inp , $res) = split('\|' , $tr_line ,2) ;
				#my $key = $inp;
				my @key_flds = split("\t", $inp);
				my $key = $key_flds[0] . "\t" . $key_flds[1] . "\t" . $key_flds[2] . "\t" . $key_flds[3] . "\t" . $key_flds[4];
				my $dist = (split('\|' , $res))[1];
				#warn "Dist = $key -- $dist\n";
				$target_info{$key} = $dist;
			}
			close $fh_ex_target;	
		}
	}
	


##	my @frqs = ("1000g_all" , "1000g_eur" , "1000g_amr" , "1000g_eas" , "1000g_afr" , "1000g_sas" , "NHLBI_all" , "NHLBI_aa" , "NHLBI_eu" , "ExAC_Freq" , "ExAC_AFR" , "ExAC_AMR",
#				 "ExAC_EAS" , "ExAC_FIN" , "ExAC_NFE" , "ExAC_OTH" , "ExAC_SAS" , "cg" , "cg_filtered" , "cgW597_AllFreq" , "cgW597_CalledFreq" ,
#				  "cg1KG436_AllFreq" , "cg1KG436_CalledFreq" );

#	my @frqs = ("1000g_all" , "1000g_eur" , "1000g_amr" , "1000g_eas" , "1000g_afr" , "1000g_sas" ,  "ExAC_Freq" , "ExAC_AFR" , "ExAC_AMR",
#				 "ExAC_EAS" , "ExAC_FIN" , "ExAC_NFE" , "ExAC_OTH" , "ExAC_SAS" ,
#				  "gnomAD_exome_ALL","gnomAD_exome_AFR","gnomAD_exome_AMR","gnomAD_exome_ASJ","gnomAD_exome_EAS","gnomAD_exome_FIN","gnomAD_exome_NFE","gnomAD_exome_OTH","gnomAD_exome_SAS",
#				  "gnomAD_genome_ALL","gnomAD_genome_AFR","gnomAD_genome_AMR","gnomAD_genome_ASJ","gnomAD_genome_EAS","gnomAD_genome_FIN","gnomAD_genome_NFE","gnomAD_genome_OTH"   );

#	my @frqs = ("A1000g_all" , "A1000g_eur" , "A1000g_amr" , "A1000g_eas" , "A1000g_afr" , "A1000g_sas" ,  "ExAC_Freq" , "ExAC_AFR" , "ExAC_AMR",
#				 "ExAC_EAS" , "ExAC_FIN" , "ExAC_NFE" , "ExAC_OTH" , "ExAC_SAS" , "gnomAD_exome_controls_AF_popmax", "gnomAD_genome_controls_AF_popmax" );

#	my @frqs = ("A1000g_all" , "A1000g_eur" , "A1000g_amr" , "A1000g_eas" , "A1000g_afr" , "A1000g_sas" ,  
#				  "gnomAD_genome211_controls_AF_popmax", "gnomAD_exome211_controls_AF_popmax" );
#
	my @frqs = ("freq_max");
				  
	my $type_col = "typeseq";
	my $spx_col = "spx_dpsi";
	my $dbsc1_col = "dbscSNV_ADA_SCORE";
	my $dbsc2_col = "dbscSNV_RF_SCORE";
	my $dist2splice_col = "distance_spliceJunction";
		
	my @samples = <$fh_smp>;
	chomp(@samples);
	
	my %var_type = ();
	my $header = <$fh_ann>;
	chomp($header);
	my @hds = split("\t" , $header);

	#store index of each header
	my %h_index = ();
	my $ind = 0;
	foreach (@hds){
		$h_index{$_} = $ind++; 
	}
	
	#my @new_hds = ("#CHROM", "start" , "end" , "ref_allele" , "alt_allele", "var_type" , "Original_VCFKEY", "MULTI_ALLELIC" , "FILTER"); # , "SNVSB" , "SNVHPOL" , "RU");
	my @new_hds = ("#CHROM", "start" , "end" , "ref_allele" , "alt_allele", "var_type");  # , "Original_VCFKEY", "MULTI_ALLELIC" , "FILTER");
	
	foreach (@extra_cols){
		push(@new_hds , $_);
	}
	
	push(@new_hds , "Original_VCFKEY");
	push(@new_hds , "MULTI_ALLELIC");
	push(@new_hds , "FILTER");
	

#56      MSSNG00131_003:NULL_CONFIG
#57      MSSNG00131_003:ML_NULL
#58      MSSNG00131_003:PP_NULL
#59      MSSNG00131_003:DNM_CONFIG
#60      MSSNG00131_003:ML_DNM
#61      MSSNG00131_003:PP_DNM
#62      MSSNG00131_003:MQ
#63      MSSNG00131_003:RD

	
	if($type_p == 2){
		push(@new_hds , "SNVSB");
		push(@new_hds , "SNVHPOL");
		push(@new_hds , "RU");
		foreach my $sm (@samples){
			my $key2go = $sm . ":" . "Zygosity";
			push(@new_hds , $key2go);
			$key2go = $sm . ":" . "Genotype";
			push(@new_hds , $key2go);
			$key2go = $sm . ":" . "DP";
			push(@new_hds , $key2go);
			$key2go = $sm . ":" . "DPF";
			push(@new_hds , $key2go);
			$key2go = $sm . ":" . "DPI";
			push(@new_hds , $key2go);
			$key2go = $sm . ":" . "AD_REF";
			push(@new_hds , $key2go);
			$key2go = $sm . ":" . "AD_ALT";
			push(@new_hds , $key2go);
			$key2go = $sm . ":" . "GQ";
			push(@new_hds , $key2go);
			$key2go = $sm . ":" . "GQX";
			push(@new_hds , $key2go);
			$key2go = $sm . ":" . "GT_PostNorm";
			push(@new_hds , $key2go);
			$key2go = $sm . ":" . "GT_PreNorm";
			push(@new_hds , $key2go);
			if(exists $dng_smps->{$sm}){
				foreach my $d_hd (@dng_cols){
					$key2go = $sm . ":" . $d_hd;
					push(@new_hds , $key2go);
				}
			}
			
		}
		
	}elsif($type_p == 1 || $type_p == 3){
		#"DP" , "FS" , "QD" , "MQ" , "MQRankSum" , "ReadPosRankSum" , "VQSLOD" , "culprit"
		push(@new_hds , "DP");
		push(@new_hds , "FS");
		push(@new_hds , "QD");
		push(@new_hds , "MQ");
		push(@new_hds , "GQ_MEAN") unless ($type_p != 1);
		push(@new_hds , "MQRankSum");
		push(@new_hds , "ReadPosRankSum");
		push(@new_hds , "VQSLOD");
		push(@new_hds , "culprit");
		foreach my $sm (@samples){
			my $key2go = $sm . ":" . "Zygosity";
			push(@new_hds , $key2go);
			$key2go = $sm . ":" . "Genotype";
			push(@new_hds , $key2go);
			$key2go = $sm . ":" . "DP";
			push(@new_hds , $key2go);
			$key2go = $sm . ":" . "AD_REF";
			push(@new_hds , $key2go);
			$key2go = $sm . ":" . "AD_ALT";
			push(@new_hds , $key2go);
			$key2go = $sm . ":" . "GQ";
			push(@new_hds , $key2go);
			$key2go = $sm . ":" . "PL";
			push(@new_hds , $key2go);
			$key2go = $sm . ":" . "GT_PostNorm";
			push(@new_hds , $key2go);
			$key2go = $sm . ":" . "GT_PreNorm";
			push(@new_hds , $key2go);
			$key2go = $sm . ":" . "PhaseGT";
			push(@new_hds , $key2go);
			$key2go = $sm . ":" . "PhaseID";
			push(@new_hds , $key2go);
			if(exists $dng_smps->{$sm}){
				foreach my $d_hd (@dng_cols){
					$key2go = $sm . ":" . $d_hd;
					push(@new_hds , $key2go);
				}
			}
			
			if($type_p == 3){
				#DP_TOTAL=430;DP_A=430;DP_C=0;DP_G=0;DP_T=0
				$key2go = $sm . ":" . "DP_TOTAL";
				push(@new_hds , $key2go);
				$key2go = $sm . ":" . "DP_A";
				push(@new_hds , $key2go);
				$key2go = $sm . ":" . "DP_C";
				push(@new_hds , $key2go);
				$key2go = $sm . ":" . "DP_G";
				push(@new_hds , $key2go);
				$key2go = $sm . ":" . "DP_T";
				push(@new_hds , $key2go);
				$key2go = $sm . ":" . "Dist_from_Target";
				push(@new_hds , $key2go);
#				$key2go = $sm . ":" . "N_count";
#				push(@new_hds , $key2go);
			}
		}		
	}elsif($type_p == 4){
		#somatic
		
		
	}elsif($type_p == 3.7) {
		#exome 3.7
		
#11      AC
#12      AF
#13      AN
#14      BaseQRankSum
#15      ClippingRankSum
#16      DB
#17      DP
#18      DS
#19      END
#20      ExcessHet
#21      FS
#22      HaplotypeScore
#23      InbreedingCoeff
#24      MLEAC
#25      MLEAF
#26      MQ
#27      MQRankSum
#28      QD
#29      RAW_MQ
#30      ReadPosRankSum
#31      SOR
#32      DP_TOTAL
#33      DP_A
#34      DP_C
#35      DP_G
#36      DP_T
#37      OVCF_KEY
#38      MULTI_ALLELIC
#39      BEEL:GTR
#40      BEEL:AD_REF
#41      BEEL:AD_ALT
#42      BEEL:DP
#43      BEEL:GQ
#44      BEEL:GT
#45      BEEL:MIN_DP
#46      BEEL:PGT
#47      BEEL:PID
#48      BEEL:PL
#49      BEEL:RGQ
#50      BEEL:SB
#51      BEEL:OGT
#52      BEEL:OZYG
#53      BEEL:Zygosity

		push(@new_hds , "BaseQRankSum");
		#push(@new_hds , "ClippingRankSum");
		push(@new_hds , "DP");
		#push(@new_hds , "ExcessHet");
		push(@new_hds , "FS");
		push(@new_hds , "MQ");
		push(@new_hds , "MQRankSum");
		push(@new_hds , "QD");
		push(@new_hds , "ReadPosRankSum");
		push(@new_hds , "SOR");
		foreach my $sm (@samples){
			my $key2go = $sm . ":" . "Zygosity";
			push(@new_hds , $key2go);
			$key2go = $sm . ":" . "Genotype";
			push(@new_hds , $key2go);
			$key2go = $sm . ":" . "DP";
			push(@new_hds , $key2go);
			$key2go = $sm . ":" . "AD_REF";
			push(@new_hds , $key2go);
			$key2go = $sm . ":" . "AD_ALT";
			push(@new_hds , $key2go);
			$key2go = $sm . ":" . "GQ";
			push(@new_hds , $key2go);
			$key2go = $sm . ":" . "PL";
			push(@new_hds , $key2go);
#			$key2go = $sm . ":" . "SB";
#			push(@new_hds , $key2go);
			$key2go = $sm . ":" . "GT_PostNorm";
			push(@new_hds , $key2go);
			$key2go = $sm . ":" . "GT_PreNorm";
			push(@new_hds , $key2go);
			$key2go = $sm . ":" . "PhaseGT";
			push(@new_hds , $key2go);
			$key2go = $sm . ":" . "PhaseID";
			push(@new_hds , $key2go);
			$key2go = $sm . ":" . "DP_TOTAL";
			push(@new_hds , $key2go);
			$key2go = $sm . ":" . "DP_A";
			push(@new_hds , $key2go);
			$key2go = $sm . ":" . "DP_C";
			push(@new_hds , $key2go);
			$key2go = $sm . ":" . "DP_G";
			push(@new_hds , $key2go);
			$key2go = $sm . ":" . "DP_T";
			push(@new_hds , $key2go);
			$key2go = $sm . ":" . "Dist_from_Target";
			push(@new_hds , $key2go);
			if(exists $dng_smps->{$sm}){
				foreach my $d_hd (@dng_cols){
					$key2go = $sm . ":" . $d_hd;
					push(@new_hds , $key2go);
				}
			}
			
		}		
		
		
	}
	
	
	
	
	my $typ_index = (exists $h_index{"typeseq"}) ? $h_index{"typeseq"} : die "column \'typeseq\' missing";
	
	
	my @ano_fields = @hds[$typ_index..$#hds];
	my @anno_fields2go = ();
	
	foreach (@ano_fields){
		my $val2go = $_;
		if($val2go ne "spx_dpsi_z"){
			push(@anno_fields2go , $val2go);
			if($val2go eq "typeseq"){
				$val2go = "typeseq";
			}
			if($val2go eq "typeseq_primary"){
				$val2go = "typeseq_priority";
			}
			push(@new_hds , $val2go);
		}
	}

	print $fh_out join("\t", @new_hds) , "\n";
	print $fh_outf join("\t", @new_hds) , "\tSUBSET_filter\n";	
	
	while(my $line = <$fh_ann>){
		chomp($line);
		my @fields = split("\t" , $line);
		my $dng_conf = 0;
		my $vcf_loc = $fields[0] . ":" . $fields[1] . ":" . $fields[3] . ":" . $fields[4];
		#my $vcf_vartype = (exists $stvc{$vcf_loc}) ? $stvc{$vcf_loc} : "NA";
		my %row_val = ();
		my $hi = 0;
		my @outline = ();
		foreach my $h (@hds){
			#my $this_val = $fields[$hi++];
#			$this_val =~ s/^\s+//; #remove leading spaces
#			$this_val =~ s/\s+$//; #remove trailing spaces
#			$this_val =~ s/^\s*(\S*(?:\s+\S+)*)\s*$/$1/;
#			if($this_val eq '.'){
#				$this_val = "NA";
#			}
			$row_val{$h} = $fields[$hi++];
		}
	
		my $m_alle = $row_val{'MULTI_ALLELIC'};
#		my $m_order = $row_val{'M_ORDER'};
#		if($m_alle == 0){
#			$m_order = 1;
#		}
		
		my $vcf_vartype = Parse::Annovar::get_annoVar_type($row_val{'ref_allele'},  $row_val{'alt_allele'});
		push(@outline , $row_val{'annovar_chr'});
		push(@outline , $row_val{'POS'});
		push(@outline , $row_val{'annovar_end'});
		push(@outline , $row_val{'REF'});
		push(@outline , $row_val{'ALT'});
		push(@outline , $vcf_vartype);
		
		#extra columns if any
		foreach my $excol (@extra_cols){
			push(@outline , $row_val{$excol});
		}
		push(@outline , $row_val{'OVCF_KEY'});
		push(@outline , $row_val{'MULTI_ALLELIC'});
		#push(@outline , $m_order);	
		push(@outline , $row_val{'FILTER'});
	
		my $ovcf_key = $row_val{'OVCF_KEY'};
#37      ref_allele
#38      alt_allele
		
		my $target_key = $row_val{'annovar_chr'} . "\t" . ($row_val{'annovar_start'} - 1) . "\t" . $row_val{'annovar_end'} . "\t"  . $row_val{'ref_allele'} . "\t" . $row_val{'alt_allele'};
		if($type_p == 2){
			push(@outline , exists($row_val{'SNVSB'})?$row_val{'SNVSB'}:"NA" );
			push(@outline , exists($row_val{'SNVHPOL'})?$row_val{'SNVHPOL'}:"NA" );
			push(@outline , exists($row_val{'RU'})?$row_val{'RU'}:"NA" );

			#sample wise
			foreach my $sm (@samples){
				my $key2go = $sm . ":" . "OZYG";
				push(@outline , $row_val{$key2go});
				$key2go = $sm . ":" . "GTR";
				push(@outline , $row_val{$key2go});
				$key2go = $sm . ":" . "DP";
				push(@outline , $row_val{$key2go});
				$key2go = $sm . ":" . "DPF";
				push(@outline , $row_val{$key2go});
				$key2go = $sm . ":" . "DPI";
				push(@outline , $row_val{$key2go});
				$key2go = $sm . ":" . "AD_REF";
				push(@outline , $row_val{$key2go});
				$key2go = $sm . ":" . "AD_ALT";
				push(@outline , $row_val{$key2go});
				$key2go = $sm . ":" . "GQ";
				push(@outline , $row_val{$key2go});
				$key2go = $sm . ":" . "GQX";
				push(@outline , $row_val{$key2go});
				$key2go = $sm . ":" . "GT";
				push(@outline , $row_val{$key2go});
				$key2go = $sm . ":" . "OGT";
				push(@outline , $row_val{$key2go});
				if(exists $dng_smps->{$sm}){
					foreach my $d_hd (@dng_cols){
						#DNG_confidence
						$key2go = $sm . ":" . $d_hd;
						push(@outline , $row_val{$key2go});
					}
					
					if($row_val{$sm . ":DNG_confidence"} eq "High"){
						$dng_conf = 1;
					}
					
				}

				
			}
						
		}elsif($type_p == 1 || $type_p == 3){
			push(@outline , exists($row_val{'DP'})?$row_val{'DP'}:"NA" );
			push(@outline , exists($row_val{'FS'})?$row_val{'FS'}:"NA" );
			push(@outline , exists($row_val{'QD'})?$row_val{'QD'}:"NA" );
			push(@outline , exists($row_val{'MQ'})?$row_val{'MQ'}:"NA" );	
			push(@outline , exists($row_val{'GQ_MEAN'})?$row_val{'GQ_MEAN'}:"NA" ) unless ($type_p != 1);								
			push(@outline , exists($row_val{'MQRankSum'})?$row_val{'MQRankSum'}:"NA" );			
			push(@outline , exists($row_val{'ReadPosRankSum'})?$row_val{'ReadPosRankSum'}:"NA" );
			push(@outline , exists($row_val{'VQSLOD'})?$row_val{'VQSLOD'}:"NA" );
			push(@outline , exists($row_val{'culprit'})?$row_val{'culprit'}:"NA" );
			
			#sample wise
			foreach my $sm (@samples){
				my $key2go = $sm . ":" . "OZYG";
				push(@outline , exists ($row_val{$key2go}) ? $row_val{$key2go} : "NA");
				$key2go = $sm . ":" . "GTR";
				#push(@outline , $row_val{$key2go});
				push(@outline , exists ($row_val{$key2go}) ? $row_val{$key2go} : "NA");
				$key2go = $sm . ":" . "DP";
				#push(@outline , $row_val{$key2go});
				push(@outline , exists ($row_val{$key2go}) ? $row_val{$key2go} : "NA");
				$key2go = $sm . ":" . "AD_REF";
				#push(@outline , $row_val{$key2go});
				push(@outline , exists ($row_val{$key2go}) ? $row_val{$key2go} : "NA");
				$key2go = $sm . ":" . "AD_ALT";
				#push(@outline , $row_val{$key2go});
				push(@outline , exists ($row_val{$key2go}) ? $row_val{$key2go} : "NA");
				$key2go = $sm . ":" . "GQ";
				#push(@outline , $row_val{$key2go});
				push(@outline , exists ($row_val{$key2go}) ? $row_val{$key2go} : "NA");
				$key2go = $sm . ":" . "PL";
				#push(@outline , $row_val{$key2go});
				push(@outline , exists ($row_val{$key2go}) ? $row_val{$key2go} : "NA");
				$key2go = $sm . ":" . "GT";
				#push(@outline , $row_val{$key2go});
				push(@outline , exists ($row_val{$key2go}) ? $row_val{$key2go} : "NA");
				$key2go = $sm . ":" . "OGT";
				#push(@outline , $row_val{$key2go});
				push(@outline , exists ($row_val{$key2go}) ? $row_val{$key2go} : "NA");
				$key2go = $sm . ":" . "PGT";
				push(@outline , exists($row_val{$key2go})?$row_val{$key2go}:"NA" );
				$key2go = $sm . ":" . "PID";
				push(@outline , exists($row_val{$key2go})?$row_val{$key2go}:"NA" );
				if(exists $dng_smps->{$sm}){
					foreach my $d_hd (@dng_cols){
						$key2go = $sm . ":" . $d_hd;
						push(@outline , exists($row_val{$key2go})?$row_val{$key2go}:"NA" );
					}
					if($row_val{$sm . ":DNG_confidence"} eq "High"){
						$dng_conf = 1;
					}
					
				}
				
				
				#new Aug 17, 2016
				#modified June 13,2017 ; expects modified vcf; #DP_TOTAL=430;DP_A=430;DP_C=0;DP_G=0;DP_T=0
				if($type_p == 3){
					$key2go = $sm . ":" . "DP_TOTAL";
					push(@outline , exists($row_val{$key2go})?$row_val{$key2go}:"NA" );
					$key2go = $sm . ":" . "DP_A";
					push(@outline , exists($row_val{$key2go})?$row_val{$key2go}:"NA" );
					$key2go = $sm . ":" . "DP_C";
					push(@outline , exists($row_val{$key2go})?$row_val{$key2go}:"NA" );
					$key2go = $sm . ":" . "DP_G";
					push(@outline , exists($row_val{$key2go})?$row_val{$key2go}:"NA" );
					$key2go =  $sm . ":" . "DP_T";
					push(@outline , exists($row_val{$key2go})?$row_val{$key2go}:"NA" );
					my $sample_target_distance = exists ($target_info{$target_key}) ? $target_info{$target_key} : "NA";
					push(@outline ,  $sample_target_distance);
					
#					my $this_read_count = (exists $cov_info{$cov_key}{'read_count'}) ? $cov_info{$cov_key}{'read_count'} : "NA"; 
#					push(@outline , $this_read_count);
#					my $this_read_countA = (exists $cov_info{$cov_key}{'A'}) ? $cov_info{$cov_key}{'A'} : "NA";
#					push(@outline , $this_read_countA);
#					my $this_read_countC = (exists $cov_info{$cov_key}{'C'}) ? $cov_info{$cov_key}{'C'} : "NA";
#					push(@outline , $this_read_countC);
#					my $this_read_countG = (exists $cov_info{$cov_key}{'G'}) ? $cov_info{$cov_key}{'G'} : "NA";
#					push(@outline , $this_read_countG);
#					my $this_read_countT = (exists $cov_info{$cov_key}{'T'}) ? $cov_info{$cov_key}{'T'} : "NA";
#					push(@outline , $this_read_countT);
#					my $this_read_countN = (exists $cov_info{$cov_key}{'N'}) ? $cov_info{$cov_key}{'N'} : "NA";
#					push(@outline , $this_read_countN);
				}
			}
		}elsif($type_p == 4){
			#somatic
		}elsif($type_p == 3.7){
			push(@outline , exists($row_val{'BaseQRankSum'})?$row_val{'BaseQRankSum'}:"NA" );
			#push(@outline , exists($row_val{'ClippingRankSum'})?$row_val{'ClippingRankSum'}:"NA" );
			push(@outline , exists($row_val{'DP'})?$row_val{'DP'}:"NA" );
			#push(@outline , exists($row_val{'ExcessHet'})?$row_val{'ExcessHet'}:"NA" );	
			push(@outline , exists($row_val{'FS'})?$row_val{'FS'}:"NA" ) ;								
			push(@outline , exists($row_val{'MQ'})?$row_val{'MQ'}:"NA" );
			push(@outline , exists($row_val{'MQRankSum'})?$row_val{'MQRankSum'}:"NA" );
			push(@outline , exists($row_val{'QD'})?$row_val{'QD'}:"NA" );
			push(@outline , exists($row_val{'ReadPosRankSum'})?$row_val{'ReadPosRankSum'}:"NA" );
			push(@outline , exists($row_val{'SOR'})?$row_val{'SOR'}:"NA" );
			foreach my $sm (@samples){
				my $key2go = $sm . ":" . "OZYG";
				push(@outline , $row_val{$key2go});
				$key2go = $sm . ":" . "GTR";
				push(@outline , $row_val{$key2go});
				$key2go = $sm . ":" . "DP";
				push(@outline , $row_val{$key2go});
				$key2go = $sm . ":" . "AD_REF";
				push(@outline , $row_val{$key2go});
				$key2go = $sm . ":" . "AD_ALT";
				push(@outline , $row_val{$key2go});
				$key2go = $sm . ":" . "GQ";
				push(@outline , $row_val{$key2go});
				$key2go = $sm . ":" . "PL";
				push(@outline , $row_val{$key2go});
#				$key2go = $sm . ":" . "SB";
#				push(@outline , $row_val{$key2go});
				$key2go = $sm . ":" . "GT";
				push(@outline , $row_val{$key2go});
				$key2go = $sm . ":" . "OGT";
				push(@outline , $row_val{$key2go});
				$key2go = $sm . ":" . "PGT";
				push(@outline , exists($row_val{$key2go})?$row_val{$key2go}:"NA" );
				$key2go = $sm . ":" . "PID";
				push(@outline , exists($row_val{$key2go})?$row_val{$key2go}:"NA" );
				if(exists $dng_smps->{$sm}){
					foreach my $d_hd (@dng_cols){
						$key2go = $sm . ":" . $d_hd;
						push(@outline , exists($row_val{$key2go})?$row_val{$key2go}:"NA" );
					}
					if($row_val{$sm . ":DNG_confidence"} eq "High"){
						$dng_conf = 1;
					}
					
				}
				$key2go = $sm . ":" . "DP_TOTAL";
				push(@outline , exists($row_val{$key2go})?$row_val{$key2go}:"NA" );
				$key2go = $sm . ":" . "DP_A";
				push(@outline , exists($row_val{$key2go})?$row_val{$key2go}:"NA" );
				$key2go = $sm . ":" . "DP_C";
				push(@outline , exists($row_val{$key2go})?$row_val{$key2go}:"NA" );
				$key2go = $sm . ":" . "DP_G";
				push(@outline , exists($row_val{$key2go})?$row_val{$key2go}:"NA" );
				$key2go =  $sm . ":" . "DP_T";
				push(@outline , exists($row_val{$key2go})?$row_val{$key2go}:"NA" );
				my $sample_target_distance = exists ($target_info{$target_key}) ? $target_info{$target_key} : "NA";
				push(@outline ,  $sample_target_distance);
			}			
		}
		
		
		foreach my $an (@anno_fields2go){
			push(@outline , $row_val{$an});
		}
		
		#filtering criteria , added May 04, 2016
		my @fil_criteria = ();
		my $typs = (exists $row_val{$type_col}) ? $row_val{$type_col} : die "missing type_col\n";
		my $spx = (exists $row_val{$spx_col}) ? $row_val{$spx_col} : die "missing $spx_col\n";
		my $dbsc1 = (exists $row_val{$dbsc1_col}) ? $row_val{$dbsc1_col} : die "missing dbsc1_col\n";
		my $dbsc2 = (exists $row_val{$dbsc2_col}) ? $row_val{$dbsc2_col} : die "missing dbsc2_col\n";
		my $dist2splice = (exists $row_val{$dist2splice_col}) ? $row_val{$dist2splice_col} : die "missing distance_splicejunction col\n";
		my $spliceAI_passed = test_spliceAI( $row_val{'spliceAI_DS_AG'} , $row_val{'spliceAI_DP_AG'}, $row_val{'spliceAI_DS_AL'}, $row_val{'spliceAI_DP_AL'}, $row_val{'spliceAI_DS_DG'}, $row_val{'spliceAI_DP_DG'} , $row_val{'spliceAI_DS_DL'},
											  $row_val{'spliceAI_DP_DL'} );

		my $passed = 0;
		
		my $typs_passed = 0;
		if($typs  =~ /upstream|downstream|splicing|exonic|UTR3|UTR5/g ){
			$typs_passed = 1;
			push(@fil_criteria , $typs);
		}elsif($spx > 2 || $spx < -2){
			$typs_passed = 1;
			push(@fil_criteria , "spidex");
		}elsif($dbsc1 > 0.6 || $dbsc2 > 0.6){
			$typs_passed = 1;
			push(@fil_criteria , "dbscSNV");
		}elsif($dist2splice ne "NA" &&  abs($dist2splice) <= $spl_thr ){
			$typs_passed = 1;
			push(@fil_criteria , "spl_distance");
		}elsif($spliceAI_passed == 1){
			$typs_passed = 1; 
			push(@fil_criteria , "spliceAI");
		}									  
		
		#if($typs =~ /upstream|downstream|splicing|exonic|UTR3|UTR5/g || $spx > 2 || $spx < -2 || $dbsc1 > 0.6 || $dbsc2 > 0.6 ){
		if($typs_passed){
			$passed = 1;
			foreach (@frqs){
				if($row_val{$_} >= 0.05){
					$passed = 0;
				}
			}
		}
		
		if($dng_conf == 1){
			$passed = 1;
			push(@fil_criteria , "DNG_confidence");
		}
		print $fh_out join("\t" , @outline) , "\n";
		if($passed == 1){
			print $fh_outf join("\t" , @outline) , "\t" , join('|' , @fil_criteria), "\n";
		}
	} #ANN end

	close $fh_ann;
	close $fh_out;
	close $fh_outf;
	$rt_val = 0;
	return $rt_val;
} #end of sub



sub test_spliceAI{
	my ($spAI_DS_AG, $spAI_DP_AG, $spAI_DS_AL, $spAI_DP_AL, $spAI_DS_DG, $spAI_DP_DG, $spAI_DS_DL, $spAI_DP_DL ) = @_;
	my $pass_spAI = 0;
	if( ($spAI_DS_AG > 0.2 && abs($spAI_DP_AG) <= 50) || ($spAI_DS_AL > 0.2 && abs($spAI_DP_AL) <= 50) || ($spAI_DS_DG > 0.2 && abs($spAI_DP_DG)<= 50) || ($spAI_DS_DL > 0.2 && abs($spAI_DP_DL)<= 50) ){
		$pass_spAI = 1;
	}
	return $pass_spAI;
}



sub format_annotationMSSNG{
	my ($ann_file , $smp_file, $type_p, $out_final,$out_filter, $add_cols, $ex_target_dist_out) = @_;
	my $rt_val = -1;
	#open(ANN, $ann_file) or die "$ann_file missing;  $!";
	my $fh_ann = IO::File->new($ann_file , q{<}) or die "$! $ann_file\n";
	#open(SMP, $smp_file) or die "$smp_file missing ;  $!";
	my $fh_smp = IO::File->new($smp_file , q{<}) or die "$! $smp_file\n";
	#open(OUT, ">$out_final") or die "error making $out_final ;  $!";
	my $fh_out = IO::File->new($out_final , q{>}) or die "$! $out_final\n";
	#open(OUTF, ">$out_filter") or die "error making $out_filter ;  $!";
	my $fh_outf = IO::File->new($out_filter , q{>}) or die "$! out_filter\n";
	



##	my @frqs = ("1000g_all" , "1000g_eur" , "1000g_amr" , "1000g_eas" , "1000g_afr" , "1000g_sas" , "NHLBI_all" , "NHLBI_aa" , "NHLBI_eu" , "ExAC_Freq" , "ExAC_AFR" , "ExAC_AMR",
#				 "ExAC_EAS" , "ExAC_FIN" , "ExAC_NFE" , "ExAC_OTH" , "ExAC_SAS" , "cg" , "cg_filtered" , "cgW597_AllFreq" , "cgW597_CalledFreq" ,
#				  "cg1KG436_AllFreq" , "cg1KG436_CalledFreq" );

	my @frqs = ("1000g_all" , "1000g_eur" , "1000g_amr" , "1000g_eas" , "1000g_afr" , "1000g_sas" ,  "ExAC_Freq" , "ExAC_AFR" , "ExAC_AMR",
				 "ExAC_EAS" , "ExAC_FIN" , "ExAC_NFE" , "ExAC_OTH" , "ExAC_SAS" ,
				  "gnomAD_exome_ALL","gnomAD_exome_AFR","gnomAD_exome_AMR","gnomAD_exome_ASJ","gnomAD_exome_EAS","gnomAD_exome_FIN","gnomAD_exome_NFE","gnomAD_exome_OTH","gnomAD_exome_SAS",
				  "gnomAD_genome_ALL","gnomAD_genome_AFR","gnomAD_genome_AMR","gnomAD_genome_ASJ","gnomAD_genome_EAS","gnomAD_genome_FIN","gnomAD_genome_NFE","gnomAD_genome_OTH"   );
				  
	my $type_col = "typeseq";
	my $spx_col = "spx_dpsi";
	my $dbsc1_col = "dbscSNV_ADA_SCORE";
	my $dbsc2_col = "dbscSNV_RF_SCORE";
		
	my @samples = <$fh_smp>;
	chomp(@samples);
	
	my %var_type = ();
	my $header = <$fh_ann>;
	chomp($header);
	my @hds = split("\t" , $header);

	#store index of each header
	my %h_index = ();
	my $ind = 0;
	foreach (@hds){
		$h_index{$_} = $ind++; 
	}
	
	#my @new_hds = ("#CHROM", "start" , "end" , "ref_allele" , "alt_allele", "var_type" , "Original_VCFKEY", "MULTI_ALLELIC" , "FILTER"); # , "SNVSB" , "SNVHPOL" , "RU");
	my @new_hds = ("#MSSNG_key", "Annovar_key" , "var_type");  # , "Original_VCFKEY", "MULTI_ALLELIC" , "FILTER");
	

	
	
	my $typ_index = (exists $h_index{"typeseq"}) ? $h_index{"typeseq"} : die "column \'typeseq\' missing";
	
	
	my @ano_fields = @hds[$typ_index..$#hds];
	my @anno_fields2go = ();
	
	foreach (@ano_fields){
		my $val2go = $_;
		if($val2go ne "spx_dpsi_z"){
			push(@anno_fields2go , $val2go);
			if($val2go eq "typeseq"){
				$val2go = "typeseq";
			}
			if($val2go eq "typeseq_primary"){
				$val2go = "typeseq_priority";
			}
			push(@new_hds , $val2go);
		}
	}

	print $fh_out join("\t", @new_hds) , "\n";
	print $fh_outf join("\t", @new_hds) , "\n";	
	
	while(my $line = <$fh_ann>){
		chomp($line);
		my @fields = split("\t" , $line);
		my $vcf_loc = $fields[0] . ":" . $fields[1] . ":" . $fields[3] . ":" . $fields[4];
		#my $vcf_vartype = (exists $stvc{$vcf_loc}) ? $stvc{$vcf_loc} : "NA";
		my %row_val = ();
		my $hi = 0;
		my @outline = ();
		foreach my $h (@hds){
			my $this_val = $fields[$hi++];
			$this_val =~ s/^\s+//; #remove leading spaces
			$this_val =~ s/\s+$//; #remove trailing spaces
			if($this_val eq '.'){
				$this_val = "NA";
			}
			$row_val{$h} = $this_val;
		}
	
#		my $m_alle = $row_val{'MULTI_ALLELIC'};
#		my $m_order = $row_val{'M_ORDER'};
#		if($m_alle == 0){
#			$m_order = 1;
#		}
		
		my $vcf_vartype = Parse::Annovar::get_annoVar_type($row_val{'ref_allele'},  $row_val{'alt_allele'});
		
		my $mssng_key = $row_val{'MSSNG_KEY'};
		my $annovar_key = $row_val{'annovar_chr'} . "-" .   $row_val{'POS'} . "-" . $row_val{'annovar_end'} . "-" . $row_val{'REF'} . "-" . $row_val{'ALT'};
		
		
		push(@outline , $mssng_key);
		push(@outline , $annovar_key);
#		push(@outline , $row_val{'annovar_end'});
#		push(@outline , $row_val{'REF'});
#		push(@outline , $row_val{'ALT'});
		push(@outline , $vcf_vartype);
		

		
		my $target_key = $row_val{'annovar_chr'} . "\t" . ($row_val{'annovar_start'} - 1) . "\t" . $row_val{'annovar_end'} . "\t"  . $row_val{'ref_allele'} . "\t" . $row_val{'alt_allele'};
		
		foreach my $an (@anno_fields2go){
			push(@outline , $row_val{$an});
		}
		
		#filtering criteria , added May 04, 2016
		my $typs = (exists $row_val{$type_col}) ? $row_val{$type_col} : die "missing $type_col\n";
		#my $spx = (exists $row_val{$spx_col}) ? $row_val{$spx_col} : die "missing $spx_col\n";
		my $dbsc1 = (exists $row_val{$dbsc1_col}) ? $row_val{$dbsc1_col} : die "missing $dbsc1_col\n";
		my $dbsc2 = (exists $row_val{$dbsc2_col}) ? $row_val{$dbsc2_col} : die "missing $dbsc2_col\n";

		my $passed = 0;
		#if($typs =~ /upstream|downstream|splicing|exonic|UTR3|UTR5/g || $spx > 2 || $spx < -2 || $dbsc1 > 0.6 || $dbsc2 > 0.6 ){
		if($typs =~ /upstream|downstream|splicing|exonic|UTR3|UTR5/g || $dbsc1 > 0.6 || $dbsc2 > 0.6 ){
			$passed = 1;
			foreach (@frqs){
				if($row_val{$_} >= 0.05){
					$passed = 0;
				}
			}
		}
		print $fh_out join("\t" , @outline) , "\n";
		if($passed == 1){
			print $fh_outf join("\t" , @outline) , "\n";
		}
		
		
		
	} #ANN end

	close $fh_ann;
	close $fh_out;
	close $fh_outf;
	$rt_val = 0;
	return $rt_val;
} #end of sub








#per row formatting
sub format_annotation_row{
	my ($rownig, $rowval , $header,  $smp_file, $type_p,  $add_cols, $fh_out, $fh_outf, $ex_target_dist_out) = @_;
	my $rt_val = -1;
	#my $fh_ann = IO::File->new($ann_file , q{<}) or die "$! $ann_file\n";
	my $fh_smp = IO::File->new($smp_file , q{<}) or die "$! $smp_file\n";
	#my $fh_out = IO::File->new($out_final , q{>}) or die "$! $out_final\n";
	#my $fh_outf = IO::File->new($out_filter , q{>}) or die "$! out_filter\n";
	
	if($type_p > 4 || $type_p < 0){
		die "Error: $type_p not allowed\n";
	}

#	#coverage file , modified Aug 17,2016
#	my %cov_info = ();
#	my @coverage_header = ();
	
	#added Jan 24,2019
	#For extra columns
	my @extra_cols = split("," , $add_cols);
	
	#target distance for exome 
	my %target_info = ();
	my $fh_ex_target;
	if(($type_p == 3 || $type_p == 3.7) && defined $ex_target_dist_out ){
		warn "exome file, loading target distances..\n";
		$fh_ex_target = IO::File->new($ex_target_dist_out , q{<}) or warn "$! $ex_target_dist_out\n";
		if(defined $fh_ex_target){
			while(my $tr_line = <$fh_ex_target>){
				chomp($tr_line);
				my ($inp , $res) = split('\|' , $tr_line ,2) ;
				my $key = $inp;
				my $dist = (split('\|' , $res))[1];
				$target_info{$key} = $dist;
			}
			close $fh_ex_target;	
		}
	}
	


##	my @frqs = ("1000g_all" , "1000g_eur" , "1000g_amr" , "1000g_eas" , "1000g_afr" , "1000g_sas" , "NHLBI_all" , "NHLBI_aa" , "NHLBI_eu" , "ExAC_Freq" , "ExAC_AFR" , "ExAC_AMR",
#				 "ExAC_EAS" , "ExAC_FIN" , "ExAC_NFE" , "ExAC_OTH" , "ExAC_SAS" , "cg" , "cg_filtered" , "cgW597_AllFreq" , "cgW597_CalledFreq" ,
#				  "cg1KG436_AllFreq" , "cg1KG436_CalledFreq" );

	my @frqs = ("1000g_all" , "1000g_eur" , "1000g_amr" , "1000g_eas" , "1000g_afr" , "1000g_sas" ,  "ExAC_Freq" , "ExAC_AFR" , "ExAC_AMR",
				 "ExAC_EAS" , "ExAC_FIN" , "ExAC_NFE" , "ExAC_OTH" , "ExAC_SAS" ,
				  "gnomAD_exome_ALL","gnomAD_exome_AFR","gnomAD_exome_AMR","gnomAD_exome_ASJ","gnomAD_exome_EAS","gnomAD_exome_FIN","gnomAD_exome_NFE","gnomAD_exome_OTH","gnomAD_exome_SAS",
				  "gnomAD_genome_ALL","gnomAD_genome_AFR","gnomAD_genome_AMR","gnomAD_genome_ASJ","gnomAD_genome_EAS","gnomAD_genome_FIN","gnomAD_genome_NFE","gnomAD_genome_OTH"   );
				  
	my $type_col = "typeseq";
	my $spx_col = "spx_dpsi";
	my $dbsc1_col = "dbscSNV_ADA_SCORE";
	my $dbsc2_col = "dbscSNV_RF_SCORE";
	
		
	my @samples = <$fh_smp>;
	chomp(@samples);
	
	my %var_type = ();
	
	#*** TBA
	#my $header = "something";
	chomp($header);
	my @hds = split("\t" , $header);

	#store index of each header
	my %h_index = ();
	my $ind = 0;
	foreach (@hds){
		$h_index{$_} = $ind++; 
	}
	
	#my @new_hds = ("#CHROM", "start" , "end" , "ref_allele" , "alt_allele", "var_type" , "Original_VCFKEY", "MULTI_ALLELIC" , "FILTER"); # , "SNVSB" , "SNVHPOL" , "RU");
	# my @new_hds = ("#CHROM", "start" , "end" , "ref_allele" , "alt_allele", "var_type");  # , "Original_VCFKEY", "MULTI_ALLELIC" , "FILTER");
	
	# foreach (@extra_cols){
		# push(@new_hds , $_);
	# }

	# push(@new_hds , "Original_VCFKEY");
	# push(@new_hds , "MULTI_ALLELIC");
	# push(@new_hds , "FILTER");
	
	
	# if($type_p == 2){
		# push(@new_hds , "SNVSB");
		# push(@new_hds , "SNVHPOL");
		# push(@new_hds , "RU");
		# foreach my $sm (@samples){
			# my $key2go = $sm . ":" . "Zygosity";
			# push(@new_hds , $key2go);
			# $key2go = $sm . ":" . "Genotype";
			# push(@new_hds , $key2go);
			# $key2go = $sm . ":" . "DP";
			# push(@new_hds , $key2go);
			# $key2go = $sm . ":" . "DPF";
			# push(@new_hds , $key2go);
			# $key2go = $sm . ":" . "DPI";
			# push(@new_hds , $key2go);
			# $key2go = $sm . ":" . "AD_REF";
			# push(@new_hds , $key2go);
			# $key2go = $sm . ":" . "AD_ALT";
			# push(@new_hds , $key2go);
			# $key2go = $sm . ":" . "GQ";
			# push(@new_hds , $key2go);
			# $key2go = $sm . ":" . "GQX";
			# push(@new_hds , $key2go);
			# $key2go = $sm . ":" . "GT_PostNorm";
			# push(@new_hds , $key2go);
			# $key2go = $sm . ":" . "GT_PreNorm";
			# push(@new_hds , $key2go);
		# }
		
	# }elsif($type_p == 1 || $type_p == 3){
		# #"DP" , "FS" , "QD" , "MQ" , "MQRankSum" , "ReadPosRankSum" , "VQSLOD" , "culprit"
		# push(@new_hds , "DP");
		# push(@new_hds , "FS");
		# push(@new_hds , "QD");
		# push(@new_hds , "MQ");
		# push(@new_hds , "GQ_MEAN") unless ($type_p != 1);
		# push(@new_hds , "MQRankSum");
		# push(@new_hds , "ReadPosRankSum");
		# push(@new_hds , "VQSLOD");
		# push(@new_hds , "culprit");
		# foreach my $sm (@samples){
			# my $key2go = $sm . ":" . "Zygosity";
			# push(@new_hds , $key2go);
			# $key2go = $sm . ":" . "Genotype";
			# push(@new_hds , $key2go);
			# $key2go = $sm . ":" . "DP";
			# push(@new_hds , $key2go);
			# $key2go = $sm . ":" . "AD_REF";
			# push(@new_hds , $key2go);
			# $key2go = $sm . ":" . "AD_ALT";
			# push(@new_hds , $key2go);
			# $key2go = $sm . ":" . "GQ";
			# push(@new_hds , $key2go);
			# $key2go = $sm . ":" . "PL";
			# push(@new_hds , $key2go);
			# $key2go = $sm . ":" . "GT_PostNorm";
			# push(@new_hds , $key2go);
			# $key2go = $sm . ":" . "GT_PreNorm";
			# push(@new_hds , $key2go);
			# $key2go = $sm . ":" . "PhaseGT";
			# push(@new_hds , $key2go);
			# $key2go = $sm . ":" . "PhaseID";
			# push(@new_hds , $key2go);
			# if($type_p == 3){
				# #DP_TOTAL=430;DP_A=430;DP_C=0;DP_G=0;DP_T=0
				# $key2go = $sm . ":" . "DP_TOTAL";
				# push(@new_hds , $key2go);
				# $key2go = $sm . ":" . "DP_A";
				# push(@new_hds , $key2go);
				# $key2go = $sm . ":" . "DP_C";
				# push(@new_hds , $key2go);
				# $key2go = $sm . ":" . "DP_G";
				# push(@new_hds , $key2go);
				# $key2go = $sm . ":" . "DP_T";
				# push(@new_hds , $key2go);
				# $key2go = $sm . ":" . "Dist_from_Target";
				# push(@new_hds , $key2go);
				# $key2go = $sm . ":" . "N_count";
				# push(@new_hds , $key2go);
			# }
		# }		
	 # }elsif($type_p == 4){
		# #somatic
		
		
	# }elsif($type_p == 3.7) {
		# #exome 3.7
		
# #11      AC
# #12      AF
# #13      AN
# #14      BaseQRankSum
# #15      ClippingRankSum
# #16      DB
# #17      DP
# #18      DS
# #19      END
# #20      ExcessHet
# #21      FS
# #22      HaplotypeScore
# #23      InbreedingCoeff
# #24      MLEAC
# #25      MLEAF
# #26      MQ
# #27      MQRankSum
# #28      QD
# #29      RAW_MQ
# #30      ReadPosRankSum
# #31      SOR
# #32      DP_TOTAL
# #33      DP_A
# #34      DP_C
# #35      DP_G
# #36      DP_T
# #37      OVCF_KEY
# #38      MULTI_ALLELIC
# #39      BEEL:GTR
# #40      BEEL:AD_REF
# #41      BEEL:AD_ALT
# #42      BEEL:DP
# #43      BEEL:GQ
# #44      BEEL:GT
# #45      BEEL:MIN_DP
# #46      BEEL:PGT
# #47      BEEL:PID
# #48      BEEL:PL
# #49      BEEL:RGQ
# #50      BEEL:SB
# #51      BEEL:OGT
# #52      BEEL:OZYG
# #53      BEEL:Zygosity

		# push(@new_hds , "BaseQRankSum");
		# #push(@new_hds , "ClippingRankSum");
		# push(@new_hds , "DP");
		# #push(@new_hds , "ExcessHet");
		# push(@new_hds , "FS");
		# push(@new_hds , "MQ");
		# push(@new_hds , "MQRankSum");
		# push(@new_hds , "QD");
		# push(@new_hds , "ReadPosRankSum");
		# push(@new_hds , "SOR");
		# foreach my $sm (@samples){
			# my $key2go = $sm . ":" . "Zygosity";
			# push(@new_hds , $key2go);
			# $key2go = $sm . ":" . "Genotype";
			# push(@new_hds , $key2go);
			# $key2go = $sm . ":" . "DP";
			# push(@new_hds , $key2go);
			# $key2go = $sm . ":" . "AD_REF";
			# push(@new_hds , $key2go);
			# $key2go = $sm . ":" . "AD_ALT";
			# push(@new_hds , $key2go);
			# $key2go = $sm . ":" . "GQ";
			# push(@new_hds , $key2go);
			# $key2go = $sm . ":" . "PL";
			# push(@new_hds , $key2go);
			# $key2go = $sm . ":" . "SB";
			# push(@new_hds , $key2go);
			# $key2go = $sm . ":" . "GT_PostNorm";
			# push(@new_hds , $key2go);
			# $key2go = $sm . ":" . "GT_PreNorm";
			# push(@new_hds , $key2go);
			# $key2go = $sm . ":" . "PhaseGT";
			# push(@new_hds , $key2go);
			# $key2go = $sm . ":" . "PhaseID";
			# push(@new_hds , $key2go);
			# $key2go = $sm . ":" . "DP_TOTAL";
			# push(@new_hds , $key2go);
			# $key2go = $sm . ":" . "DP_A";
			# push(@new_hds , $key2go);
			# $key2go = $sm . ":" . "DP_C";
			# push(@new_hds , $key2go);
			# $key2go = $sm . ":" . "DP_G";
			# push(@new_hds , $key2go);
			# $key2go = $sm . ":" . "DP_T";
			# push(@new_hds , $key2go);
			# $key2go = $sm . ":" . "Dist_from_Target";
			# push(@new_hds , $key2go);
		# }		
		
		
	# }
	
	
	my $typ_index = (exists $h_index{"typeseq"}) ? $h_index{"typeseq"} : die "column \'typeseq\' missing";
	
	
	my @ano_fields = @hds[$typ_index..$#hds];
	my @anno_fields2go = ();
	
	foreach (@ano_fields){
		my $val2go = $_;
		if($val2go ne "spx_dpsi_z"){
			push(@anno_fields2go , $val2go);
			if($val2go eq "typeseq"){
				$val2go = "typeseq";
			}
			if($val2go eq "typeseq_primary"){
				$val2go = "typeseq_priority";
			}
			#push(@new_hds , $val2go);
		}
	}

#	print $fh_out join("\t", @new_hds) , "\n";
#	print $fh_outf join("\t", @new_hds) , "\n";	
	
	#while(my $line = <$fh_ann>){
		chomp($rowval);
		my @fields = split("\t" , $rowval);
		my $vcf_loc = $fields[0] . ":" . $fields[1] . ":" . $fields[3] . ":" . $fields[4];
		#my $vcf_vartype = (exists $stvc{$vcf_loc}) ? $stvc{$vcf_loc} : "NA";
		my %row_val = ();
		my $hi = 0;
		my @outline = ();
		foreach my $h (@hds){
			my $this_val = $fields[$hi++];
#			$this_val =~ s/^\s+//; #remove leading spaces
#			$this_val =~ s/\s+$//; #remove trailing spaces
			$this_val =~ s/^\s*(\S*(?:\s+\S+)*)\s*$/$1/;
			if($this_val eq '.'){
				$this_val = "NA";
			}
			$row_val{$h} = $this_val;
		}
	
		my $m_alle = $row_val{'MULTI_ALLELIC'};
#		my $m_order = $row_val{'M_ORDER'};
#		if($m_alle == 0){
#			$m_order = 1;
#		}
		
		my $vcf_vartype = Parse::Annovar::get_annoVar_type($row_val{'ref_allele'},  $row_val{'alt_allele'});
		push(@outline , $row_val{'annovar_chr'});
		push(@outline , $row_val{'POS'});
		push(@outline , $row_val{'annovar_end'});
		push(@outline , $row_val{'REF'});
		push(@outline , $row_val{'ALT'});
		push(@outline , $vcf_vartype);
		
		#extra columns if any
		foreach my $excol (@extra_cols){
			push(@outline , $row_val{$excol});
		}
		push(@outline , $row_val{'OVCF_KEY'});
		push(@outline , $row_val{'MULTI_ALLELIC'});
		#push(@outline , $m_order);	
		push(@outline , $row_val{'FILTER'});
	
		my $ovcf_key = $row_val{'OVCF_KEY'};
#37      ref_allele
#38      alt_allele
		
		my $target_key = $row_val{'annovar_chr'} . "\t" . ($row_val{'annovar_start'} - 1) . "\t" . $row_val{'annovar_end'} . "\t"  . $row_val{'ref_allele'} . "\t" . $row_val{'alt_allele'};
		if($type_p == 2){
			push(@outline , exists($row_val{'SNVSB'})?$row_val{'SNVSB'}:"NA" );
			push(@outline , exists($row_val{'SNVHPOL'})?$row_val{'SNVHPOL'}:"NA" );
			push(@outline , exists($row_val{'RU'})?$row_val{'RU'}:"NA" );

			#sample wise
			foreach my $sm (@samples){
				my $key2go = $sm . ":" . "OZYG";
				push(@outline , $row_val{$key2go});
				$key2go = $sm . ":" . "GTR";
				push(@outline , $row_val{$key2go});
				$key2go = $sm . ":" . "DP";
				push(@outline , $row_val{$key2go});
				$key2go = $sm . ":" . "DPF";
				push(@outline , $row_val{$key2go});
				$key2go = $sm . ":" . "DPI";
				push(@outline , $row_val{$key2go});
				$key2go = $sm . ":" . "AD_REF";
				push(@outline , $row_val{$key2go});
				$key2go = $sm . ":" . "AD_ALT";
				push(@outline , $row_val{$key2go});
				$key2go = $sm . ":" . "GQ";
				push(@outline , $row_val{$key2go});
				$key2go = $sm . ":" . "GQX";
				push(@outline , $row_val{$key2go});
				$key2go = $sm . ":" . "GT";
				push(@outline , $row_val{$key2go});
				$key2go = $sm . ":" . "OGT";
				push(@outline , $row_val{$key2go});
			}
						
		}elsif($type_p == 1 || $type_p == 3){
			push(@outline , exists($row_val{'DP'})?$row_val{'DP'}:"NA" );
			push(@outline , exists($row_val{'FS'})?$row_val{'FS'}:"NA" );
			push(@outline , exists($row_val{'QD'})?$row_val{'QD'}:"NA" );
			push(@outline , exists($row_val{'MQ'})?$row_val{'MQ'}:"NA" );	
			push(@outline , exists($row_val{'GQ_MEAN'})?$row_val{'GQ_MEAN'}:"NA" ) unless ($type_p != 1);								
			push(@outline , exists($row_val{'MQRankSum'})?$row_val{'MQRankSum'}:"NA" );			
			push(@outline , exists($row_val{'ReadPosRankSum'})?$row_val{'ReadPosRankSum'}:"NA" );
			push(@outline , exists($row_val{'VQSLOD'})?$row_val{'VQSLOD'}:"NA" );
			push(@outline , exists($row_val{'culprit'})?$row_val{'culprit'}:"NA" );
			
			#sample wise
			foreach my $sm (@samples){
				my $key2go = $sm . ":" . "OZYG";
				push(@outline , exists ($row_val{$key2go}) ? $row_val{$key2go} : "NA");
				$key2go = $sm . ":" . "GTR";
				#push(@outline , $row_val{$key2go});
				push(@outline , exists ($row_val{$key2go}) ? $row_val{$key2go} : "NA");
				$key2go = $sm . ":" . "DP";
				#push(@outline , $row_val{$key2go});
				push(@outline , exists ($row_val{$key2go}) ? $row_val{$key2go} : "NA");
				$key2go = $sm . ":" . "AD_REF";
				#push(@outline , $row_val{$key2go});
				push(@outline , exists ($row_val{$key2go}) ? $row_val{$key2go} : "NA");
				$key2go = $sm . ":" . "AD_ALT";
				#push(@outline , $row_val{$key2go});
				push(@outline , exists ($row_val{$key2go}) ? $row_val{$key2go} : "NA");
				$key2go = $sm . ":" . "GQ";
				#push(@outline , $row_val{$key2go});
				push(@outline , exists ($row_val{$key2go}) ? $row_val{$key2go} : "NA");
				$key2go = $sm . ":" . "PL";
				#push(@outline , $row_val{$key2go});
				push(@outline , exists ($row_val{$key2go}) ? $row_val{$key2go} : "NA");
				$key2go = $sm . ":" . "GT";
				#push(@outline , $row_val{$key2go});
				push(@outline , exists ($row_val{$key2go}) ? $row_val{$key2go} : "NA");
				$key2go = $sm . ":" . "OGT";
				#push(@outline , $row_val{$key2go});
				push(@outline , exists ($row_val{$key2go}) ? $row_val{$key2go} : "NA");
				$key2go = $sm . ":" . "PGT";
				push(@outline , exists($row_val{$key2go})?$row_val{$key2go}:"NA" );
				$key2go = $sm . ":" . "PID";
				push(@outline , exists($row_val{$key2go})?$row_val{$key2go}:"NA" );
				
				#new Aug 17, 2016
				#modified June 13,2017 ; expects modified vcf; #DP_TOTAL=430;DP_A=430;DP_C=0;DP_G=0;DP_T=0
				if($type_p == 3){
					$key2go = $sm . ":" . "DP_TOTAL";
					push(@outline , exists($row_val{$key2go})?$row_val{$key2go}:"NA" );
					$key2go = $sm . ":" . "DP_A";
					push(@outline , exists($row_val{$key2go})?$row_val{$key2go}:"NA" );
					$key2go = $sm . ":" . "DP_C";
					push(@outline , exists($row_val{$key2go})?$row_val{$key2go}:"NA" );
					$key2go = $sm . ":" . "DP_G";
					push(@outline , exists($row_val{$key2go})?$row_val{$key2go}:"NA" );
					$key2go =  $sm . ":" . "DP_T";
					push(@outline , exists($row_val{$key2go})?$row_val{$key2go}:"NA" );
					my $sample_target_distance = exists ($target_info{$target_key}) ? $target_info{$target_key} : "NA";
					push(@outline ,  $sample_target_distance);
					
#					my $this_read_count = (exists $cov_info{$cov_key}{'read_count'}) ? $cov_info{$cov_key}{'read_count'} : "NA"; 
#					push(@outline , $this_read_count);
#					my $this_read_countA = (exists $cov_info{$cov_key}{'A'}) ? $cov_info{$cov_key}{'A'} : "NA";
#					push(@outline , $this_read_countA);
#					my $this_read_countC = (exists $cov_info{$cov_key}{'C'}) ? $cov_info{$cov_key}{'C'} : "NA";
#					push(@outline , $this_read_countC);
#					my $this_read_countG = (exists $cov_info{$cov_key}{'G'}) ? $cov_info{$cov_key}{'G'} : "NA";
#					push(@outline , $this_read_countG);
#					my $this_read_countT = (exists $cov_info{$cov_key}{'T'}) ? $cov_info{$cov_key}{'T'} : "NA";
#					push(@outline , $this_read_countT);
#					my $this_read_countN = (exists $cov_info{$cov_key}{'N'}) ? $cov_info{$cov_key}{'N'} : "NA";
#					push(@outline , $this_read_countN);
				}
			}
		}elsif($type_p == 4){
			#somatic
		}elsif($type_p == 3.7){
			push(@outline , exists($row_val{'BaseQRankSum'})?$row_val{'BaseQRankSum'}:"NA" );
			#push(@outline , exists($row_val{'ClippingRankSum'})?$row_val{'ClippingRankSum'}:"NA" );
			push(@outline , exists($row_val{'DP'})?$row_val{'DP'}:"NA" );
			#push(@outline , exists($row_val{'ExcessHet'})?$row_val{'ExcessHet'}:"NA" );	
			push(@outline , exists($row_val{'FS'})?$row_val{'FS'}:"NA" ) ;								
			push(@outline , exists($row_val{'MQ'})?$row_val{'MQ'}:"NA" );
			push(@outline , exists($row_val{'MQRankSum'})?$row_val{'MQRankSum'}:"NA" );
			push(@outline , exists($row_val{'QD'})?$row_val{'QD'}:"NA" );
			push(@outline , exists($row_val{'ReadPosRankSum'})?$row_val{'ReadPosRankSum'}:"NA" );
			push(@outline , exists($row_val{'SOR'})?$row_val{'SOR'}:"NA" );
			foreach my $sm (@samples){
				my $key2go = $sm . ":" . "OZYG";
				push(@outline , $row_val{$key2go});
				$key2go = $sm . ":" . "GTR";
				push(@outline , $row_val{$key2go});
				$key2go = $sm . ":" . "DP";
				push(@outline , $row_val{$key2go});
				$key2go = $sm . ":" . "AD_REF";
				push(@outline , $row_val{$key2go});
				$key2go = $sm . ":" . "AD_ALT";
				push(@outline , $row_val{$key2go});
				$key2go = $sm . ":" . "GQ";
				push(@outline , $row_val{$key2go});
				$key2go = $sm . ":" . "PL";
				push(@outline , $row_val{$key2go});
#				$key2go = $sm . ":" . "SB";
#				push(@outline , $row_val{$key2go});
				$key2go = $sm . ":" . "GT";
				push(@outline , $row_val{$key2go});
				$key2go = $sm . ":" . "OGT";
				push(@outline , $row_val{$key2go});
				$key2go = $sm . ":" . "PGT";
				push(@outline , exists($row_val{$key2go})?$row_val{$key2go}:"NA" );
				$key2go = $sm . ":" . "PID";
				push(@outline , exists($row_val{$key2go})?$row_val{$key2go}:"NA" );
				$key2go = $sm . ":" . "DP_TOTAL";
				push(@outline , exists($row_val{$key2go})?$row_val{$key2go}:"NA" );
				$key2go = $sm . ":" . "DP_A";
				push(@outline , exists($row_val{$key2go})?$row_val{$key2go}:"NA" );
				$key2go = $sm . ":" . "DP_C";
				push(@outline , exists($row_val{$key2go})?$row_val{$key2go}:"NA" );
				$key2go = $sm . ":" . "DP_G";
				push(@outline , exists($row_val{$key2go})?$row_val{$key2go}:"NA" );
				$key2go =  $sm . ":" . "DP_T";
				push(@outline , exists($row_val{$key2go})?$row_val{$key2go}:"NA" );
				my $sample_target_distance = exists ($target_info{$target_key}) ? $target_info{$target_key} : "NA";
				push(@outline ,  $sample_target_distance);
			}			
		}
		
		
		foreach my $an (@anno_fields2go){
			push(@outline , $row_val{$an});
		}
		
		#filtering criteria , added May 04, 2016
		my $typs = (exists $row_val{$type_col}) ? $row_val{$type_col} : die "missing $type_col\n";
		#my $spx = (exists $row_val{$spx_col}) ? $row_val{$spx_col} : die "missing $spx_col\n";
		my $dbsc1 = (exists $row_val{$dbsc1_col}) ? $row_val{$dbsc1_col} : die "missing $dbsc1_col\n";
		my $dbsc2 = (exists $row_val{$dbsc2_col}) ? $row_val{$dbsc2_col} : die "missing $dbsc2_col\n";

		my $passed = 0;
		#if($typs =~ /upstream|downstream|splicing|exonic|UTR3|UTR5/g || $spx > 2 || $spx < -2 || $dbsc1 > 0.6 || $dbsc2 > 0.6 ){
		if($typs =~ /upstream|downstream|splicing|exonic|UTR3|UTR5/g || $dbsc1 > 0.6 || $dbsc2 > 0.6 ){
			$passed = 1;
			foreach (@frqs){
				if($row_val{$_} >= 0.05){
					$passed = 0;
				}
			}
		}
		#print "OUT = " , join("\t" , @outline) , "\n";
		print { $$fh_out } join("\t" , @outline) , "\n";
		if($passed == 1){
			print { $$fh_outf } join("\t" , @outline) , "\n";
		}
		
		
		
#	} #ANN end

#	close $fh_out;
#	close $fh_outf;
	$rt_val = 0;
	return $rt_val;
} #end of sub





#Validation for 1st 2M entries
sub run_validation{
	my ($ann_file, $err_file, $no_bgzip) = @_;
	my $fh_ann = undef;
	
	if($no_bgzip != 0){
		$fh_ann = IO::File->new($ann_file , q{<}) or die "$! $ann_file\n";
	}else{
		open($fh_ann, "gunzip -c $ann_file |") or die "gunzip $ann_file: $!";
	}
	my $fh_out = IO::File->new($err_file , q{>}) or die "$! $err_file\n";
	 
	my %validation = ();
	my $ann_header = <$fh_ann>;
	chomp($ann_header);
	my @hds = split("\t" , $ann_header);
	
	my $line_count = 0;
	while(my $line = <$fh_ann> ){
		chomp($line);
		my @fields = split("\t" , $line);
		my $fld_index = 0;
		foreach my $f (@fields){
			if(exists $validation{$hds[$fld_index]}){
				$validation{$hds[$fld_index]} = ($f eq "0" || $f eq "NA") ?  $validation{$hds[$fld_index]} : $f;
			}else{
				$validation{$hds[$fld_index]} = $f;
			}
			$fld_index++;
		}
		$line_count++;
		if($line_count >=200000){
			last;
		}
	}
	close $fh_ann;
	#report
	foreach my $vk (keys %validation){
		if( $validation{$vk} eq "NA" ){
			print  "Warning! $vk has only 0 or  NA values\n";
			print $fh_out "Warning! $vk has only 0 or NA values\n";
		}elsif($validation{$vk} eq "0"){
			print $fh_out "Warning! $vk has only 0 or NA values\n";
		}
	}
	close $fh_out;
	return 1;
}

sub generate_sub_wind_inputs{
	warn "********using sub with run type *********\n";
	my ($in, $out_sub, $out_win, $out_cho, $tmp_dir, $run_type) = @_;
	print "$in, $out_sub, $out_win,$out_cho, $tmp_dir, $run_type\n";
	my %fh;
	my %fhw;
	my %chr_info = ();
	#open(IN, $in) or die "$!\n";
	my $fh_in = IO::File->new($in , q{<}) or die "$! $in\n";
	#open(SUB, ">$out_sub") or die "$!\n";
	my $fh_sub = IO::File->new($out_sub , q{>}) or die "$! $out_sub\n";
	#open(WIN, ">$out_win") or die "$!\n";
	my $fh_win = IO::File->new($out_win , q{>}) or die "$! $out_win\n";
#	open(CHO, ">$out_cho") or die "$!\n";
	my $fh_cho = IO::File->new($out_cho , q{>}) or die "$! $out_cho\n";
    my %processed = {};
    
	while(my $line = <$fh_in>){
		chomp($line);
		if($line =~ m/^#/){ next; }
		$processed{'total'} = (exists $processed{'total'}) ? ($processed{'total'} + 1) : 1;
		my @fields = split("\t", $line);
		my $ref = $fields[3];
		my $obs = $fields[4];
		my $vcf_chr = $fields[0];
		my $vcf_chr0 = $vcf_chr ;
		my $pos_start = $fields[1];
		my $pos_end = $fields[2];
		my $this_ovcf_key = "NA";
		my $multi_index = 0;
#		if($run_type eq "u"){
#			$this_ovcf_key = get_ovcfkey($fields[12]);
#			$multi_index = get_multiIndex($fields[12]);
#		}
		$vcf_chr =~ s/chr//;
		#print $fh_sub $vcf_chr0, "\t" , $pos_start , "\t" , $pos_end , "\t" , $ref , "\t" , $obs , "\t" , $this_ovcf_key , "\t" , $multi_index ,  "\n";
		print $fh_sub $vcf_chr0, "\t" , $pos_start , "\t" , $pos_end , "\t" , $ref , "\t" , $obs , "\n";
		my $win_line = make_window($line,7);
		print $fh_win $win_line, "\n";
#		my $wno_line = make_region($line,2500);
#		print WNO $wno_line , "\n";
		my $chr_filename = $tmp_dir . "/chr_" .  $vcf_chr . ".txt";
		#$chr_info{$vcf_chr0} = $chr_filename;
        if (not $fh{$chr_filename}) {
           open $fh{$chr_filename}, '>', $chr_filename or die "Unable to open '$chr_filename' for output: $!";
           warn "$chr_filename created\n";
        }
        print { $fh{$chr_filename} } $vcf_chr0, "\t" , $pos_start , "\t" , $pos_end , "\t" , $ref , "\t" , $obs ,  "\n"  ;
        $processed{$vcf_chr0} = (exists $processed{$vcf_chr0}) ? ($processed{$vcf_chr0} + 1) : 1;
        		

		my $chrw_filename = $tmp_dir . "/chrw_" .  $vcf_chr . ".txt";
		$chr_info{$vcf_chr0} = $chr_filename . "\t" . $chrw_filename;
#        if (not $fhw{$chrw_filename}) {
#           open $fhw{$chrw_filename}, '>', $chrw_filename or die "Unable to open '$chr_filename' for output: $!";
#           warn "$chrw_filename created\n";
#        }
#        print { $fhw{$chrw_filename} } $win_line , "\n"  ;		
	}

	foreach (keys %chr_info){
		print $fh_cho $_ , "\t" , $chr_info{$_} , "\n" unless ($_ =~ /#/ || $_ eq "MT" );
	}

	close $fh_sub;
	close $fh_win;
	close $fh_in;
	close $fh_cho;	
	return \%processed;
}


sub get_ovcfkey{
	my ($this_info) = @_;
	my $ret = -1;
	if( $this_info =~ m/OVCF_KEY=(.+?);/ || $this_info =~ m/OVCF_KEY=(.+)/ ){
		$ret = $1;
	}else{
		die "Critical error: OVCF_KEY not found in line $this_info\n";
	}
	
	return $ret;
}

sub get_multiIndex{
	my ($this_info) = @_;
	my $ret = -1;
	if( $this_info =~ m/MULTI_ALLELIC=(.+?);/ || $this_info =~ m/MULTI_ALLELIC=(.+)/ ){
		$ret = $1;
	}else{
		die "Critical error: OVCF_KEY not found in line $this_info\n";
	}
	
	return $ret;
}
sub make_window{
	my ($r, $span) = @_;
	my $pos1 = (split("\t", $r))[1];
	my $pos2 = (split("\t", $r))[2];
	my $new_pos1 = $pos1 - $span;
	my $new_pos1_go = ($new_pos1<=0)? 1: $new_pos1;
	my $new_pos2 = $pos2 + $span;
	my $no_of_req_bases = ($new_pos2 - $new_pos1_go) + 1; 
	my $al1 = 'C'x $no_of_req_bases; 
	my $al2 = "A";
	my $new_line = (split("\t", $r))[0] . "\t" . $new_pos1_go . "\t" . $new_pos2 . "\t" . $al1 . "\t" . $al2 . "\t" . (split("\t", $r))[5] ;
	return($new_line);
}

sub make_region{
	my ($r, $span) = @_;
	my $pos1 = (split("\t", $r))[1];
	my $pos2 = (split("\t", $r))[2];
	my $new_pos1 = $pos1 - $span;
        $new_pos1 = ($new_pos1<=0)?1:$new_pos1;
	my $new_pos2 = $pos2 + $span;
	my $no_of_req_bases = ($new_pos2 - $new_pos1) + 1; 
#	my $al1 = 'C'x $no_of_req_bases; 
#	my $al2 = "A";
	my $new_line = (split("\t", $r))[0] . "\t" . $new_pos1 . "\t" . $new_pos2 . "\t" . $pos1 . "\t" . $pos2  ;
	return($new_line);
}

#Not completed
sub post_process_DNG{
	my ($ann_file, $dng_vcf , $orginal_vcf, $ann_dir , $work_dir, $genome) = @_;
	my $base_dng_vcgf = basename($dng_vcf);
	my $base_ann_file = basename($ann_file);
	my $ann_bed = $work_dir .  "/" . $orginal_vcf . ".bedops.bed";
	get_bed($ann_file , $ann_bed);
	
	
	
}

sub get_bed{
	my ($a_file, $b_file) = @_;
	#open(AFL, $a_file) or die "$! $a_file\n";
	my $fh_afl = IO::File->new($a_file , q{<}) or die "$! $a_file\n";
	#open(BFL, ">$b_file") or die "$! $b_file\n";
	my $fh_bfl = IO::File->new($b_file , q{>}) or die "$! $b_file\n";
	
	my $this_hds = <$fh_afl>;
	chomp($this_hds);
	my $idx = 0;
	my %hd_idx = ();
	foreach ( split("\t" , $this_hds) ){
		$hd_idx{$_} = $idx++;
	}
	
	while(my $line = <$fh_afl>){
		chomp($line);
		my @fields = split("\t" , $line);
		my $idx_chr = (exists $hd_idx{'annovar_chr'}) ? $hd_idx{'annovar_chr'} : die "index error for annovar_chr\n ";
		my $idx_start = (exists $hd_idx{'annovar_start'}) ? $hd_idx{'annovar_start'} : die "index error for annovar_start\n ";
		my $idx_end = (exists $hd_idx{'annovar_end'}) ? $hd_idx{'annovar_end'} : die "index error for annovar_end\n ";
		my $idx_ref = (exists $hd_idx{'ref_allele'}) ? $hd_idx{'ref_allele'} : die "index error for ref_allele\n ";
		my $idx_alt = (exists $hd_idx{'alt_allele'}) ? $hd_idx{'alt_allele'} : die "index error for alt_allele\n ";
		print $fh_bfl $fields[$idx_chr] , "\t" , ($fields[$idx_start] - 1) , "\t" , $fields[$idx_end] , "\t" ,   $fields[$idx_chr] . ":" . $fields[$idx_start] . ":" . $fields[$idx_end] . ":" . $fields[$idx_ref] . ":" . $fields[$idx_alt] , "\n";
	}
	close $fh_afl;
	close $fh_bfl;
	return 0;
}

#for exome targets, fixed for running multiple files on same directory
sub add_exome_target_dist{
	my ($fin_file, $target_bed, $odir) = @_;
	
	#Convert target bed to sorted one
	warn "Sorting target bed..\n";
	my $base_target = basename($target_bed);
	my $base_fin = basename($fin_file);
	my $target_bed_sorted = $odir . "/" . $base_fin . "_" .  $base_target . "_sorted.bed";
	my $target_bed_sort_cmd =  "cat $target_bed | "  . q(awk '{print $1"\t"($2)"\t"$3"\t"$1":"$2":"$3}') .  " \| sort-bed - > $target_bed_sorted ";  
	my $ex_stat = system($target_bed_sort_cmd);
	if($ex_stat !=0){
		die "$target_bed_sort_cmd :  \n Error running bedop sorting ...;";
	}
	
	#convert fin to sorted bed and run bedops
	warn "running bedops..\n";
	my @files2del = ();
	
	my $target_bedout = $odir . "/" . $base_fin . "_targets.bed";
	my $fin_out = $odir . "/" . $base_fin . "_withTargetInfo.tsv";
	open(TOUT, ">$fin_out") or die "error making $fin_out $!"; 
	my $bedops_cmd = "cat $fin_file | "  . q(awk '{print $1"\t"($2-1)"\t"$3"\t"$1":"$2":"$3":"$4":"$5}') .  " \| sort-bed - \| closest-features --closest  --dist - $target_bed_sorted > $target_bedout";
	warn "$bedops_cmd\n"; 
	$ex_stat = system($bedops_cmd);
	if($ex_stat !=0){
		die "$bedops_cmd : \n  Error running bedops;";
	}
	
	#chr1    11849446        11849447        chr1:11849447:11849447:A:G|chr1 11849423        11849543        chr1:11849423:11849543|0
	#chr1    11851002        11851003        chr1:11851003:11851003:G:C|chr1 11850759        11850935        chr1:11850759:11850935|-68
	
	warn "Parsing bedops output..\n";
	#parse bedops out
	my %target_info = ();
	open(TAR , $target_bedout) or die "$!";
	while(my $line = <TAR>){
		chomp($line);
		my ($inp , $res) = split('\|' , $line ,2) ;
		my $key = (split("\t" , $inp))[3];
		my $dist = (split('\|' , $res))[1];
		$target_info{$key} = $dist;
	}
	close TAR;
	open(INP , $fin_file) or die "$!";
	my $header = <INP>;
	chomp($header);
	my %hds_a = ();
	my $hi = 0;
	my @hds = split("\t" , $header);
	foreach (split("\t" , $header)){
		$hds_a{$_} = $hi++;
	}
	my $typeseq_hi = (exists $hds_a{'typeseq'}) ? $hds_a{'typeseq'}  : die "typeseq column error in Annotated file\n" ;
	my $name_gt_postnorm = $hds[($typeseq_hi-1)];
	my $sample_name = (split(":" , $name_gt_postnorm))[0];
	
	my @newheader_A = @hds[0..($typeseq_hi-1)];
	my @newheader_B = @hds[($typeseq_hi)..(scalar(@hds)-1)];
	print TOUT join("\t" , @newheader_A) ,  "\t" . $sample_name . ":Dist_from_Target\t" ,  join("\t" , @newheader_B) , "\n";
	
	my $na_targets_passed = 0;
	while(my $line = <INP>){
		chomp($line);
		my @fields = split("\t" , $line);
		my $key = $fields[0] . ":" . $fields[1] . ":" . $fields[2] . ":" . $fields[3] . ":" . $fields[4];
		my $targ_val = (exists $target_info{$key}) ?   $target_info{$key} : die "Targer value for $key is missing\n";
		if($targ_val ne "NA"){
			$na_targets_passed = 1;
		}
		my @fields_A =  @fields[0..($typeseq_hi-1)];
		my @fields_B =  @fields[($typeseq_hi)..(scalar(@fields)-1)];
		print TOUT join("\t" , @fields_A) , "\t" , $targ_val , "\t" , join("\t" , @fields_B) , "\n";
	} 
	unlink($target_bedout);
	unlink($target_bed_sorted);
	warn "Done! Result file : $fin_out\n";
	if($na_targets_passed == 0){
		warn "Warning! All variant has missing target distance(NAs).. (Incompatible target file ?) \n";
	}
}


#added may 15,2017
sub vcf_bed_annotation{
	my ($in_file, $dp_file, $out_v2b_prefix) = @_;
	my $db_bed = $out_v2b_prefix . "_tmp";
	my $ext_st = -8888;
	#open(DPB, ">$db_bed") or die "$!"; 
	my $fh_dpb = IO::File->new($db_bed , q{>}) or die "$! $db_bed\n";
	#convert dp2bed
	#open(COV, $dp_file) or die "$!";
	my $fh_cov= IO::File->new($dp_file , q{>}) or die "$! $dp_file\n";
	
	my $c_header = <$fh_cov>;
	chomp($c_header);
	while(my $line = <$fh_cov>){
		chomp($line);
		my @fields = split("\t" , $line);
		my ($ch, $pos) = split(':' , $fields[0]);
		my @cov2go = ();
		my $read_count = $fields[1];
		foreach (split(' ' , $fields[4])){
			my ($base , $val) = split(':' , $_);
			push(@cov2go , $_) unless ($val == 0);
		}
		print $fh_dpb $ch , "\t" , ($pos - 1 ) , "\t" , $pos , "\t" , $line , "\n"; 
	}
	close $fh_dpb;
	close $fh_cov;
	warn "dp to bed done\n";
	#sort db
	my $anno_out = $out_v2b_prefix; # . "_bed_out.bed"; 
	my $anno_cmd = "vcf2bed < $in_file \| bedmap --delim '\|\|' --skip-unmapped --exact  --echo-map - $db_bed > $anno_out";
	warn "vcf2bed --> $anno_cmd\n";
	$ext_st = system($anno_cmd);
	return $ext_st;
}

#bedmap annovar

sub make_subset_file{
	my ($ann_file , $out_filter, $frq_thr) = @_;
	my $fh_ann = IO::File->new($ann_file , q{<}) or die "Error reading input file \" $ann_file \" ($!)\n";
	#my $fh_smp = IO::File->new($smp_file , q{<}) or die "$! $smp_file\n";
	my $fh_outf = IO::File->new($out_filter , q{>}) or die "Error making output file \" $out_filter \" ($!) \n";
	my @frqs = ("1000g_all" , "1000g_eur" , "1000g_amr" , "1000g_eas" , "1000g_afr" , "1000g_sas" ,  "ExAC_Freq" , "ExAC_AFR" , "ExAC_AMR",
				 "ExAC_EAS" , "ExAC_FIN" , "ExAC_NFE" , "ExAC_OTH" , "ExAC_SAS" ,
				  "gnomAD_exome_ALL","gnomAD_exome_AFR","gnomAD_exome_AMR","gnomAD_exome_ASJ","gnomAD_exome_EAS","gnomAD_exome_FIN","gnomAD_exome_NFE","gnomAD_exome_OTH","gnomAD_exome_SAS",
				  "gnomAD_genome_ALL","gnomAD_genome_AFR","gnomAD_genome_AMR","gnomAD_genome_ASJ","gnomAD_genome_EAS","gnomAD_genome_FIN","gnomAD_genome_NFE","gnomAD_genome_OTH"   );

	my @frqsOld = ("1000g_all" , "1000g_eur" , "1000g_amr" , "1000g_eas" , "1000g_afr" , "1000g_sas" , "NHLBI_all" , "NHLBI_aa" , "NHLBI_eu" , "ExAC_Freq" , "ExAC_AFR" , "ExAC_AMR",
				 "ExAC_EAS" , "ExAC_FIN" , "ExAC_NFE" , "ExAC_OTH" , "ExAC_SAS" , "cg" , "cg_filtered" , "cgW597_AllFreq" , "cgW597_CalledFreq" ,
				  "cg1KG436_AllFreq" , "cg1KG436_CalledFreq" );
				  
	my $type_col = "typeseq";
	my $spx_col = "spx_dpsi";
	my $dbsc1_col = "dbscSNV_ADA_SCORE";
	my $dbsc2_col = "dbscSNV_RF_SCORE";
#	my @samples = <$fh_smp>;
#	chomp(@samples);

	my $header = <$fh_ann>;
	chomp($header);
	my @hds = split("\t" , $header);

	#store index of each header
	my %h_index = ();
	my $ind = 0;
	my $anno_ver = 1;
	foreach (@hds){
		$h_index{$_} = $ind++;
		if($_ eq "cg_filtered"){
			warn "Detected old version of annotation!\n";
			warn "Following frequecy columns are expected:\n" ,  join("\n" , @frqsOld) , "\n\n";
			$anno_ver = 0;
		} 
	}
	my $typ_index = (exists $h_index{"typeseq"}) ? $h_index{"typeseq"} : die "column \'typeseq\' missing";
	print $fh_outf $header , "\n";
	my $fil_count = 0;
	while(my $line = <$fh_ann>){
		chomp($line);
		my @fields = split("\t" , $line);
		my $passed = 0;
		my $typs = $fields[$h_index{"typeseq"}];
		my $spx = $fields[$h_index{$spx_col}];
		my $dbsc1 = $fields[$h_index{$dbsc1_col}];
		my $dbsc2 = $fields[$h_index{$dbsc2_col}];
		
		my @act_freq = ($anno_ver == 1) ? @frqs : @frqsOld;
		if($typs =~ /upstream|downstream|splicing|exonic|UTR3|UTR5/g || $spx > 2 || $spx < -2 || $dbsc1 > 0.6 || $dbsc2 > 0.6 ){
			$passed = 1;
			foreach (@act_freq){
				if($fields[$h_index{$_}] >= $frq_thr){
					$passed = 0;
				}
			}
		}
		if($passed == 1){
			$fil_count++;
			print $fh_outf $line , "\n";
		}
	}
	print STDOUT $fil_count , " entries written to $out_filter\n";
}

#with splicing threshold 
sub make_subset_file_ext{
	my ($ann_file , $out_filter, $frq_thr, $spl_thr) = @_;
	my $fh_ann = IO::File->new($ann_file , q{<}) or die "Error reading input file \" $ann_file \" ($!)\n";
	#my $fh_smp = IO::File->new($smp_file , q{<}) or die "$! $smp_file\n";
	my $fh_outf = IO::File->new($out_filter , q{>}) or die "Error making output file \" $out_filter \" ($!) \n";
	my @frqs = ("1000g_all" , "1000g_eur" , "1000g_amr" , "1000g_eas" , "1000g_afr" , "1000g_sas" ,  "ExAC_Freq" , "ExAC_AFR" , "ExAC_AMR",
				 "ExAC_EAS" , "ExAC_FIN" , "ExAC_NFE" , "ExAC_OTH" , "ExAC_SAS" ,
				  "gnomAD_exome_ALL","gnomAD_exome_AFR","gnomAD_exome_AMR","gnomAD_exome_ASJ","gnomAD_exome_EAS","gnomAD_exome_FIN","gnomAD_exome_NFE","gnomAD_exome_OTH","gnomAD_exome_SAS",
				  "gnomAD_genome_ALL","gnomAD_genome_AFR","gnomAD_genome_AMR","gnomAD_genome_ASJ","gnomAD_genome_EAS","gnomAD_genome_FIN","gnomAD_genome_NFE","gnomAD_genome_OTH"   );

	my @frqsOld = ("1000g_all" , "1000g_eur" , "1000g_amr" , "1000g_eas" , "1000g_afr" , "1000g_sas" , "NHLBI_all" , "NHLBI_aa" , "NHLBI_eu" , "ExAC_Freq" , "ExAC_AFR" , "ExAC_AMR",
				 "ExAC_EAS" , "ExAC_FIN" , "ExAC_NFE" , "ExAC_OTH" , "ExAC_SAS" , "cg" , "cg_filtered" , "cgW597_AllFreq" , "cgW597_CalledFreq" ,
				  "cg1KG436_AllFreq" , "cg1KG436_CalledFreq" );
				  
	my $type_col = "typeseq";
	my $spx_col = "spx_dpsi";
	my $dbsc1_col = "dbscSNV_ADA_SCORE";
	my $dbsc2_col = "dbscSNV_RF_SCORE";
	#distance_spliceJunction
	my $distspl_col = "distance_spliceJunction";
#	my @samples = <$fh_smp>;
#	chomp(@samples);
    warn "**Making SUBSET file ...\n";
	my $header = <$fh_ann>;
	chomp($header);
	my @hds = split("\t" , $header);

	#store index of each header
	my %h_index = ();
	my $ind = 0;
	my $anno_ver = 1;
	foreach (@hds){
		$h_index{$_} = $ind++;
		if($_ eq "cg_filtered"){
			warn "Detected old version of annotation!\n";
			warn "Following frequecy columns are expected:\n" ,  join("\n" , @frqsOld) , "\n\n";
			$anno_ver = 0;
		} 
	}
	my $typ_index = (exists $h_index{"typeseq"}) ? $h_index{"typeseq"} : die "column \'typeseq\' missing";
	my $distspl_index = (exists $h_index{$distspl_col}) ? $h_index{$distspl_col} : die "column \' $distspl_col \' missing";
	print $fh_outf $header , "\n";
	my $fil_count = 0;
	while(my $line = <$fh_ann>){
		chomp($line);
		my @fields = split("\t" , $line);
		my $passed = 0;
		my $typs = $fields[$h_index{"typeseq"}];
		my $spx = $fields[$h_index{$spx_col}];
		my $dbsc1 = $fields[$h_index{$dbsc1_col}];
		my $dbsc2 = $fields[$h_index{$dbsc2_col}];
		my $this_distSp = $fields[$distspl_index];
		
		my @act_freq = ($anno_ver == 1) ? @frqs : @frqsOld;
		#print $this_distSp , " -- " ,  abs($this_distSp) , "\n";
		if( ($typs =~ /upstream|downstream|splicing|exonic|UTR3|UTR5/g || $spx > 2 || $spx < -2 || $dbsc1 > 0.6 || $dbsc2 > 0.6) || ($this_distSp ne "NA" &&  abs($this_distSp) <= $spl_thr ) ){
			$passed = 1;
			foreach (@act_freq){
				if($fields[$h_index{$_}] >= $frq_thr){
					$passed = 0;
				}
			}
		}
		if($passed == 1){
			$fil_count++;
			print $fh_outf $line , "\n";
		}
	}
	print STDOUT $fil_count , " entries written to $out_filter\n";
}

sub make_Rare_file{
	my ($ann_file , $out_filter, $frq_thr) = @_;
	my $fh_ann = IO::File->new($ann_file , q{<}) or die "Error reading input file \" $ann_file \" ($!)\n";
	#my $fh_smp = IO::File->new($smp_file , q{<}) or die "$! $smp_file\n";
	my $fh_outf = IO::File->new($out_filter , q{>}) or die "Error making output file \" $out_filter \" ($!) \n";
	my @frqs = ("1000g_all" , "1000g_eur" , "1000g_amr" , "1000g_eas" , "1000g_afr" , "1000g_sas" ,  "ExAC_Freq" , "ExAC_AFR" , "ExAC_AMR",
				 "ExAC_EAS" , "ExAC_FIN" , "ExAC_NFE" , "ExAC_OTH" , "ExAC_SAS" ,
				  "gnomAD_exome_ALL","gnomAD_exome_AFR","gnomAD_exome_AMR","gnomAD_exome_ASJ","gnomAD_exome_EAS","gnomAD_exome_FIN","gnomAD_exome_NFE","gnomAD_exome_OTH","gnomAD_exome_SAS",
				  "gnomAD_genome_ALL","gnomAD_genome_AFR","gnomAD_genome_AMR","gnomAD_genome_ASJ","gnomAD_genome_EAS","gnomAD_genome_FIN","gnomAD_genome_NFE","gnomAD_genome_OTH"   );

	my @frqsOld = ("1000g_all" , "1000g_eur" , "1000g_amr" , "1000g_eas" , "1000g_afr" , "1000g_sas" , "NHLBI_all" , "NHLBI_aa" , "NHLBI_eu" , "ExAC_Freq" , "ExAC_AFR" , "ExAC_AMR",
				 "ExAC_EAS" , "ExAC_FIN" , "ExAC_NFE" , "ExAC_OTH" , "ExAC_SAS" , "cg" , "cg_filtered" , "cgW597_AllFreq" , "cgW597_CalledFreq" ,
				  "cg1KG436_AllFreq" , "cg1KG436_CalledFreq" );
				  
	my $type_col = "typeseq";
	my $spx_col = "spx_dpsi";
	my $dbsc1_col = "dbscSNV_ADA_SCORE";
	my $dbsc2_col = "dbscSNV_RF_SCORE";
#	my @samples = <$fh_smp>;
#	chomp(@samples);

	my $header = <$fh_ann>;
	chomp($header);
	my @hds = split("\t" , $header);

	#store index of each header
	my %h_index = ();
	my $ind = 0;
	my $anno_ver = 1;
	foreach (@hds){
		$h_index{$_} = $ind++;
		if($_ eq "cg_filtered"){
			warn "Detected old version of annotation!\n";
			warn "Following frequecy columns are expected:\n" ,  join("\n" , @frqsOld) , "\n\n";
			$anno_ver = 0;
		} 
	}
	my $typ_index = (exists $h_index{"typeseq"}) ? $h_index{"typeseq"} : die "column \'typeseq\' missing";
	print $fh_outf $header , "\n";
	my $fil_count = 0;
	while(my $line = <$fh_ann>){
		chomp($line);
		my @fields = split("\t" , $line);
		my $passed = 0;
		my $typs = $fields[$h_index{"typeseq"}];
		my $spx = $fields[$h_index{$spx_col}];
		my $dbsc1 = $fields[$h_index{$dbsc1_col}];
		my $dbsc2 = $fields[$h_index{$dbsc2_col}];
		
		my @act_freq = ($anno_ver == 1) ? @frqs : @frqsOld;
		#if($typs =~ /upstream|downstream|splicing|exonic|UTR3|UTR5/g || $spx > 2 || $spx < -2 || $dbsc1 > 0.6 || $dbsc2 > 0.6 ){
			$passed = 1;
			foreach (@act_freq){
				if( $fields[$h_index{$_}] >= $frq_thr){
					$passed = 0 unless $fields[$h_index{$_}] eq "NA";
				}
			}
		#}
		if($passed == 1){
			$fil_count++;
			print $fh_outf $line , "\n";
		}
	}
	print STDOUT $fil_count , " entries written to $out_filter\n";
}

1;
