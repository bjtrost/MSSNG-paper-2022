#!/usr/bin/perl
use strict;
use IO::File;
use IO::Zlib;

our $AUTHOR =   '$Author: Thomas Nalpathamkalam <thomas.nalpathamkalam@sickkids.ca> $';

## Module to process and extract information from vcf files, called by the main script 

package Vcfparse;
use File::Basename;


sub new{
	my $class = shift;
	my $filename = shift;
	my $self = {};
	bless ($self, $class);
	$self->{"name"} = $filename;
	$self->{"samples"} = undef;
	$self->{"chr"} = undef;
	$self->{"version"} = "Vcfparse8";
	return $self; 
}

sub get_moduleInfo{
	my $self = shift;
	my @out = ();
	push(@out , $self->{"version"});
	push(@out , $self->{"name"});
	return @out;
}

sub get_samples{
	my $self = shift;
	my $vcffile = $self->{"name"};
	#my $sample_command = "vcf-query -l  $vcffile";
	my $sample_command = "bcftools query -l  $vcffile";
	print STDOUT $sample_command , "\n";
	my @y = `$sample_command`;
	chomp(@y);
	$self->{"samples"} = join("," , @y);
	return @y;
}


##obse
sub get_header{
	my $self = shift;
	my $vcffile = $self->{"name"};
	my $vcfh_command = "bcftools view -h $vcffile";
	my @ex_hds = `$vcfh_command`;
	chomp(@ex_hds);
	return @ex_hds;
}


sub get_zygosity{
	my $self = shift;
	my $gtr = shift;
	my $this_zygosity = "NA";
	if ( $gtr =~ m#^(\d)/(\d)# ) {
		if ( $1 == $2 ) {
			$this_zygosity = ($1==0 && $2 == 0) ? "hom-ref" : "hom-alt";
		#	$this_zygosity = "hom";
		}
		elsif ( $1 == 0 || $2 == 0 ) {
			$this_zygosity = "ref-alt";
		}
		else {
				$this_zygosity = "alt-alt";
		}
		}else {
				$this_zygosity = "unknown";
		}
		return $this_zygosity;
	
}

sub get_ads{
	my $self = shift;
	my $als = shift;
	my $gtr = shift;
	my $ad = shift;
	my %sm_parsed = ();
	my @proper_als = ((split("," , $als))[0]);
	my $ref_al = (split("," , $als))[0];
	my $ref_ad = (split("," , $ad))[0];
	
	my $this_zygosity = -1;
	my @ad_pattern = (-1,-1);
	if ( $gtr =~ m#^(\d)/(\d)# ) {
		if ( $1 == $2 ) {
			$this_zygosity = ($1 == 0) ? 0 : 1 ;
			@ad_pattern = ($1 == 0) ? (0,0) : (0,1);
			my $this_allele = $ref_al . ":" . (split("," , $als))[$1];
			$sm_parsed{$this_allele} = $ref_ad . "\t" . (split("," , $ad))[$ad_pattern[0]];
			$this_allele = $ref_al . ":" . (split("," , $als))[$2];
			$sm_parsed{$this_allele} = $ref_ad . "\t" . (split("," , $ad))[$ad_pattern[1]];
		}
		elsif ( $1 == 0 || $2 == 0 ) {
			$this_zygosity = 2;
			@ad_pattern = ($1 == 0) ? (0,1) : (1,0);
			my $this_allele = $ref_al . ":" . (split("," , $als))[$1];
			$sm_parsed{$this_allele} = $ref_ad . "\t" . (split("," , $ad))[$ad_pattern[0]];
			$this_allele = $ref_al . ":" . (split("," , $als))[$2];
			$sm_parsed{$this_allele} = $ref_ad . "\t" . (split("," , $ad))[$ad_pattern[1]];
		}
		else {
			$this_zygosity = 3;
			@ad_pattern = (0,1,2);
			my $this_allele = $ref_al . ":" . (split("," , $als))[$1];
			$sm_parsed{$this_allele} = $ref_ad . "\t" . (split("," , $ad))[$ad_pattern[1]];
			 $this_allele = $ref_al . ":" . (split("," , $als))[$2];
			$sm_parsed{$this_allele} = $ref_ad . "\t" . (split("," , $ad))[$ad_pattern[2]];
		}
		}else {
			$this_zygosity = -1;
		}
		#no chance of -1
		return %sm_parsed;
}

##############updates Jul 30,2019#############


sub reformat_vcf_simple{
	#formatting without merging text information
	my $self = shift;
	my $vcffile = $self->{"name"};
	my $infos = shift;
	my $formats = shift;
	my $parsed_out = shift;
	my $denovo_gear = shift;
	my $annovar_home = shift;
	my @info_filters   = ();
	my @format_filters = ();
	my $info_fls = $infos;
	my %rfmt_stats = ();
	my %count_alts = ();
	my $vdir  = dirname($parsed_out);
	my $prefix = basename($vcffile);


	print STDOUT "Converting to annovar format...\n";
	my $this_annovar_in = $vdir .  "/" . $prefix . "_annovar.in.txt";
	my $c2annovar_log = $vdir .  "/" . $prefix . "_c2annovar.log";
	my $c2a_command = "perl $annovar_home\/convert2annovar.pl -format vcf4 -allsample -withfreq -include -outfile $this_annovar_in $vcffile 2> $c2annovar_log";
	print STDOUT "Convert2annovar: $c2a_command\n";
	my $ex_st = system($c2a_command);
	if ( $ex_st == 0 ) {
		print STDOUT "** convert2annovar is finished  with exit code: $?\n";
	}
	else {
		warn "convert2annovar is failed with exit code: $?\n";
		usage("Critical:convert2annovar is failed with exit code: $?");
	}
	
	#check any warning
	my $warn_c2annovar =  `grep "WARN" $c2annovar_log`;
	my $fh_ann = IO::File->new($this_annovar_in , q{<}) or die "$! $this_annovar_in\n";
	my $fh_vcf = IO::File->new($vcffile , q{<}) or die "$! $vcffile\n";
	my $fh_rout = IO::File->new($parsed_out , q{>}) or die "$! $parsed_out\n";	
	print STDOUT "Printing headers for input file..\n";
	while ( my $line = <$fh_vcf> ) {
		chomp($line);
		if ( $line =~ m/^#/ ) {
			if ( $line =~ m/^#CHROM/ ) {
				print $fh_rout "#ann_chr\tann_start\tann_end\tann_ref\tann_alt\t" ,  $line, "\n";
			}
		}else{
			last;
		}
	}
	

	print STDOUT "Populating input file..\n";
	my $fh_ann_count = 0;
	
	my $this_ref = ".";
	my %count_zyg = ();
	
	while ( my $line = <$fh_ann>){
	 	chomp($line);
	 	$fh_ann_count++;	
		my @fields = split( "\t", $line );
		my $alt = $fields[12];
		my $line_to_go =
		    $fields[0] . "\t"
			  . $fields[1] . "\t"
			  . $fields[2] . "\t"
			  . $fields[3] . "\t"
			  . $fields[4]  . "\t"
			  . join( "\t", @fields[ 8 .. ( scalar @fields ) - 1 ] ) ;
				  
			  print $fh_rout   $line_to_go, "\n";				  
	 	
	}  #my $line = <$fh_ann>
	
	$rfmt_stats{'var-counts'} = $fh_ann_count;
	$rfmt_stats{'c2annovar_warn'} = $warn_c2annovar;
	foreach (keys %count_zyg){
		$rfmt_stats{$_} = $count_zyg{$_};
	}

	foreach (keys %count_alts){
		$rfmt_stats{$_} = $count_alts{$_};
	}	

	close $fh_rout;
	close $fh_vcf;

	return \%rfmt_stats;	
	
} #end of sub


sub reformat_vcf{
	my $self = shift;
	my $vcffile = $self->{"name"};
	my $infos = shift;
	my $formats = shift;
	my $parsed_out = shift;
	my $denovo_gear = shift;
	my $annovar_home = shift;
	my $dng_meta = shift;
	my @info_filters   = ();
	my @format_filters = ();
	my $info_fls = $infos;
	my %rfmt_stats = ();
	my %count_alts = ();
	
	my $vdir  = dirname($parsed_out);
	my $prefix = basename($vcffile);
	my $tmp_out = $vdir .  "/" . $prefix . "_parsed.tsv";

	my $parse_command =
		  "bcftools query -f \'\%CHROM\\t\%POS\\t\%REF\\t\%ALT\\t\%QUAL\\t%FILTER";
	
		foreach my $inf (split("," , $infos)) {
			push( @info_filters, "\\t\%INFO\/$inf" );
		}
		
		push( @format_filters, "\%TGT" ) unless ($denovo_gear == 1);
		foreach my $inf (split(",",$formats)) {
			push( @format_filters, "\\t\%$inf" );
		}
		
		my $format_fls = join( ",", @format_filters );

		$parse_command =
		    $parse_command
		  . join( "", @info_filters )
		  . "\[\\t\{"
		  . join( "", @format_filters )
		  . "}\]\n\' $vcffile > $tmp_out";

		print STDOUT "Parse multi-sample vcf file: \n";
		print STDOUT "$parse_command\n";
		my $ex_st = system($parse_command);
		if ( $ex_st == 0 ) {
			print STDOUT "** parse-vcf is finished  with exit code: $?\n";
		}
		else {
			warn "parse-vcf is failed with exit code: $?\n";
			usage("Critical:parse-vcf is failed with exit code: $?");
		}

		my $org_vcf = $vcffile;
		my $fil_vcf = $tmp_out;
		my $rfmt_out = $parsed_out;
		my @all_samples = $self->get_samples();
		
		my @info_fields = $self->parse_INFO_fields();
		print "INFOs: " , join("," , @info_fields) , "\n";
		my @format_fields = ();
		@format_fields = ("GTR") unless ($denovo_gear == 1);
		push (@format_fields, $self-> parse_FORMAT_fields());

		my $fh_vcf = IO::File->new($org_vcf , q{<}) or die "$! $org_vcf\n";
		my $fh_fil = IO::File->new($fil_vcf , q{<}) or die "$! $fil_vcf\n";
		

		print  "Scanning for multi-allelic sites...\n";
		my %mult_alle_vars	= ();
		my $count_malleles = 0;
		while ( my $line = <$fh_fil> ) {
			chomp($line);
			my @fields = split( "\t", $line );
			my @fields_tmp = @fields[ 6 .. ( scalar @fields ) - 1 ];
			my @infos = @fields_tmp[ 0 .. ( scalar @info_fields ) - 1 ];
		
			my %info_hash = ();
			while ( my ( $index, $value ) = each @info_fields ) {
				$info_hash{$value} = $infos[$index];
			}
		
			my $mal_flag = $info_hash{'MULTI_ALLELIC'};
			if($mal_flag != 0){
				my $ovcf_key = $info_hash{'OVCF_KEY'};
				$count_malleles++;
				#print "MALF: $mal_flag , $ovcf_key \n";
				$mult_alle_vars{$ovcf_key} = (exists $mult_alle_vars{$ovcf_key}) ? $mult_alle_vars{$ovcf_key} . "\n" . ($mal_flag . "\t" . $line)  : ($mal_flag . "\t" . $line);
			}
		}
		close $fh_fil;		
		
		$fh_fil = IO::File->new($fil_vcf , q{<}) or die "$! $fil_vcf\n";
		
		print STDOUT "Converting to annovar format...\n";
		my $this_annovar_in = $vdir .  "/" . $prefix . "_annovar.in.txt";
		my $c2annovar_log = $vdir .  "/" . $prefix . "_c2annovar.log";
		my $c2a_command = "perl $annovar_home\/convert2annovar.pl -format vcf4 -allsample -withfreq -include -outfile $this_annovar_in $org_vcf 2> $c2annovar_log";
		print STDOUT "Convert2annovar: $c2a_command\n";
		$ex_st = system($c2a_command);
		if ( $ex_st == 0 ) {
			print STDOUT "** convert2annovar is finished  with exit code: $?\n";
		}
		else {
			warn "convert2annovar is failed with exit code: $?\n";
			usage("Critical:convert2annovar is failed with exit code: $?");
		}
	
		#check any warning
		my $warn_c2annovar =  `grep "WARN" $c2annovar_log`;
		my $fh_ann = IO::File->new($this_annovar_in , q{<}) or die "$! $this_annovar_in\n";


		#dng if available
		my %dng_info = ();
		my @dng_smps = ();
		my %dng_chksmp = ();
		my @dng_cols = ("NULL_CONFIG","ML_NULL","PP_NULL","DNM_CONFIG","ML_DNM","PP_DNM","MQ","RD","pair_null_code","pair_denovo_code");
		if($dng_meta ne "NA"){
			
			
			open(DMTA, $dng_meta) or die "$! dng_meta : $dng_meta\n";
			while(my $dmt_line = <DMTA>){
				chomp($dmt_line);
				my ($dn_smp, $dn_fl) = split("\t" , $dmt_line);
				push(@dng_smps , $dn_smp);
				$dng_chksmp{$dn_smp} = 1;
				open(TDG, $dn_fl) or die "$!  $dn_fl \n";
				my $dn_hd = <TDG>;
				my %dn_hds_idx = ();
				my $hds_i = 0;
				foreach (split("\t" , $dn_hd)){
					$dn_hds_idx{$_}  = $hds_i++;
				}
				
				while(my $tdg_line = <TDG>){
					chomp($tdg_line);
					my @tdg_fields = split("\t" , $tdg_line);
					my $tdg_key = $tdg_fields[0] . ":" . $tdg_fields[1] . ":" . $tdg_fields[3] . ":" . $tdg_fields[4];
					foreach my $this_dn_col (@dng_cols){
						my $sm_key = $dn_smp . ":" . $this_dn_col;
						$dng_info{$dn_smp}{$tdg_key}{$sm_key} = (exists $dn_hds_idx{$sm_key}) ? $tdg_fields[$dn_hds_idx{$sm_key}] : "NA";
					}
				}
				close TDG;
				
			}
			
		} #NA
		

		
		print STDOUT "Printing headers for input file..\n";
		my $fh_rout = IO::File->new($rfmt_out , q{>}) or die "$! $rfmt_out\n";
		my @new_head = ();	
		my %dng_tracking = ();	
		while ( my $line = <$fh_vcf> ) {
			chomp($line);
			if ( $line =~ m/^#/ ) {
			###############
				if ( $line =~ m/^#CHROM/ ) {
					my $samples = ( split( "\tFORMAT\t", $line ) )[1];
					
					$samples =~ s/\s+/,/g;
					@all_samples = split( ",", $samples );
					print $fh_rout "#ann_chr\tann_start\tann_end\tann_ref\tann_alt\t" , $line, "\t"; ;
					
					push( @new_head, join( "\t", @info_fields ) );
					my $smp_track_count = 0;
					foreach my $k (@all_samples) {
						$smp_track_count++;
						my $new_h = join( "\t", @format_fields ) . "\tZygosity";
						$new_h =~ s/AD/AD_REF\tAD_ALT/;
						foreach my $nh ( split( "\t", $new_h ) ) {
##							my $format_prefix = ( exists $fam{$k} ) ? $fam{$k} : $k;
							my $h = $k . ":" . $nh;
							push( @new_head, $h );
						}

						if(exists $dng_chksmp{$k}){
							$dng_tracking{$smp_track_count} = 1;
							foreach my $d_col (@dng_cols){
								push( @new_head, $k . ":" . $d_col );
							}
							push(@new_head, $k . ":DNG_confidence" )
						}

						
					}
					print $fh_rout  join( "\t", @new_head ), "\n";

				}

			###############
			}else{
				last;
			}
		
		}
		
		
		#print $fh_rout  join( "\t", @new_head ), "\n";
		


		print STDOUT "Populating input file..\n";
		
		my $fh_ann_count = 0;
		
		my $this_ref = ".";
		my %count_zyg = ();
		
		while (  (my $line = <$fh_ann>) and (my $line_fil = <$fh_fil>) ) {
			 	chomp($line);
			 	$fh_ann_count++;	
				my @fields = split( "\t", $line );
				my $alt = $fields[12];
				my @fields_fil = split( "\t", $line_fil );
				my $alt_fil = $fields_fil[3];
				my @sample_values = ();
				
				#making sure that both files are in sync
				#MSSNG_KEY=chr10-40999016-40999017-G-C;OVCF_KEY=chr10:40999017:G:C;MULTI_ALLELIC=0
				my $line_ovcf = -1;
				if($fields[15] =~ m/OVCF_KEY=(.+?);/ ){
					$line_ovcf = $1;
				}else{
					usage("Critical: Cannot parse OVCF_KEY from line: $line\n");
				}
				
				

				my $line_key = $fields[8] . ":" . $fields[9] . ":"  . $fields[11] . ":" . $fields[12] . ":" . $line_ovcf;
				
				
				my @fields_fil_tmp = @fields_fil[ 6 .. ( scalar @fields_fil ) - 1 ];
				my $field_fil_tab  = join( "\t", @fields_fil_tmp );

				my @sample_formats = ();
				my %info_hash = ();
				my @infos0 = @fields_fil_tmp[ 0 .. ( scalar @info_fields ) - 1 ];
				
				my @infos = ();
				#format fields now , added Sep13,2019
				#$this_val =~ s/^\s+//; #remove leading spaces
				#$this_val =~ s/\s+$//; #remove trailing spaces
				foreach my $this_inf (@infos0){
					$this_inf =~ s/^\s+//;
					$this_inf =~ s/\s+$//;
					$this_inf = ($this_inf eq '.') ? "NA" :  $this_inf;
					push(@infos, $this_inf);
				}
				
				
				while ( my ( $index, $value ) = each @info_fields ) {
					$info_hash{$value} = $infos[$index];
				}
				
				my $mal_flag = $info_hash{'MULTI_ALLELIC'};
				
				my $line_fil_key    =  $fields_fil[0] . ":" . $fields_fil[1] . ":"  . $fields_fil[2] . ":" . $fields_fil[3] . ":" . $info_hash{'OVCF_KEY'};
				
				#print "*** $line_fil_key --- $line_key\n";
				
				if($line_key ne $line_fil_key){
					usage ("Critical: lines are not in sync: ( ** $line_key  ===  $line_fil_key)\nAnnovar = $line\nParsed = $line_fil");
				}
				
				
				

				while ( $field_fil_tab =~ m/\{(.+?)\}/g ) {
					push( @sample_formats, $1 );
				}
				
				my $alt_count = 1;

				if($alt_fil eq "*"){
					$count_alts{$alt_fil} = (exists $count_alts{$alt_fil}) ? ($count_alts{$alt_fil} + 1) : 1;
				}
				
				my $smp_track_count = 0;
				my $dng_line_key = $fields[8] . ":" . $fields[9] . ":"  . $fields[11] . ":" . $fields[12];
				my $this_annovar_type = get_annoVar_type($fields[11] , $fields[12]);
				my @dng_smps_tmp = @dng_smps;
				foreach my $this_format (@sample_formats) {
					$smp_track_count++;
					$this_format =~ s/^\s+//;    #remove leading spaces
					$this_format =~ s/\s+$//;    #remove trailing spaces
					my @format_vals = split( "\t", $this_format );
					my %format_hash = ();
					while ( my ( $index, $value ) = each @format_fields ) {
						$format_hash{$value} = $format_vals[$index];
					}
					
					my $this_gtr      = $format_hash{'GTR'};
					my $this_gt       = $format_hash{'GT'};
					my $this_ad       = $format_hash{'AD'};
					my $this_ozyg       = $format_hash{'OZYG'};
					my $this_ogt       = $format_hash{'OGT'};
					my $this_zygosity = "NA";
					
					my $gtr2go = $this_gtr;
					
					if($mal_flag != 0){
						my $ovcf_key = $info_hash{'OVCF_KEY'};
						$this_ref = $fields_fil[2] unless ($mal_flag > 1);
						my %this_refs = ();
						my %this_alts = ();
						my @multi_vals = split("\n" , $mult_alle_vars{$ovcf_key} );
						$this_alts{0} = $this_ref;
						foreach(@multi_vals){
							my $this_order = (split("\t" , $_))[0];
							$this_alts{$this_order} = (split("\t" , $_))[4];
							$this_refs{$this_order} = (split("\t" , $_))[3];
						}
						
						if ( $this_ogt =~ m#^(\d)/(\d)# ) {
							my $gt1 = $1;
							my $gt2 = $2;
							my $new_ogtr = $this_alts{$1} . "/" . $this_alts{$2};
							if($gt1 == 0 && $gt2 == 0){
								$this_refs{$2} . "/" . $this_refs{$2};
							}elsif($gt1 == 0){
								$new_ogtr = $this_refs{$2} . "/" . $this_alts{$2};
							}elsif($gt2 == 0){
								$new_ogtr = $this_alts{$2} . "/" . $this_refs{$2};
							}else{
								$new_ogtr = $this_alts{$1} . "/" . $this_alts{$2};
							}
					    	$gtr2go  = $new_ogtr;
						
						}#$this_ogt
						
						
						
					} #$mal_flag
					
					if ( $this_gt =~ m#^(\d)/(\d)# ) {
						if ( $1 == $2 ) {
							$this_zygosity = ($1 == 0 && $2 == 0) ? "hom-ref" : "hom-alt";
						}elsif ( $1 == 0 || $2 == 0 ) {
							$this_zygosity = "ref-alt";
						}else{
							$this_zygosity = "alt-alt";
						}
						
					}elsif($this_gt eq "0" || $this_gt eq "1"){
						$this_zygosity = "hap";
					}else{
						$this_zygosity = "unknown";
					}
					$count_zyg{$this_zygosity} = (exists $count_zyg{$this_zygosity} ) ? ($count_zyg{$this_zygosity} + 1) : 1;
					my @ads = split( ",", $this_ad );
					my $ad_for_this_alt = ( $alt_count < ( scalar @ads ) ) ? $ads[$alt_count] : "NA";
					my $ad_for_this_ref = $ads[0];
					$format_hash{'AD'} = $ad_for_this_ref . "\t" . $ad_for_this_alt;
					my $this_alt_fract = (($ad_for_this_ref + $ad_for_this_alt) != 0) ? $ad_for_this_alt / ($ad_for_this_ref + $ad_for_this_alt) : 0;;
					my @tmp_out = ();
					my $this_GQ = -1;
					while ( my ( $index, $value ) = each @format_fields ) {
						push( @tmp_out, $format_hash{$value} );
						if($value eq "GQ"){
							$this_GQ = $format_hash{$value};
						}
					}
					push( @sample_values, join( "\t", @tmp_out ) );
					push( @sample_values, $this_zygosity );
					if(exists $dng_tracking{$smp_track_count}){
						my $this_dng_smp = shift @dng_smps_tmp;
						my $this_ppdnm = (exists $dng_info{$this_dng_smp}{$dng_line_key}{$this_dng_smp . ":PP_DNM"}) ? $dng_info{$this_dng_smp}{$dng_line_key}{$this_dng_smp . ":PP_DNM"} : "NA";
						#print "$this_dng_smp , $dng_line_key , this_ppdnm = $this_ppdnm , $this_alt_fract  , $this_GQ, $this_annovar_type\n";
						my $dng_conf = get_dng_tag($this_ppdnm , $this_alt_fract , $this_GQ , $this_annovar_type);
						foreach my $dcol (@dng_cols){
							my $d_hd = $this_dng_smp . ":" . $dcol;
							push(@sample_values , (exists $dng_info{$this_dng_smp}{$dng_line_key}{$d_hd}) ? $dng_info{$this_dng_smp}{$dng_line_key}{$d_hd} : "NA" );
						}
						push(@sample_values , $dng_conf);
						
					} #$smp_track_count
					
					
				}
				
				
				my @sample_valuesFx = ();
				foreach my $sv (@sample_values){
					$sv =~ s/^\s+//;
					$sv =~ s/\s+$//;
					$sv = ($sv eq '.') ? "NA" :  $sv;
					push(@sample_valuesFx, $sv);
				}
				
#				my @dng_values = ();
#				foreach my $sm (@dng_smps){
#					foreach my $dcol (@dng_cols){
#						my $d_hd = $sm . ":" . $dcol;
#						#print "$dng_line_key , $sm , $dcol , $d_hd \n ";
#						push(@dng_values , (exists $dng_info{$sm}{$dng_line_key}{$d_hd}) ? $dng_info{$sm}{$dng_line_key}{$d_hd} : "NA" );
#					}
#				}
#				
				

				my $line_to_go =
				    $fields[0] . "\t"
				  . $fields[1] . "\t"
				  . $fields[2] . "\t"
				  . $fields[3] . "\t"
				  . $fields[4]  . "\t"
#				  . $fields[13]  . "\t"
#				  . $fields[14];
				  
				  . join( "\t", @fields[ 8 .. ( scalar @fields ) - 1 ] ) ;
				  
				  print $fh_rout $line_to_go, "\t ", join( "\t", @infos ) . "\t" . join( "\t", @sample_valuesFx) ,   "\n";				  
		}
		

		$rfmt_stats{'var-counts'} = $fh_ann_count;
		$rfmt_stats{'multi_allele'} = $count_malleles;
		$rfmt_stats{'c2annovar_warn'} = $warn_c2annovar;
		foreach (keys %count_zyg){
			$rfmt_stats{$_} = $count_zyg{$_};
		}
	
		foreach (keys %count_alts){
			$rfmt_stats{$_} = $count_alts{$_};
		}	
	
		close $fh_rout;
		close $fh_vcf;
	
		return \%rfmt_stats;	
	
} #end of sub





###########################################

sub get_dng_tag{
	my ($sub_ppdnm , $sub_alt_fract , $sub_GQ , $sub_annovar_type) = @_;
	my $sub_conf = "NA";
	my @failed = ();
	my $conf_high = 1;
	if($sub_ppdnm ne "NA" ){
		if(lc($sub_annovar_type) eq "snp"){
			push(@failed , ($sub_ppdnm < 0.95) ? "PP_DNM<0.95": 1 );
			push(@failed , ($sub_GQ != 99) ? "GQ<99": 1 );
			push(@failed , ($sub_alt_fract < 0.3) ? "alt_frac<0.3": 1 );
			my $tag = join("," , @failed);
			$sub_conf = ($tag eq "1,1,1") ? "High" : "Low(" . $tag . ")";  
		}elsif(lc($sub_annovar_type) ne "unknown"){
			push(@failed , ($sub_GQ < 90) ? "GQ<90": 1 );
			push(@failed , ($sub_alt_fract < 0.3) ? "alt_frac<0.3": 1 );
			my $tag = join("," , @failed);
			$sub_conf = ($tag eq "1,1") ? "High" : "Low(" . $tag . ")";  
		}else{
			$sub_conf = "Low(unknown vartype)";
		}
		
		
	}
	return $sub_conf;	
}


#Old reformatting
sub reformat_vcf_dng{
	my $self = shift;
	my $vcffile = $self->{"name"};
	my $infos = shift;
	my $formats = shift;
	my $parsed_out = shift;
	my $denovo_gear = shift;
	my @info_filters   = ();
	my @format_filters = ();
	my $info_fls = $infos;
	
	my $vdir  = dirname($parsed_out);
	my $prefix = basename($vcffile);
	my $tmp_out = $vdir .  "/" . $prefix . "_parsed.tsv";
	my $parse_command =
		  "bcftools query -f \'\%CHROM\\t\%POS\\t\%REF\\t\%ALT\\t\%QUAL\\t%FILTER";
	
		foreach my $inf (split("," , $infos)) {
			push( @info_filters, "\\t\%INFO\/$inf" );
		}
        
		push( @format_filters, "\%TGT" ) unless ($denovo_gear == 1);
		foreach my $inf (split(",",$formats)) {
			push( @format_filters, "\\t\%$inf" );
		}
		#push( @format_filters, "\%GTR" );
		my $format_fls = join( ",", @format_filters );
		$parse_command =
		    $parse_command
		  . join( "", @info_filters )
		  . "\[\\t\{"
		  . join( "", @format_filters )
		  . "}\]\n\' $vcffile > $tmp_out";
		print STDOUT "Parse multi-sample vcf file: \n";
		print STDOUT "$parse_command\n";
		my $ex_st = system($parse_command);
		if ( $ex_st == 0 ) {
			print STDOUT "** parse-vcf is finished  with exit code: $?\n";
		}
		else {
			warn "parse-vcf is failed with exit code: $?\n";
			usage("Critical:parse-vcf is failed with exit code: $?");
		}
		
		#reformat2
#		my (
#		$org_vcf, $fil_vcf,   $rel_info, $smp_file,
#		$info_ls, $format_ls, $rfmt_out, $vhd
#	) = @_;
		
		my $org_vcf = $vcffile;
		my $fil_vcf = $tmp_out;
		my $rel_info = "rel";
		my $rfmt_out = $parsed_out;
		my @all_samples = $self->get_samples();
		#$info_ls, $format_ls,
		my @info_fields = $self->parse_INFO_fields();
		print "INFOs: " , join("," , @info_fields) , "\n";
		my @format_fields = ();
		@format_fields = ("GTR") unless ($denovo_gear == 1);
		push (@format_fields, $self-> parse_FORMAT_fields());

		#open( VCF,  $org_vcf )  or die "$!";
		my $fh_vcf = IO::File->new($org_vcf , q{<}) or die "$! $org_vcf\n";
		#open( FIL,  $fil_vcf )  or die "$!";
		my $fh_fil = IO::File->new($fil_vcf , q{<}) or die "$! $fil_vcf\n";
		#open( FAM,  $rel_info ) or warn "Pedigree info file is missing (OK!)\n";
		my $fh_fam = IO::File->new($rel_info , q{<}) or warn "Pedigree info file is missing (OK!) \n";
		
#		print join("," , @all_samples) , " in ref2 \n"; 
		my %fam        = ();
		my %fil        = ();
		my %filinfo = ();
		my %sample_ext = ();
		
		while ( my $line = <$fh_fam> ) {
			chomp($line);
			$line =~ s/\r//;
			my $id     = ( split( "\t", $line ) )[1];
			my $family = ( split( "\t", $line ) )[0];
			my $sex    = ( split( "\t", $line ) )[3];
			my $rel    = ( split( "\t", $line ) )[2];
			my $rel_to_go = substr( $rel, 0, 1 );
			$fam{$id} = $rel_to_go . "_" . $id;
		}
		#close $fh_fam;
	
	print  "Scanning for multi-allelic sites...\n";
	my %mult_alle_vars	= ();
	while ( my $line = <$fh_fil> ) {
		chomp($line);
		my @fields = split( "\t", $line );
		my @fields_tmp = @fields[ 6 .. ( scalar @fields ) - 1 ];
		my @infos = @fields_tmp[ 0 .. ( scalar @info_fields ) - 1 ];
		
		my %info_hash = ();
		while ( my ( $index, $value ) = each @info_fields ) {
			$info_hash{$value} = $infos[$index];
		}
		
		my $mal_flag = $info_hash{'MULTI_ALLELIC'};
		if($mal_flag != 0){
			my $ovcf_key = $info_hash{'OVCF_KEY'};
			#print "MALF: $mal_flag , $ovcf_key \n";
			$mult_alle_vars{$ovcf_key} = (exists $mult_alle_vars{$ovcf_key}) ? $mult_alle_vars{$ovcf_key} . "\n" . ($mal_flag . "\t" . $line)  : ($mal_flag . "\t" . $line);
		}
	}
	close $fh_fil;
	
	$fh_fil = IO::File->new($fil_vcf , q{<}) or die "$! $fil_vcf\n";
	#
	my $this_ref = ".";
	while ( my $line = <$fh_fil> ) {
		chomp($line);
		my @fields = split( "\t", $line );
		my $key =
		  $fields[0] . ":" . $fields[1] . ":" . $fields[2] . ":" . $fields[3];
		my $alt        = $fields[3];
		my @fields_tmp = @fields[ 6 .. ( scalar @fields ) - 1 ];
		my $field_tab  = join( "\t", @fields_tmp );

		my @sample_formats = ();
		my @infos = @fields_tmp[ 0 .. ( scalar @info_fields ) - 1 ];
		
		my %info_hash = ();
		while ( my ( $index, $value ) = each @info_fields ) {
			$info_hash{$value} = $infos[$index];
		}
		my $mal_flag = $info_hash{'MULTI_ALLELIC'};
		
		while ( $field_tab =~ m/\{(.+?)\}/g ) {
			#	print $key , "\t" , $1 , "\n";
			#	my $infos = $1;
			push( @sample_formats, $1 );
		}

		my $alt_count = 1;
		foreach my $this_alt ( split( ",", $alt ) ) {
			my $new_key =
			    $fields[0] . ":"
			  . $fields[1] . ":"
			  . $fields[2] . ":"
			  . $this_alt;
			my @sample_values = ();
			
			foreach my $this_format (@sample_formats) {
				$this_format =~ s/^\s+//;    #remove leading spaces
				$this_format =~ s/\s+$//;    #remove trailing spaces
				my @format_vals = split( "\t", $this_format );
				my %format_hash = ();
				while ( my ( $index, $value ) = each @format_fields ) {
					$format_hash{$value} = $format_vals[$index];
				}



#			my ($this_gtr, $this_gt, $this_ab , $this_ad , $this_dp , $this_gq , $this_pl) = split("\t" , $this_format);
				my $this_gtr      = $format_hash{'GTR'};
				my $this_gt       = $format_hash{'GT'};
				my $this_ad       = $format_hash{'AD'};
				my $this_ozyg       = $format_hash{'OZYG'};
				my $this_ogt       = $format_hash{'OGT'};
				my $this_zygosity = "NA";
				
				#gtr2go updated July 2017
				my $gtr2go = $this_gtr;
				
				if($mal_flag != 0){
					my $ovcf_key = $info_hash{'OVCF_KEY'};
					$this_ref = $fields[2] unless ($mal_flag > 1);
					my %this_refs = ();
					my %this_alts = ();
					my @multi_vals = split("\n" , $mult_alle_vars{$ovcf_key} );
					$this_alts{0} = $this_ref;
					foreach(@multi_vals){
						my $this_order = (split("\t" , $_))[0];
						#1       chr1    29845252        C       CT      1871.19 PASS 
						#2       chr1    29845252        CTCT    C       1871.19 PASS    1       0.5 
						$this_alts{$this_order} = (split("\t" , $_))[4];
						$this_refs{$this_order} = (split("\t" , $_))[3];
					}
					
					#OGT
					if ( $this_ogt =~ m#^(\d)/(\d)# ) {
						my $gt1 = $1;
						my $gt2 = $2;
						my $new_ogtr = $this_alts{$1} . "/" . $this_alts{$2};
						if($gt1 == 0 && $gt2 == 0){
							$this_refs{$2} . "/" . $this_refs{$2};
						}elsif($gt1 == 0){
							$new_ogtr = $this_refs{$2} . "/" . $this_alts{$2};
						}elsif($gt2 == 0){
							$new_ogtr = $this_alts{$2} . "/" . $this_refs{$2};
						}else{
							$new_ogtr = $this_alts{$1} . "/" . $this_alts{$2};
						}
					    $gtr2go  = $new_ogtr;
					}
					
				}
				
				#print "GTR2go: " , $new_key , "\t" , $mal_flag , "\t" , $gtr2go ,  "\t" , $this_gt , "\t" , $this_gtr ,  "\n" if($mal_flag != 0);
				
				
				#my $gtr2go = ($this_ozyg eq "alt-alt") ? (split('/', $this_gtr))[1] : $this_gtr;

		#			warn "GT, GTR = $this_gt , $this_gtr\n";
				if ( $this_gt =~ m#^(\d)/(\d)# ) {
					if ( $1 == $2 ) {
						#$this_zygosity = "hom";
						$this_zygosity = ($1 == 0 && $2 == 0) ? "hom-ref" : "hom-alt";
					}
					elsif ( $1 == 0 || $2 == 0 ) {
						$this_zygosity = "ref-alt";
					}
					else {
						$this_zygosity = "alt-alt";
					}
				}elsif($this_gt eq "0" || $this_gt eq "1"){
					$this_zygosity = "hap";
				}else {
					$this_zygosity = "unknown";
				}
				
				

				my @ads = split( ",", $this_ad );
				my $ad_for_this_alt =
				  ( $alt_count < ( scalar @ads ) ) ? $ads[$alt_count] : "NA";
				my $ad_for_this_ref = $ads[0];
				$format_hash{'AD'} = $ad_for_this_ref . "\t" . $ad_for_this_alt;
				$format_hash{'GTR'} = $gtr2go;
				my @tmp_out = ();
				while ( my ( $index, $value ) = each @format_fields ) {
					push( @tmp_out, $format_hash{$value} );
				}
				push( @sample_values, join( "\t", @tmp_out ) );
				push( @sample_values, $this_zygosity );

			}
			$alt_count++;
			$fil{$new_key} =
			  join( "\t", @infos ) . "\t" . join( "\t", @sample_values );

		}

	}
	close $fh_fil;

	foreach my $s (@all_samples) {
		$sample_ext{$s} = ( exists $fam{$s} ) ? $fam{$s} : $s;
	}
	
		#open( ROUT, ">$rfmt_out" ) or die "error $rfmt_out $!";
		my $fh_rout = IO::File->new($rfmt_out , q{>}) or die "$! $rfmt_out\n";		
	while ( my $line = <$fh_vcf> ) {
		chomp($line);
		if ( $line =~ m/^#/ ) {
			###############
			if ( $line =~ m/^#CHROM/ ) {
				my $samples = ( split( "\tFORMAT\t", $line ) )[1];
				$samples =~ s/\s+/,/g;
				@all_samples = split( ",", $samples );
				foreach my $k (@all_samples) {
					my $k_val = ( exists $fam{$k} ) ? $fam{$k} : $k;
					my $new_id = $k_val;
					$line =~ s/$k/$new_id/;
				}
				print $fh_rout $line, "\t";
				my @new_head = ();
				push( @new_head, join( "\t", @info_fields ) );
				foreach my $k (@all_samples) {
					my $new_h = join( "\t", @format_fields ) . "\tZygosity";
					$new_h =~ s/AD/AD_REF\tAD_ALT/;
					foreach my $nh ( split( "\t", $new_h ) ) {
						my $format_prefix = ( exists $fam{$k} ) ? $fam{$k} : $k;
						my $h = $format_prefix . ":" . $nh;
						push( @new_head, $h );
					}
				}
				print $fh_rout join( "\t", @new_head ), "\n";

			}

			###############
		}
		else {
			my @fields = split( "\t", $line );
			my $alt = $fields[4];
			foreach my $this_alt ( split( ",", $alt ) ) {
				my $new_key =
				    $fields[0] . ":"
				  . $fields[1] . ":"
				  . $fields[3] . ":"
				  . $this_alt;
				my $line_to_go =
				    $fields[0] . "\t"
				  . $fields[1] . "\t"
				  . $fields[2] . "\t"
				  . $fields[3] . "\t"
				  . $this_alt . "\t"
				  . join( "\t", @fields[ 5 .. ( scalar @fields ) - 1 ] );

				#		print $fil{$new_key} , "\n";
#				print $line_to_go , "\n";
#				print $fil{$new_key} , "\n";
				print $fh_rout $line_to_go, "\t ", $fil{$new_key}, "\n";
			}

		}

	}
		

	close $fh_rout;
	close $fh_vcf;
	return 0;
	
}



sub parse_INFO_fields{
	my $self = shift;
	my $vfile = $self->{"name"};
	#open( TMPV, $vfile );
	my $fh_tmpv = IO::Zlib->new($vfile, "rb") or  usage( "$! $vfile");
	if(not defined $fh_tmpv){
		usage( "$! error opening $vfile\n");
	}
	
#	my $fh_tmpv = IO::File->new($vfile , q{<}) or usage( "$! $vfile");
	my @infos = ();
	while ( my $line = <$fh_tmpv> ) {
		chomp($line);
		if ( $line =~ m/^#CHROM/ ) { last; }
		else {

# ##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
			if ( $line =~ m/^##INFO=.+ID=(.+?),.+/ ) {
				push( @infos, $1 );
			}

		}
	}
	close $fh_tmpv;
	return @infos;
}

sub parse_FORMAT_fields {
	my $self = shift;
	my $vfile = $self->{"name"};
	#open( TMPV, $vfile );
	#my $fh_tmpv = IO::File->new($vfile , q{<}) or usage( "$! $vfile");
	my $fh_tmpv = IO::Zlib->new($vfile, "rb") or  usage( "$! $vfile");
	if(not defined $fh_tmpv){
		usage( "$! error opening $vfile\n");
	}
	my @infos = ();
	while ( my $line = <$fh_tmpv> ) {
		chomp($line);
		if ( $line =~ m/^#CHROM/ ) { last; }
		else {

# ##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
			if ( $line =~ m/^##FORMAT=.+ID=(.+?),.+/ ) {
				push( @infos, $1 );
			}

		}
	}
	close $fh_tmpv;
	return @infos;
}

sub chrstat{
	my $self = shift;
	my $vfile = $self->{"name"};
	# ##contig=<ID=chr1,assembly=b37,length=249250621>
	#open( TMPV, $vfile );
	my $fh_tmpv = IO::File->new($vfile , q{<}) or usage( "$! $vfile");
	my $ischr = 0;
	while ( my $line = <$fh_tmpv> ) {
		chomp($line);
		if ( $line =~ m/^#CHROM/ ) { last; }
		else {

			if ( $line =~ m/^##contig=.+ID=(.+?),.+/ ) {
				$ischr = ($1 =~ m/^chr/) ? 1 : 0;
			}

		}
	}
	close $fh_tmpv;
	$self->{"chr"} = $ischr;
	return $ischr;
}


sub reformat_vcf_old{
	my $self = shift;
	my $vcffile = $self->{"name"};
	my $infos = shift;
	my $formats = shift;
	my $parsed_out = shift;
	my $denovo_gear = shift;
	my @info_filters   = ();
	my @format_filters = ();
	my $info_fls = $infos;
	
	my $vdir  = dirname($parsed_out);
	my $prefix = basename($vcffile);
	my $tmp_out = $vdir .  "/" . $prefix . "_parsed.tsv";
	my $parse_command =
		  "bcftools query -f \'\%CHROM\\t\%POS\\t\%REF\\t\%ALT\\t\%QUAL\\t%FILTER";
	
		foreach my $inf (split("," , $infos)) {
			push( @info_filters, "\\t\%INFO\/$inf" );
		}
        
		push( @format_filters, "\%TGT" ) unless ($denovo_gear == 1);
		foreach my $inf (split(",",$formats)) {
			push( @format_filters, "\\t\%$inf" );
		}
		#push( @format_filters, "\%GTR" );
		my $format_fls = join( ",", @format_filters );
		$parse_command =
		    $parse_command
		  . join( "", @info_filters )
		  . "\[\\t\{"
		  . join( "", @format_filters )
		  . "}\]\n\' $vcffile > $tmp_out";
		print STDOUT "Parse multi-sample vcf file: \n";
		print STDOUT "$parse_command\n";
		my $ex_st = system($parse_command);
		if ( $ex_st == 0 ) {
			print STDOUT "** parse-vcf is finished  with exit code: $?\n";
		}
		else {
			warn "parse-vcf is failed with exit code: $?\n";
			usage("Critical:parse-vcf is failed with exit code: $?");
		}
		
		#reformat2
#		my (
#		$org_vcf, $fil_vcf,   $rel_info, $smp_file,
#		$info_ls, $format_ls, $rfmt_out, $vhd
#	) = @_;
	
		
		
		
		my $org_vcf = $vcffile;
		my $fil_vcf = $tmp_out;
		my $rel_info = "rel";
		my $rfmt_out = $parsed_out;
		my @all_samples = $self->get_samples();
		#$info_ls, $format_ls,
		my @info_fields = $self->parse_INFO_fields();
		print join("," , @info_fields) , "\n";
		my @format_fields = ();
		@format_fields = ("GTR") unless ($denovo_gear == 1);
		push (@format_fields, $self-> parse_FORMAT_fields());

		#open( VCF,  $org_vcf )  or die "$!";
		my $fh_vcf = IO::File->new($org_vcf , q{<}) or die "$! $org_vcf\n";
		#open( FIL,  $fil_vcf )  or die "$!";
		my $fh_fil = IO::File->new($fil_vcf , q{<}) or die "$! $fil_vcf\n";
		#open( FAM,  $rel_info ) or warn "Pedigree info file is missing (OK!)\n";
		my $fh_fam = IO::File->new($rel_info , q{<}) or warn "Pedigree info file is missing (OK!) \n";
		
#		print join("," , @all_samples) , " in ref2 \n"; 
		my %fam        = ();
		my %fil        = ();
		my %filinfo = ();
		my %sample_ext = ();
		
		while ( my $line = <$fh_fam> ) {
			chomp($line);
			$line =~ s/\r//;
			my $id     = ( split( "\t", $line ) )[1];
			my $family = ( split( "\t", $line ) )[0];
			my $sex    = ( split( "\t", $line ) )[3];
			my $rel    = ( split( "\t", $line ) )[2];
			my $rel_to_go = substr( $rel, 0, 1 );
			$fam{$id} = $rel_to_go . "_" . $id;
		}
		close FAM;
		
	while ( my $line = <$fh_fil> ) {
		chomp($line);
		my @fields = split( "\t", $line );
		my $key =
		  $fields[0] . ":" . $fields[1] . ":" . $fields[2] . ":" . $fields[3];
		my $alt        = $fields[3];
		my @fields_tmp = @fields[ 6 .. ( scalar @fields ) - 1 ];
		my $field_tab  = join( "\t", @fields_tmp );

		my @sample_formats = ();
		my @infos = @fields_tmp[ 0 .. ( scalar @info_fields ) - 1 ];
		while ( $field_tab =~ m/\{(.+?)\}/g ) {

			#	print $key , "\t" , $1 , "\n";
			#	my $infos = $1;
			push( @sample_formats, $1 );
		}

		my $alt_count = 1;
		foreach my $this_alt ( split( ",", $alt ) ) {
			my $new_key =
			    $fields[0] . ":"
			  . $fields[1] . ":"
			  . $fields[2] . ":"
			  . $this_alt;
			my @sample_values = ();
			foreach my $this_format (@sample_formats) {
				$this_format =~ s/^\s+//;    #remove leading spaces
				$this_format =~ s/\s+$//;    #remove trailing spaces
				my @format_vals = split( "\t", $this_format );
				my %format_hash = ();
				while ( my ( $index, $value ) = each @format_fields ) {
					$format_hash{$value} = $format_vals[$index];
				}

#			my ($this_gtr, $this_gt, $this_ab , $this_ad , $this_dp , $this_gq , $this_pl) = split("\t" , $this_format);
				my $this_gtr      = $format_hash{'GTR'};
				my $this_gt       = $format_hash{'GT'};
				my $this_ad       = $format_hash{'AD'};
				my $this_ozyg       = $format_hash{'OZYG'};
				my $this_zygosity = "NA";
				my $gtr2go = ($this_ozyg eq "alt-alt") ? (split('/', $this_gtr))[1] : $this_gtr;

		#			warn "GT, GTR = $this_gt , $this_gtr\n";
				if ( $this_gt =~ m#^(\d)/(\d)# ) {
					if ( $1 == $2 ) {
						#$this_zygosity = "hom";
						$this_zygosity = ($1 == 0 && $2 == 0) ? "hom-ref" : "hom-alt";
					}
					elsif ( $1 == 0 || $2 == 0 ) {
						$this_zygosity = "ref-alt";
					}
					else {
						$this_zygosity = "alt-alt";
					}
				}elsif($this_gt eq "0" || $this_gt eq "1"){
					$this_zygosity = "hap";
				}else {
					$this_zygosity = "unknown";
				}

				my @ads = split( ",", $this_ad );
				my $ad_for_this_alt =
				  ( $alt_count < ( scalar @ads ) ) ? $ads[$alt_count] : "NA";
				my $ad_for_this_ref = $ads[0];
				$format_hash{'AD'} = $ad_for_this_ref . "\t" . $ad_for_this_alt;
				$format_hash{'GTR'} = $gtr2go;
				my @tmp_out = ();
				while ( my ( $index, $value ) = each @format_fields ) {
					push( @tmp_out, $format_hash{$value} );
				}
				push( @sample_values, join( "\t", @tmp_out ) );
				push( @sample_values, $this_zygosity );

			}
			$alt_count++;
			$fil{$new_key} =
			  join( "\t", @infos ) . "\t" . join( "\t", @sample_values );

		}

	}
	close $fh_fil;

	foreach my $s (@all_samples) {
		$sample_ext{$s} = ( exists $fam{$s} ) ? $fam{$s} : $s;
	}
	
		#open( ROUT, ">$rfmt_out" ) or die "error $rfmt_out $!";
		my $fh_rout = IO::File->new($rfmt_out , q{>}) or die "$! $rfmt_out\n";		
	while ( my $line = <$fh_vcf> ) {
		chomp($line);
		if ( $line =~ m/^#/ ) {
			###############
			if ( $line =~ m/^#CHROM/ ) {
				my $samples = ( split( "\tFORMAT\t", $line ) )[1];
				$samples =~ s/\s+/,/g;
				@all_samples = split( ",", $samples );
				foreach my $k (@all_samples) {
					my $k_val = ( exists $fam{$k} ) ? $fam{$k} : $k;
					my $new_id = $k_val;
					$line =~ s/$k/$new_id/;
				}
				print $fh_rout $line, "\t";
				my @new_head = ();
				push( @new_head, join( "\t", @info_fields ) );
				foreach my $k (@all_samples) {
					my $new_h = join( "\t", @format_fields ) . "\tZygosity";
					$new_h =~ s/AD/AD_REF\tAD_ALT/;
					foreach my $nh ( split( "\t", $new_h ) ) {
						my $format_prefix = ( exists $fam{$k} ) ? $fam{$k} : $k;
						my $h = $format_prefix . ":" . $nh;
						push( @new_head, $h );
					}
				}
				print $fh_rout join( "\t", @new_head ), "\n";

			}

			###############
		}
		else {
			my @fields = split( "\t", $line );
			my $alt = $fields[4];
			foreach my $this_alt ( split( ",", $alt ) ) {
				my $new_key =
				    $fields[0] . ":"
				  . $fields[1] . ":"
				  . $fields[3] . ":"
				  . $this_alt;
				my $line_to_go =
				    $fields[0] . "\t"
				  . $fields[1] . "\t"
				  . $fields[2] . "\t"
				  . $fields[3] . "\t"
				  . $this_alt . "\t"
				  . join( "\t", @fields[ 5 .. ( scalar @fields ) - 1 ] );

				#		print $fil{$new_key} , "\n";
				print $fh_rout $line_to_go, "\t ", $fil{$new_key}, "\n";
			}

		}

	}
		

	close $fh_rout;
	close $fh_vcf;
	return 0;
	
}


#
#sub reformat_vcf2 {
#	my (
#		$org_vcf, $fil_vcf,   $rel_info, $smp_file,
#		$info_ls, $format_ls, $rfmt_out, $vhd
#	) = @_;
#	my @info_fields   = split( ",", $info_ls );
#	my @format_fields = split( ",", $format_ls );
#
#	print STDOUT "REFORMAT(sub)---\n";
#	open( VCF,  $org_vcf )  or die "$!";
#	open( FIL,  $fil_vcf )  or die "$!";
#	open( FAM,  $rel_info ) or warn "$! rel_info is missing\n";
#	open( SAMP, $smp_file ) or die "$! sample file ($smp_file)";
#	my @all_samples = <SAMP>;
#	close SAMP;
#
#	my %fam        = ();
#	my %fil        = ();
#	my %sample_ext = ();
#
#	while ( my $line = <FAM> ) {
#		chomp($line);
#		$line =~ s/\r//;
#		my $id     = ( split( "\t", $line ) )[1];
#		my $family = ( split( "\t", $line ) )[0];
#		my $sex    = ( split( "\t", $line ) )[3];
#		my $rel    = ( split( "\t", $line ) )[2];
#		my $rel_to_go = substr( $rel, 0, 1 );
#		$fam{$id} = $rel_to_go . "_" . $id;
#	}
#	close FAM;
#
#	while ( my $line = <FIL> ) {
#		chomp($line);
#		my @fields = split( "\t", $line );
#		my $key =
#		  $fields[0] . ":" . $fields[1] . ":" . $fields[2] . ":" . $fields[3];
#		my $alt        = $fields[3];
#		my @fields_tmp = @fields[ 6 .. ( scalar @fields ) - 1 ];
#		my $field_tab  = join( "\t", @fields_tmp );
#
#		my @sample_formats = ();
#		my @infos = @fields_tmp[ 0 .. ( scalar @info_fields ) - 1 ];
#		while ( $field_tab =~ m/\{(.+?)\}/g ) {
#			push( @sample_formats, $1 );
#		}
#
#		my $alt_count = 1;
#		foreach my $this_alt ( split( ",", $alt ) ) {
#			my $new_key =
#			    $fields[0] . ":"
#			  . $fields[1] . ":"
#			  . $fields[2] . ":"
#			  . $this_alt;
#			my @sample_values = ();
#			foreach my $this_format (@sample_formats) {
#				$this_format =~ s/^\s+//;    #remove leading spaces
#				$this_format =~ s/\s+$//;    #remove trailing spaces
#				my @format_vals = split( "\t", $this_format );
#				my %format_hash = ();
#				while ( my ( $index, $value ) = each @format_fields ) {
#					$format_hash{$value} = $format_vals[$index];
#				}
#
##			my ($this_gtr, $this_gt, $this_ab , $this_ad , $this_dp , $this_gq , $this_pl) = split("\t" , $this_format);
#				my $this_gtr      = $format_hash{'GTR'};
#				my $this_gt       = $format_hash{'GT'};
#				my $this_ad       = $format_hash{'AD'};
#				my $this_zygosity = "NA";
#
#				#	warn "GTR = $this_gt\n";
#				if ( $this_gtr =~ m#^(\d)/(\d)# ) {
#					if ( $1 == $2 ) {
#						$this_zygosity = "hom";
#					}
#					elsif ( $1 == 0 || $2 == 0 ) {
#						$this_zygosity = "ref-alt";
#					}
#					else {
#						$this_zygosity = "alt-alt";
#					}
#				}
#				else {
#					$this_zygosity = "unknown";
#				}
#
#				my @ads = split( ",", $this_ad );
#				my $ad_for_this_alt =
#				  ( $alt_count < ( scalar @ads ) ) ? $ads[$alt_count] : "NA";
#				my $ad_for_this_ref = $ads[0];
#				$format_hash{'AD'} = $ad_for_this_ref . "\t" . $ad_for_this_alt;
#				my @tmp_out = ();
#				while ( my ( $index, $value ) = each @format_fields ) {
#					push( @tmp_out, $format_hash{$value} );
#				}
#				push( @sample_values, join( "\t", @tmp_out ) );
#				push( @sample_values, $this_zygosity );
#
#			}
#			$alt_count++;
#			$fil{$new_key} =
#			  join( "\t", @infos ) . "\t" . join( "\t", @sample_values );
#
#		}
#
#	}
#	close FIL;
#
#	foreach my $s (@all_samples) {
#		$sample_ext{$s} = ( exists $fam{$s} ) ? $fam{$s} : $s;
#	}
#
#	open( ROUT, ">$rfmt_out" ) or die "error $rfmt_out $!";
#	open( HDR,  $vhd )         or die "error $vhd $!";
#	my @all_hdr = <HDR>;
#	close HDR;
#	chomp(@all_hdr);
#	foreach (@all_hdr) {
#		print ROUT $_, "\n" unless ( $_ =~ m/^#CHROM.+/ );
#	}
#
#	while ( my $line = <VCF> ) {
#		chomp($line);
#		if ( $line =~ m/^#/ ) {
#			###############
#			if ( $line =~ m/^#CHROM/ ) {
#				my $samples = ( split( "\tFORMAT\t", $line ) )[1];
#				$samples =~ s/\s+/,/g;
#				@all_samples = split( ",", $samples );
#				foreach my $k (@all_samples) {
#					my $k_val = ( exists $fam{$k} ) ? $fam{$k} : $k;
#					my $new_id = $k_val;
#					$line =~ s/$k/$new_id/;
#				}
#				print ROUT $line, "\t";
#				my @new_head = ();
#				push( @new_head, join( "\t", @info_fields ) );
#				foreach my $k (@all_samples) {
#					my $new_h = join( "\t", @format_fields ) . "\tZygosity";
#					$new_h =~ s/AD/AD_REF\tAD_ALT/;
#					foreach my $nh ( split( "\t", $new_h ) ) {
#						my $format_prefix = ( exists $fam{$k} ) ? $fam{$k} : $k;
#						my $h = $format_prefix . ":" . $nh;
#						push( @new_head, $h );
#					}
#				}
#				print ROUT join( "\t", @new_head ), "\n";
#
#			}
#
#			###############
#		}
#		else {
#			my @fields = split( "\t", $line );
#			my $alt = $fields[4];
#			foreach my $this_alt ( split( ",", $alt ) ) {
#				my $new_key =
#				    $fields[0] . ":"
#				  . $fields[1] . ":"
#				  . $fields[3] . ":"
#				  . $this_alt;
#				my $line_to_go =
#				    $fields[0] . "\t"
#				  . $fields[1] . "\t"
#				  . $fields[2] . "\t"
#				  . $fields[3] . "\t"
#				  . $this_alt . "\t"
#				  . join( "\t", @fields[ 5 .. ( scalar @fields ) - 1 ] );
#
#				#		print $fil{$new_key} , "\n";
#				print ROUT $line_to_go, "\t ", $fil{$new_key}, "\n";
#			}
#
#		}
#
#	}
#	close ROUT;
#	close VCF;
#	return 0;
#}

#$line =~ m/^##FORMAT=.+ID=(.+?),.+/
# ##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
sub add_orig_keys{
	#warn "***********UPDATED:************** \n"; 
	my $self = shift;
	my $vcffile = $self->{"name"};
	my $vcfheader = shift;
	my $outvcf = shift;
	my @all_samples = $self->get_samples();
	# ##INFO=<ID=OLD_MULTIALLELIC,Number=1,Type=String,Description="Original chr:pos:ref:alt encoding">
	# ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
	#open(OUTV, ">$outvcf") or die "$!";
	my $fh_outv = IO::File->new($outvcf , q{>}) or  usage("$! $outvcf\n");
	#open(HDR, $vcfheader ) or die "$vcfheader";
	my $fh_hdr = IO::File->new($vcfheader , q{<}) or usage("$! $vcfheader\n");
	
	my $eff_number_of_vars = 0;
	while(my $line = <$fh_hdr>){
		chomp($line);
		#check proper AD
		if($line =~ m/^##FORMAT=.+ID=AD,Number=(.).+/){
			my $ad_number = $1;
			if($ad_number eq '.'){
				warn "\n" , $line , "\n";
				warn "Warning! AD = $ad_number ,   \'AD NUMBER\' not correctly tagged (expected \'R\' or \'G\' instead of \'\.\') \n\n"; 
			}
		}
		
		if($line =~ m/^#CHROM/){
			print $fh_outv "##INFO=<ID=OVCF_KEY,Number=1,Type=String,Description=\"Original var key chr:pos:ref:alt encoding, before Normalization\">\n";
			print $fh_outv "##INFO=<ID=MULTI_ALLELIC,Number=A,Type=Integer,Description=\"Flags or index for multialleleic ,before decomposing\">\n";
			print $fh_outv "##FORMAT=<ID=OGT,Number=1,Type=String,Description=\"Old Genotype, before decomposing\">\n";
			print $fh_outv "##FORMAT=<ID=OZYG,Number=1,Type=String,Description=\"Zygosity, before decomposing\">\n";
			#print $fh_outv "##INFO=<ID=M_ORDER,Number=A,Type=String,Description=\"multi-alleleic order; before Normalization\">\n";
		}
		print $fh_outv $line , "\n";
	}
	close $fh_hdr;
	
	my $fh_vcf = undef;
	if ($vcffile =~ /\.gz$/) {
			warn "gzipped input file\n";
			open($fh_vcf, "gunzip -c $vcffile |") || die "can’t open pipe to $vcffile";
			#$fh_vcf = IO::File->new("gunzip -c $vcffile |" , q{<}) or usage("$! cannnot load input vcf $vcffile\n");
		}else{
			warn "$vcffile normal input file\n";
			#open(VCF, $vcffile) or die "$! $vcffile\n";
			$fh_vcf = IO::File->new($vcffile , q{<}) or usage("$! cannnot load input vcf $vcffile\n");
		}	
	#open(VCF, $vcffile ) or die "$! $vcffile\n";
	while(my $line = <$fh_vcf>){
		chomp($line);
		if($line =~ m/^#/){next;}else{
			my @fields = split("\t", $line);
			my $vcf_key =  $fields[0] . ":" . $fields[1] . ":" . $fields[3] . ":" . $fields[4];
			my $m_alleleic = 0;
			my $count_comma = 0;
			if($fields[4] =~ m/,/){
				my $field_4 = $fields[4];
				$count_comma = $field_4 =~ tr/,//;
				my @tmp_morder = (1 .. ($count_comma+1));
				#print $vcf_key , "\t" , $count_comma , "\t" , join("," , @tmp_morder) , "\n" ;
				$m_alleleic = join("," , @tmp_morder);
			}
			$eff_number_of_vars = $eff_number_of_vars + 1 + $count_comma;
			$fields[7] =~ s/\.//  if ($fields[7] eq "\.");
			my $new_info_field =  $fields[7] . ";OVCF_KEY=" . $vcf_key . ";MULTI_ALLELIC=" . $m_alleleic ;
			my @format_specs =  split(":" , $fields[8]);
			#my $new_format_spec = $fields[8] . ":OGT:OZYG";
			my $new_format_spec = "OGT:OZYG:" . $fields[8]; 
			my $sm_index = 9;
			my %sample_value = ();
			my @new_sample_value = ();
			foreach my $sm (@all_samples){
				my @this_samaple_values = split(":" , $fields[$sm_index++]);
				while( my( $index, $value ) = each @format_specs ){
					$sample_value{$sm}{$value} = $this_samaple_values[$index];
					if($value eq "GT"){
						my $this_gt = $this_samaple_values[$index];
						my $this_zyg = "unknown";
						if ( $this_gt =~ m#^(\d)/(\d)# ) {
						 if ( $1 == $2 ) {
						 	$this_zyg = ($1 == 0 && $2 == 0) ? "hom-ref" : "hom-alt";
							#$this_zyg = "hom";
						 }
						 elsif ( $1 == 0 || $2 == 0 ) {
							$this_zyg = "ref-alt";
						 }
						 else {
							$this_zyg = "alt-alt";
						 }
					 }elsif($this_gt eq "0" || $this_gt eq "1"){
						$this_zyg = "hap";
					 }else {
						 $this_zyg = "unknown";
					 }
						$sample_value{$sm}{'zyg'} = $this_zyg;
					}
				}
				#push(@new_sample_value , join(":" , @this_samaple_values) . ":"  . $sample_value{$sm}{'GT'} . ":" . $sample_value{$sm}{'zyg'});
				push(@new_sample_value ,  $sample_value{$sm}{'GT'} . ":" . $sample_value{$sm}{'zyg'} . ":" .  join(":" , @this_samaple_values) );
			}
			print $fh_outv $fields[0] . "\t" .  $fields[1] . "\t" .  $fields[2] . "\t" . $fields[3] . "\t" . $fields[4] . "\t" . $fields[5] . "\t" . $fields[6] . "\t" . $new_info_field . "\t" . $new_format_spec . "\t" . join("\t" , @new_sample_value) , "\n";
		}
	}
	close $fh_vcf;
	close $fh_outv;
	return $eff_number_of_vars;
	
}

sub add_multiallele_order{
	my $self = shift;
	my $vcffile = $self->{"name"};
	my $vcfheader = shift;
	my $outvcf = shift;
	my @all_samples = $self->get_samples();
	# ##INFO=<ID=OLD_MULTIALLELIC,Number=1,Type=String,Description="Original chr:pos:ref:alt encoding">
	# ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#	open(OUTV, ">$outvcf") or die "$!";
	my $fh_outv = IO::File->new($outvcf , q{>}) or usage("$! $outvcf\n");
#	open(HDR, $vcfheader ) or die "$vcfheader";
	my $fh_hdr = IO::File->new($vcfheader , q{<}) or usage("$! $vcfheader\n");
	while(my $line = <$fh_hdr>){
		chomp($line);
		
		if($line =~ m/^#CHROM/){
			print $fh_outv "##INFO=<ID=M_ORDER,Number=1,Type=String,Description=\"multi-alleleic order; before Normalization\">\n";
		}
		print $fh_outv $line , "\n";
	}
	close $fh_hdr;
	#open(VCF, $vcffile ) or die "$! $vcffile\n";
	my $fh_vcf = IO::File->new($vcffile , q{<}) or usage("$! $vcffile\n");
	
	while(my $line = <$fh_vcf>){
		chomp($line);
		if($line =~ m/^#/){next;}else{
			    my @fields = split("\t", $line);
				my $info = $fields[7];
				my $this_alt_allele = $fields[4];
				my $is_multi = 0;
				my $ovcf_key = "0";
				my $m_order = 0;
                if($info =~ /.*MULTI_ALLELIC=(.+)/){
					$is_multi = $1;
				}
                if($info =~ /.*OVCF_KEY=(.+?);/){
					$ovcf_key = $1;
					if ($ovcf_key eq "0"){
						die "Cannot parse vcf key , $ovcf_key\n";
					}
				}
				if($is_multi == 1){
					#multi-alleleic , OVCF_KEY=2:175614908:CAAA:CA,C;MULTI_ALLELELIC=1
					my @ovcfs = split(":", $ovcf_key);
					my @this_multis = split("," , $ovcfs[3]);
					while( my( $index, $value ) = each @this_multis ) {
						if($value eq $this_alt_allele){
							$m_order = $index + 1;
						}
					}
				}
				#add order
				$fields[7] =  $fields[7] . ";M_ORDER=" . $m_order ;
				print $fh_outv join("\t" , @fields) , "\n";
 		}
	}
	close $fh_vcf;
	close $fh_outv;
}

sub fix_multiallele_var{
	my $self = shift;
	my $vcffile = $self->{"name"};
	my $smps = shift;
	my $outvcf = shift;
	my @samples = split("," , $smps); 
#	open(VCF, $vcffile ) or die "$! $vcffile\n";
	my $fh_vcf = IO::File->new($vcffile , q{<}) or usage("$! $vcffile\n");
#	open(OUTV, ">$outvcf") or die "$!";
	my $fh_outv = IO::File->new($outvcf , q{>}) or usage("$! $outvcf\n");
	my $header = <$fh_vcf>;
	my @hds = split("\t" , $header);
	chomp($header);
	my %mult_alle_vars = ();
	
	while(my $line = <$fh_vcf>){
		chomp($line);
		my @fields = split("\t", $line);
		my %tmp_row = ();
		while( my( $index, $value ) = each @hds ) {
			 $tmp_row{$value} = $fields[$index];
		}
		
		my $mal_flag = $tmp_row{'MULTI_ALLELIC'};
		if($mal_flag != 0){
			#print "$mal_flag\n";
			my $ovcf_key = $tmp_row{'OVCF_KEY'};
			my $m_order = $tmp_row{'MULTI_ALLELIC'};
#			print "$line\n";
#			print "$m_order\n";
			$mult_alle_vars{$ovcf_key} = (exists $mult_alle_vars{$ovcf_key}) ? $mult_alle_vars{$ovcf_key} . "\n" . ($m_order . "\t" . $line)  : ($m_order . "\t" . $line);
		}else{
			#regular
		}
	}
	close $fh_vcf;
#	open($fh_vcf, $vcffile ) or die "$! $vcffile\n";
	$fh_vcf = IO::File->new($vcffile , q{<}) or usage("$! $vcffile\n");
	
	my $header = <$fh_vcf>;
	my @hds = split("\t" , $header);
	chomp($header);
	
	#print header
	print $fh_outv $header , "\n";
	my $this_ref = ".";
	my $eff_number_of_vars = 0;
	while(my $line = <$fh_vcf>){
		$eff_number_of_vars++;
		chomp($line);
		my @fields = split("\t", $line);
		my %tmp_row = ();
		while( my( $index, $value ) = each @hds ) {
			 $tmp_row{$value} = $fields[$index];
		}
		my $mal_flag = $tmp_row{'MULTI_ALLELIC'};
		if($mal_flag != 0){
			my $ovcf_key = $tmp_row{'OVCF_KEY'};
			my $m_order = $tmp_row{'MULTI_ALLELIC'};
			$this_ref = $tmp_row{'REF'} unless ($m_order > 1);
			my %this_refs = ();
			my %this_alts = ();
			my @multi_vals = split("\n" , $mult_alle_vars{$ovcf_key} );
			$this_alts{0} = $this_ref;
			foreach(@multi_vals){
				my $this_order = (split("\t" , $_))[0];
				$this_alts{$this_order} = (split("\t" , $_))[5];
				$this_refs{$this_order} = (split("\t" , $_))[4];
			}
			
#			foreach (keys %this_alts){
#				print $_ , " : " , $this_alts{$_} , "\n";
#			}
#			 
			
			#GGC_1_S1:OGT
			#GGC_1_S1:GTR
			foreach my $s (@samples){
				my $this_ogt =  $tmp_row{$s . ":OGT"};
				#print "$this_ogt\n";
				if ( $this_ogt =~ m#^(\d)/(\d)# ) {
					my $gt1 = $1;
					my $gt2 = $2;
					#print "$gt1 == $gt2\n";
					my $new_ogtr = $this_alts{$1} . "/" . $this_alts{$2};
					if($gt1 == 0 && $gt2 == 0){
						$this_refs{$2} . "/" . $this_refs{$2};
					}elsif($gt1 == 0){
						$new_ogtr = $this_refs{$2} . "/" . $this_alts{$2};
					}elsif($gt2 == 0){
						$new_ogtr = $this_alts{$2} . "/" . $this_refs{$2};
					}else{
						$new_ogtr = $this_alts{$1} . "/" . $this_alts{$2};
					}
					$tmp_row{$s . ":GTR"} = $new_ogtr;
				}
				 
			}
			#print line;
			my @outs = ();
			foreach my $h (@hds){
				push(@outs , $tmp_row{$h})
			}
			print $fh_outv join("\t" , @outs) , "\n";
			
		}else{
			print $fh_outv $line , "\n";
		}
		
	}
	return $eff_number_of_vars;
}


#using normalized vcf file(beta , don't use)
sub fix_multiallele_var_vcf{
	my $self = shift;
	my $vcffile = $self->{"name"};
	my $smps = shift;
	my $outvcf = shift;
	my @samples = split("," , $smps); 
#	open(VCF, $vcffile ) or die "$! $vcffile\n";
	my $fh_vcf = IO::File->new($vcffile , q{<}) or usage("$! $vcffile\n");
#	open(OUTV, ">$outvcf") or die "$!";
	my $fh_outv = IO::File->new($outvcf , q{>}) or usage("$! $outvcf\n");
#	my $header = <$fh_vcf>;
#	my @hds = split("\t" , $header);
#	chomp($header);
	my %mult_alle_vars = ();
	
	while(my $line = <$fh_vcf>){
		chomp($line);
		if($line =~ /^#/){
			next;
		}
		my @fields = split("\t", $line);
		
#		my %tmp_row = ();
#		while( my( $index, $value ) = each @hds ) {
#			 $tmp_row{$value} = $fields[$index];
#		}
		
		##CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  Venter-blood-V4
		
		my $this_info = $fields[7];
		my $mal_flag = -1;
		my $ovcf_key = -1;
		if($this_info =~ /.*MULTI_ALLELIC=(.+)/){
				$mal_flag = $1;
		}
        if($this_info =~ /.*OVCF_KEY=(.+?);/){
				$ovcf_key = $1;
				if ($ovcf_key == -1 || $mal_flag == -1 ){
					die "Cannot parse vcf key or multiallele flag, $ovcf_key\n";
				}
		}
		
		
		if($mal_flag != 0){
			$mult_alle_vars{$ovcf_key} = (exists $mult_alle_vars{$ovcf_key}) ? $mult_alle_vars{$ovcf_key} . "\n" . ($mal_flag . "\t" . $line)  : ($mal_flag . "\t" . $line);
		}else{
			#regularget_annoVar_type
		}
	}
	close $fh_vcf;
#	open($fh_vcf, $vcffile ) or die "$! $vcffile\n";
	$fh_vcf = IO::File->new($vcffile , q{<}) or usage("$! $vcffile\n");
#	my $header = <$fh_vcf>;
#	my @hds = split("\t" , $header);
#	chomp($header);
	
	#print header
	#print $fh_outv $header , "\n";
	my $this_ref = ".";
	my $eff_number_of_vars = 0;
	my $header = undef;
	my @hds = ();
	while(my $line = <$fh_vcf>){
		
		chomp($line);
		if($line =~ /^#/){
			print $fh_outv "$line\n";
			if($line =~/^#CHROM/){
				$header = $line;
				@hds = split("\t" , $header);
			}
			next;
		}
		$eff_number_of_vars++;
		my @fields = split("\t", $line);
		
		my %tmp_row = ();
		while( my( $index, $value ) = each @hds ) {
			 $tmp_row{$value} = $fields[$index];
		}
		
		my $this_info = $fields[7];
		my $mal_flag = -1;
		my $ovcf_key = -1;
		if($this_info =~ /.*MULTI_ALLELIC=(.+)/){
				$mal_flag = $1;
		}
        if($this_info =~ /.*OVCF_KEY=(.+?);/){
				$ovcf_key = $1;
				if ($ovcf_key == -1 || $mal_flag == -1 ){
					die "Cannot parse vcf key or multiallele flag, $ovcf_key\n";
				}
		}
		
		
		if($mal_flag != 0){
			$this_ref = $fields[3] unless ($mal_flag > 1);
			my %this_refs = ();
			my %this_alts = ();
			my @multi_vals = split("\n" , $mult_alle_vars{$ovcf_key} );
			$this_alts{0} = $this_ref;
			foreach(@multi_vals){
				my $this_order = (split("\t" , $_))[0];
				$this_alts{$this_order} = (split("\t" , $_))[5];
				$this_refs{$this_order} = (split("\t" , $_))[4];
			}
			
#			foreach (keys %this_alts){
#				print $_ , " : " , $this_alts{$_} , "\n";
#			}
#			 
			
			#GGC_1_S1:OGT
			#GGC_1_S1:GTR
			my $format_ref = $fields[8];
			print "FORMATREF=$format_ref\n";
			#OGT:OZYG:GT:AD:DP:GQ:PL

			my $ogt_index = -1;
			my %format_index = ();
			my @fmt_refs = split(":" , $format_ref);
			while( my( $index, $value ) = each @fmt_refs ) {
				 print "val = $value\n";
				 $format_index{$value} = $index;
			}
			
			foreach my $s (@samples){
				my $this_fmt_val = (exists $tmp_row{$s} ) ? $tmp_row{$s} : die "Sample info for $s missing\n";
				my @fmt_vals = split(':' , $this_fmt_val);
				my $this_ogt = (exists $format_index{'OGT'}) ?    $fmt_vals[$format_index{'OGT'}] : die "OGT is missing for sample $s\n" ;
				
				#print "$this_ogt\n";
				if ( $this_ogt =~ m#^(\d)/(\d)# ) {
					my $gt1 = $1;
					my $gt2 = $2;
					#print "$gt1 == $gt2\n";
					my $new_ogtr = $this_alts{$1} . "/" . $this_alts{$2};
					if($gt1 == 0 && $gt2 == 0){
						$this_refs{$2} . "/" . $this_refs{$2};
					}elsif($gt1 == 0){
						$new_ogtr = $this_refs{$2} . "/" . $this_alts{$2};
					}elsif($gt2 == 0){
						$new_ogtr = $this_alts{$2} . "/" . $this_refs{$2};
					}else{
						$new_ogtr = $this_alts{$1} . "/" . $this_alts{$2};
					}
					
					$fmt_vals[$format_index{'GTR'}] = $new_ogtr;
					$tmp_row{$s} = join(':' , @fmt_vals);
					#$tmp_row{$s . ":GTR"} = $new_ogtr;
				}
				 
			}
			#print line;
			my @outs = ();
			foreach my $h (@hds){
				push(@outs , $tmp_row{$h})
			}
			print $fh_outv join("\t" , @outs) , "\n";
			
		}else{
			print $fh_outv $line , "\n";
		}
	}
	return $eff_number_of_vars;
}



sub split2_small{
	my $self = shift;
	my $vcffile = $self->{"name"};
	my $outdir = shift;	
	my $thre = shift;	
	my $base_vcf = basename($vcffile);
	
	my $header = "#NA";
	open(VCF, $vcffile) or die "$!";
	my $set = 1;
	my $tmp_file = $outdir . "/\." . $base_vcf . "_" . $set;  
	my @split_files = ();
    open(TOUT, ">$tmp_file") or die "$!";
    push(@split_files , $tmp_file);
	my $count = 1;
	my @outs = ();
	while(my $line = <VCF>){
		chomp($line);
		if($line =~ m/^#/){
			$header = $line;
		}else{
			push(@outs , $line);
			$count++; 
		}
		
		if($count == $thre){
			print TOUT $header , "\n";
			print TOUT join("\n" , @outs) , "\n";
			close TOUT;
			$count = 1;
			$set++;
			@outs = ();
			$tmp_file = $outdir . "/\." . $base_vcf . "_" . $set;
			push(@split_files , $tmp_file);
			open(TOUT, ">$tmp_file") or die "$!";
		}
	}
	print TOUT $header , "\n";
	print TOUT join("\n" , @outs) , "\n";
	close TOUT;
	return @split_files;
}

sub merge2_one{
	my $self = shift;
	my $vcffile = $self->{"name"};
	my $outdir = shift;
	my $split_file = shift;	
	my $base_vcf = basename($vcffile);
	my $out_file = $outdir . "/" . $base_vcf . "_annovar.in.txt"; 
    open(AOUT, ">$out_file") or die "$!";
	open(SPN, $split_file) or die "$!";
	my @outs = ();
	while(my $line = <SPN>){
		chomp($line);
		my ($org_file, $this_file) = split("\t" , $line);
		warn "MERGE = $this_file\n";
		open(TIN, $this_file) or die "$!";
		while(my $tin_line = <TIN>){
			chomp($tin_line);
			if($tin_line =~ m/^#/){}else{
				print AOUT $tin_line , "\n";
			}
		}
		close TIN;
		unlink $this_file;
		unlink $org_file;
	}
	close SPN;
	close AOUT;
	return $out_file;
}

#output for help and errors
use Carp;
sub usage {
	my $error = shift;
	if($error){
		croak "\n\nERROR: $error\n\n" ;
	}else{
		#pod2usage();
		exit(-1);
	}
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

1;
