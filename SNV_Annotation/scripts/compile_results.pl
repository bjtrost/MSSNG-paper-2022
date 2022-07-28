#!/tools/perl/5.10.0/bin/perl
use strict;

our $AUTHOR =   '$Author: Thomas Nalpathamkalam <thomas.nalpathamkalam@sickkids.ca> $';

## Module to compile and make final outputs from various annovar and custom modules.
## Called by the main script. 
## Loads all intermediate files to memory and compare with main input file.

use Getopt::Std;
$Getopt::Long::ignorecase = 0;    #case sensitive
use Data::Dumper qw(Dumper);
use List::Util qw( min max sum reduce);
use Scalar::Util qw(looks_like_number);
use List::MoreUtils qw(uniq);
use File::Basename;
use FindBin;
use lib "$FindBin::Bin";
use PipelineExtra;
use Time::Piece;
use IO::File;
use POSIX;


my %options = ();
getopts( "m:t:v:", \%options );
&usage0("Please specify pipeline_meta file\n") unless ( -s $options{m} );

my $mstart_time = localtime;

my %mta      = ();
my %exit_noz = ();
my %anno_empty = ();
my %files_parallel = ();
my %pr_names = ();
my $proc = 1;
my %gene_stats = ();
my %eff_categories = ();
my %full_summary = ();
my %validated = ();

#hard coded , to be fixed, July 20, 2016
#my $spidex_header = "gene:strand:transcript:exon_number:location:cds_type:ss_dist:raw_dpsi:wt_psi:dlogit:is_stop_gain:sequence_event_type:dpsi";
#1       861180  861181  G/A     0.5983:1.226:SAMD11:+:NM_152486:2:intronic:5pUTR:121 (old)
#chr1    925800  925801  G/A     0.5983:1.226:SAMD11:+:NM_152486:2:intronic:5pUTR:121 (old)
#chr1    861180  861181  G/A     SAMD11:NM_152486:2:121:96.065185:splice_donor_variant:3.411093 ( new)
#my $spidex_header = "gene:transcript:exon_number:ss_dist:wt_psi:sequence_event_type:dpsi" ; #(new)
#my $spidex_header = "gene:transcript:exon_number:ss_dist:wt_psi:sequence_event_type:dpsi";  #(new)
my $spidex_header = "dpsi:dpsi_z:gene:strand:transcript:exonN:seqType:effType:spliceDist"; #old
#0.3090:0.904:SAMD11:+:NM_152486:2:intronic:5pUTR:107

#SYMBOL=ARHGEF16;STRAND=+;TYPE=E;DIST=-51;DS_AG=0.0014;DS_AL=0.0000;DS_DG=0.4581;DS_DL=0.0000;DP_AG=-47;DP_AL=-2;DP_DG=-3;DP_DL=14

#v3.1
#SpliceAI=G|OR4F5|0.01|0.08|0.00|0.00|-10|26|-28|-25
#ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL

my $spliceAIsnv_header = "ALLELE:SYMBOL:DS_AG:DS_AL:DS_DG:DS_DL:DP_AG:DP_AL:DP_DG:DP_DL";
my @cgInt_header = ("AF_parents","AF_parents_cgi","AF_parents_cgi_eur","AF_parents_ilm","AF_parents_ilm_eur","Key_hg19","refMatchCallsetsCount_parents","refMatchCallsetsCount_parents_cgi","refMatchCallsetsCount_parents_cgi_eur","refMatchCallsetsCount_parents_ilm","refMatchCallsetsCount_parents_ilm_eur","cgi_max_int_freq","cgi_min_int_ref");
my @illInt_header = ("AC_ilm","AN_ilm","AF_ilm","refMatchSamples_ilm","homMatchSamples_ilm","noCallSamples_ilm", "homMatchSamples_fathers_ilm", "AC_parents_ilm","AN_parents_ilm","AF_parents_ilm","refMatchSamples_parents_ilm","homMatchSamples_parents_ilm","noCallSamples_parents_ilm","AC_parents_ilm_eur","AN_parents_ilm_eur","AF_parents_ilm_eur","refMatchSamples_parents_ilm_eur","noCallSamples_parents_ilm_eur","AC_parents_ilm_non-eur","AN_parents_ilm_non-eur","AF_parents_ilm_non-eur","refMatchSamples_parents_ilm_non-eur","noCallSamples_parents_ilm_non-eur","AC_parents_ilm-hiseqx","AN_parents_ilm-hiseqx","AF_parents_ilm-hiseqx","refMatchSamples_parents_ilm-hiseqx","noCallSamples_parents_ilm-hiseqx","AC_parents_ilm-non-hiseqx","AN_parents_ilm-non-hiseqx","AF_parents_ilm-non-hiseqx","refMatchSamples_parents_ilm-non-hiseqx","noCallSamples_parents_ilm-non-hiseqx","AC_parents_ilm-hiseqx_eur","AN_parents_ilm-hiseqx_eur","AF_parents_ilm-hiseqx_eur","refMatchSamples_parents_ilm-hiseqx_eur","noCallSamples_parents_ilm-hiseqx_eur","AC_parents_ilm-hiseqx_non-eur","AN_parents_ilm-hiseqx_non-eur","AF_parents_ilm-hiseqx_non-eur","refMatchSamples_parents_ilm-hiseqx_non-eur","noCallSamples_parents_ilm-hiseqx_non-eur","ilm_max_int_freq","ilm_min_int_ref"); 

my $pat_weld = qr/W_AllFreq=(.+?);W_CalledFreq=(.+?);.+W_11s=(.+?);W_Hs=(.+?);W_Ls=(.+?);/;
my $pat_cg2 = qr/AllFreq=(.+?);CalledFreq=(.+?);.+11s=(.+?);Hs=(.+?);Ls=(.+?);/;
my $pat_mssng = qr/MSSNG_KEY=(.+?);.+/;


#121     spx_gene
#122     spx_transcript
#123     spx_exon_number
#124     spx_ss_dist
#125     spx_wt_psi
#126     spx_sequence_event_type
#127     spx_dpsi


#my $spidex_header = "dpsi:dpsi_z:gene:strand:transcript:exonN:seqType:effType:spliceDist"; #old


#for bitwise impact

my %EFFECTS  = (
	 'frameshift' => 0,
	 'stop_gain' => 3,
	 'splice_site' => 6,
	 'lof' => 9,
	 'missense' => 12,
	 'other' => 15,
	 'reg_dec_exon' => 18,
	 'reg_inc_exon' => 21,
	 'utr' => 24,
	 'non_coding' => 27
);

#EFFECTS = [
#        [:frameshift,   0, 'Frameshift'],
#        [:stop_gain,    3, 'Stop Gain'],
#        [:splice_site,  6, 'Splice Site'],
#        [:lof,          9, 'LOF'],
#        [:missense,     12, 'Missense'],
#        [:other,        15, 'Other'],
#        [:reg_dec_exon, 18, 'Predicted Splicing'],
#        [:reg_inc_exon, 21, 'Splicing Reg Pos'],
#        [:utr,          24, 'UTR'],
#        [:non_coding,   27, 'Non-coding RNA gene']


my %EFFECTS_str = ();

foreach (keys %EFFECTS){
	$EFFECTS_str{$EFFECTS{$_}} = $_;
}


my %IMPACTS  = (
	 'high' => 0,
	 'medium' => 1,
	 'low' => 2,
);

my %IMPACTS_str = ();
foreach (keys %IMPACTS){
	$IMPACTS_str{$IMPACTS{$_}} = $_;
}


my %EXONIC = (
    'exonic' => 1,
    'exonic;splicing' => 1,
    'splicings' => 1,
);


my %NON_CODING = (
    'ncRNA_exonic' => 1,
    'ncRNA_splicing' => 1,
    'ncRNA_exonic;ncRNA_splicing' => 1,
);

my %FRAMESHIFT = (
    'frameshift deletion' => 1,
    'frameshift insertion' => 1,
    'frameshift substitution' => 1,
);

my %stop_gain = (
    'stopgain SNV' => 1,
    'stopgain' => 1,
    'frameshift substitution' => 1,
);

my %other = (
    'nonframeshift deletion' => 1,
    'nonframeshift insertion' => 1,
    'nonframeshift substitution' => 1,
    'stoploss SNV' => 1,
    'stoploss' => 1,
);


#gnomeAD_exome_
my %gnADgn_reqs = (
    'gnomAD_genome_AF' => 'gnomAD_genome_ALL',
    'gnomAD_genome_AF_afr' => 'gnomAD_genome_AFR',
    'gnomAD_genome_AF_amr' => 'gnomAD_genome_AMR',
    'gnomAD_genome_AF_asj' => 'gnomAD_genome_ASJ',
    'gnomAD_genome_AF_eas' => 'gnomAD_genome_EAS',
    'gnomAD_genome_AF_fin' => 'gnomAD_genome_FIN',
    'gnomAD_genome_AF_nfe' => 'gnomAD_genome_NFE',
    'gnomAD_genome_AF_oth' => 'gnomAD_genome_OTH',
    
);


my %gnADex_reqs = (
    'gnomAD_exome_AF' => 'gnomAD_exome_ALL',
    'gnomAD_exome_AF_afr' => 'gnomAD_exome_AFR',
    'gnomAD_exome_AF_amr' => 'gnomAD_exome_AMR',
    'gnomAD_exome_AF_asj' => 'gnomAD_exome_ASJ',
    'gnomAD_exome_AF_eas' => 'gnomAD_exome_EAS',
    'gnomAD_exome_AF_fin' => 'gnomAD_exome_FIN',
    'gnomAD_exome_AF_nfe' => 'gnomAD_exome_NFE',
    'gnomAD_exome_AF_oth' => 'gnomAD_exome_OTH',
    'gnomAD_exome_AF_sas' => 'gnomAD_exome_SAS',
    
);

my @gnADgn30_reqs = ("gnomAD_genome30_FILTER", "gnomAD_genome30_AC","gnomAD_genome30_AC_female","gnomAD_genome30_AC_male","gnomAD_genome30_AF","gnomAD_genome30_AF_afr", "gnomAD_genome30_AF_nfe", "gnomAD_genome30_AF_ami","gnomAD_genome30_AF_amr","gnomAD_genome30_AF_asj","gnomAD_genome30_AF_eas","gnomAD_genome30_AF_female","gnomAD_genome30_AF_fin","gnomAD_genome30_AF_male","gnomAD_genome30_AF_oth","gnomAD_genome30_AF_raw","gnomAD_genome30_AF_sas","gnomAD_genome30_AN","gnomAD_genome30_AN_female","gnomAD_genome30_AN_male","gnomAD_genome30_faf95_adj","gnomAD_genome30_faf95_afr","gnomAD_genome30_faf95_amr","gnomAD_genome30_faf95_eas","gnomAD_genome30_faf95_nfe","gnomAD_genome30_faf95_sas","gnomAD_genome30_nhomalt","gnomAD_genome30_nhomalt_female","gnomAD_genome30_nhomalt_male");
my @gnADgn211_reqs = ("gnomAD_genome211_FILTER", "gnomAD_genome211_AC","gnomAD_genome211_AC_female","gnomAD_genome211_AC_male","gnomAD_genome211_AF","gnomAD_genome211_AF_afr", "gnomAD_genome211_AF_nfe", "gnomAD_genome211_AF_amr","gnomAD_genome211_AF_asj","gnomAD_genome211_AF_eas","gnomAD_genome211_AF_female","gnomAD_genome211_AF_fin","gnomAD_genome211_AF_male","gnomAD_genome211_AF_oth","gnomAD_genome211_AF_raw","gnomAD_genome211_AF_sas","gnomAD_genome211_AN","gnomAD_genome211_AN_female","gnomAD_genome211_AN_male","gnomAD_genome211_faf95_afr","gnomAD_genome211_faf95_amr","gnomAD_genome211_faf95_eas","gnomAD_genome211_faf95_nfe","gnomAD_genome211_faf95_sas","gnomAD_genome211_nhomalt","gnomAD_genome211_nhomalt_female","gnomAD_genome211_nhomalt_male", "gnomAD_genome211_controls_popmax");
my @gnADex211_reqs = ("gnomAD_exome211_FILTER", "gnomAD_exome211_AC","gnomAD_exome211_AC_female","gnomAD_exome211_AC_male","gnomAD_exome211_AF","gnomAD_exome211_AF_afr", "gnomAD_exome211_AF_nfe", "gnomAD_exome211_AF_amr","gnomAD_exome211_AF_asj","gnomAD_exome211_AF_eas","gnomAD_exome211_AF_female","gnomAD_exome211_AF_fin","gnomAD_exome211_AF_male","gnomAD_exome211_AF_oth","gnomAD_exome211_AF_raw","gnomAD_exome211_AF_sas","gnomAD_exome211_AN","gnomAD_exome211_AN_female","gnomAD_exome211_AN_male","gnomAD_exome211_faf95_afr","gnomAD_exome211_faf95_amr","gnomAD_exome211_faf95_eas","gnomAD_exome211_faf95_nfe","gnomAD_exome211_faf95_sas","gnomAD_exome211_nhomalt","gnomAD_exome211_nhomalt_female","gnomAD_exome211_nhomalt_male", "gnomAD_exome211_controls_popmax");


#hard-coded, make this 1 if you want to map the columns based on above hash
my $gnAD_fix = 0;

my @VALID_CHR = 1 .. 22;
#@a = map { (my $s = $_) =~ s/^_+//; $s } @a;
#@VALID_CHR =  map { "chr" . $_ } @VALID_CHR;

foreach (1 .. 22){
	push(@VALID_CHR , "chr" . $_);
}

push(@VALID_CHR , ("X", "Y", "MT" , "chrM", "chrX" , "chrY") );  # add more if you've db for special chromosomes

my %valid_chrs = ();
$valid_chrs{$_}++ for (@VALID_CHR);

my @spx_valna = ();
foreach (split(":" , $spidex_header)){
	push(@spx_valna , "NA");
}

my @spliceAIsnv_valna = ();
foreach (split(":" , $spliceAIsnv_header)){
	push(@spliceAIsnv_valna , "NA");
}

my $mr_tmp_dir = "/tmp";
if(defined $options{t}){
       $mr_tmp_dir = $options{t};
}
my $mr_core = 2;

my $vcfout = 0;
if(defined $options{v}){
	$vcfout =  $options{v};
}
warn "VCF = $vcfout  -- $options{v} \n";
#open( MTA, $options{m} );
my $fh_mta = IO::File->new($options{m} , q{<}) or die "$! meta file\n";

#merge meta file for stats, added oct16,2016

my $merge_metafile = $mr_tmp_dir . "/final_merge.meta";
my $fh_mrg_mta = IO::File->new($merge_metafile , q{>}) or die "$! $merge_metafile\n";

while ( my $line = <$fh_mta> ) {
	chomp($line);
	my ( $key, $val ) = split( "=",  $line );
	my ( $mle, $st )  = split( '->', $key );
	if ( $st eq "exit_status" ) {
		#usage0("Module: $mle  has non-zero exit status\n") unless ($val == 0);
		$exit_noz{$mle} = 1 unless ( $val == 0 );
	}
	if ( $st eq "out_state" ) {
		$anno_empty{$mle} = 1 if ( $val ne "non-zero" );
	}
	$mta{$key} = $val;
}

if ( scalar( keys %exit_noz ) != 0 ) {
	my @erron_mod = ();
	foreach my $nzk (keys %exit_noz){
		push(@erron_mod, ">>" . $nzk);
		my $this_err = $nzk . "->error_file";
		my $this_err_file = $mta{$this_err};
		my $fh_mod_err = IO::File->new($this_err_file , q{<}) or die "$! Cannot find $this_err_file\n";
		push(@erron_mod, <$fh_mod_err>);
	}
	
	usage0( "Following module(s) has non-zero exit status\n"
		  . join("\n", keys %exit_noz) . "\n\n" 	
		  . join( "\n", ( @erron_mod ) )
		  . "\n" );

}

#output config file
my $ouput_config_file = ( exists $mta{'pipeline->output_config'} ) ? $mta{'pipeline->output_config'}  : -1;
if ( $ouput_config_file == -1 ) {
	&usage0("pipeline_meta: Missing Output config file");
}

my $pipeline_tmp_dir =
  ( exists $mta{'pipeline->tmp_dir'} )
  ? $mta{'pipeline->tmp_dir'}
  : -1;    # &usage("Missing temp directory  \n") ;
if ( $pipeline_tmp_dir == -1 ) {
	&usage0("pipeline_meta: Missing temp directory");
}

my $pipeline_work_dir =
  ( exists $mta{'pipeline->work_dir'} )
  ? $mta{'pipeline->work_dir'}
  : -1;    # &usage("Missing temp directory  \n") ;
if ( $pipeline_work_dir == -1 ) {
	&usage0("pipeline_meta: Missing work directory info");
}

my $pipeline_chrInfo =   ( exists $mta{'pipeline->chr_info'} )
  ? $mta{'pipeline->chr_info'}
  : -1;    # &usage("Missing temp directory  \n") ;
if ( $pipeline_chrInfo == -1 ) {
	&usage0("pipeline_meta: Missing chrInfo file");
}

#open(PCHR, $pipeline_chrInfo) or die "$! $pipeline_chrInfo\n";
my $fh_pchr = IO::File->new($pipeline_chrInfo , q{<}) or die "$! $pipeline_chrInfo\n";

my @pipeline_valid_chrs = ();
my @pipeline_valid_chrsAll = ();

while(my $pline = <$fh_pchr>){
	chomp($pline);
	my( $p_chr, $p_chrF, $p_chrW )  = split("\t" , $pline);
	push(@pipeline_valid_chrsAll , $p_chr);	
	if(exists $valid_chrs{$p_chr}){
		push(@pipeline_valid_chrs , $p_chr);
	}
}
close $fh_pchr;


#pipeline->annovar_in
my $annovar_in =
  ( exists $mta{'pipeline->annovar_in'} )
  ? $mta{'pipeline->annovar_in'}
  : -1;    #  &usage("Missing annovar-input file  \n")  ;
if ( $annovar_in == -1 ) {
	&usage0("pipeline_meta: Missing annovar-input file");
}
my $annovar_in_basename = basename($annovar_in);

#pipeline->org_input
my $orgvcf_in =
  ( exists $mta{'pipeline->org_input'} )
  ? $mta{'pipeline->org_input'}
  : -1;    #  &usage("Missing annovar-input file  \n")  ;
if ( $orgvcf_in == -1 ) {
	&usage0("pipeline_meta: Missing org_input file");
}

my $orgvcf_in_basename = basename($orgvcf_in);


my $warn_file           =
  $pipeline_tmp_dir . "/warnings.txt";
open( my $wfh, '>', "$warn_file" );
if( scalar (keys %anno_empty) != 0){
	print $wfh "Following annotations are missing: \n " , join( "," ,  sort (keys %anno_empty)) , "\n";
}

#Gene->var_output
my $refgene_var =
  ( exists $mta{'Gene->var_output'} )
  ? $mta{'Gene->var_output'}
  : -1;    # &usage2("Missing gene-variant annotation file \n")  ;
if ( $refgene_var == -1 ) {
	&usage2("pipeline_meta: Missing gene-variant annotation file");
}

if ( $mta{'Gene->exit_status'} != 0 ) {
	usage("Gene annotation has non-zero exit status\n");
}else{
	$files_parallel{$mta{'Gene->var_output'}} = "g";
	$pr_names{$mta{'Gene->var_output'}} = "Gene-var";
}

# unless ( -s $options{x} );
my $refgene_vex =
  ( exists $mta{'Gene->exonic_output'} )
  ? $mta{'Gene->exonic_output'}
  : -1;    # &usage2("Missing gene-exonic variant annotation file \n")  ;
if ( $refgene_vex == -1 ) {
	&usage2("pipeline_meta: Missing gene-exonic variant annotation file ");
}else{
	$files_parallel{$mta{'Gene->exonic_output'}} = "x";
	$pr_names{$mta{'Gene->exonic_output'}} = "Gene-exonic";
}

#&usage2("Missing UCSC gene-variant annotation file \n") unless ( -s $options{N} );
#Gene-UCSC->var_output
my $ucgene_var =
  ( exists $mta{'UCSC-Gene->var_output'} )
  ? $mta{'UCSC-Gene->var_output'}
  : -1;    # &usage2("Missing UCSC gene-variant annotation file \n")  ;
if ( $ucgene_var == -1 ) {
	&usage2("pipeline_meta: Missing UCSC gene-variant annotation file ");
}

if ( $mta{'UCSC-Gene->exit_status'} != 0 ) {
	usage("UCSC-Gene annotation has non-zero exit status");
}else{
	$files_parallel{$mta{'UCSC-Gene->var_output'}} = "e";
	$pr_names{$mta{'UCSC-Gene->var_output'}} = "UCSC-Gene-var";
}

#&usage2("Missing UCSC gene-exonic variant annotation file \n")
# unless ( -s $options{X} );
my $ucgene_vex =
  ( exists $mta{'UCSC-Gene->exonic_output'} )
  ? $mta{'UCSC-Gene->exonic_output'}
  : -1;    #  &usage2("Missing UCSC gene-exonic variant annotation file \n")  ;
if ( $ucgene_vex == -1 ) {
	&usage2(
		"pipeline_meta: Missing UCSC gene-exonic variant annotation file");
}else{
	$files_parallel{$mta{'UCSC-Gene->exonic_output'}} = "ux";
	$pr_names{$mta{'UCSC-Gene->exonic_output'}} = "UCSC-Gene-exonic";
}

#refgene-relaxed
my $gcode_var =
  ( exists $mta{'Gene-Rlx->var_output'} )
  ? $mta{'Gene-Rlx->var_output'}
  : -1;    # &usage2("Missing UCSC gene-variant annotation file \n")  ;
if ( $gcode_var == -1 ) {
	&usage2("pipeline_meta: Missing Gene-Rlx gene-variant annotation file ");
}

if ( $mta{'Gene-Rlx->exit_status'} != 0 ) {
	usage("Gene-Rlx annotation has non-zero exit status");
}else{
	$files_parallel{$mta{'Gene-Rlx->var_output'}} = "e";
	$pr_names{$mta{'Gene-Rlx->var_output'}} = "Gene-Rlx";
}

#ambigeous sites
#my $db_hg39Ambg =
#  ( exists $mta{'pipeline->hg38_ambiguous_sites'} )
#  ? $mta{'pipeline->hg38_ambiguous_sites'}
#  : -1;    # &usage2("Missing UCSC gene-variant annotation file \n")  ;
#if ( $db_hg39Ambg == -1 ) {
#	&usage("pipeline_meta: Missing hg38_ambiguous_sites file  \n");
#}



#&usage2("Missing UCSC gene-exonic variant annotation file \n")
# unless ( -s $options{X} );
my $gcode_vex =
  ( exists $mta{'Gene-Rlx->exonic_output'} )
  ? $mta{'Gene-Rlx->exonic_output'}
  : -1;    #  &usage2("Missing UCSC gene-exonic variant annotation file \n")  ;
if ( $gcode_vex == -1 ) {
	&usage2(
		"pipeline_meta: Missing Gene-Rlx-Gene exonic variant annotation file ");
}else{
	$files_parallel{$mta{'Gene-Rlx->exonic_output'}} = "ux";
	$pr_names{$mta{'Gene-Rlx->exonic_output'}} = "Gene-Rlx-exonic";
}



#refgene-select
my $gsel_var =
  ( exists $mta{'Gene-Select->var_output'} )
  ? $mta{'Gene-Select->var_output'}
  : -1;    # &usage2("Missing UCSC gene-variant annotation file \n")  ;
if ( $gsel_var == -1 ) {
	&usage("pipeline_meta: Missing Gene-Select gene-variant annotation file ");
}

if ( $mta{'Gene-Select->exit_status'} != 0 ) {
	usage("Gene-Select annotation has non-zero exit status");
}else{
	$files_parallel{$mta{'Gene-Select->var_output'}} = "e";
	$pr_names{$mta{'Gene-Select->var_output'}} = "Gene-Select";
}


#&usage2("Missing UCSC gene-exonic variant annotation file \n")
# unless ( -s $options{X} );
my $gsel_vex =
  ( exists $mta{'Gene-Select->exonic_output'} )
  ? $mta{'Gene-Select->exonic_output'}
  : -1;    #  &usage2("Missing UCSC gene-exonic variant annotation file \n")  ;
if ( $gsel_vex == -1 ) {
	&usage(
		"pipeline_meta: Missing Gene-Select exonic variant annotation file ");
}else{
	$files_parallel{$mta{'Gene-Select->exonic_output'}} = "ux";
	$pr_names{$mta{'Gene-Select->exonic_output'}} = "Gene-Select-exonic";
}


#Gene-Appris
my $gappris_var =
  ( exists $mta{'Gene-Appris->var_output'} )
  ? $mta{'Gene-Appris->var_output'}
  : -1;    # &usage2("Missing UCSC gene-variant annotation file \n")  ;
if ( $gappris_var == -1 ) {
	&usage2("pipeline_meta: Missing Gene-Appris gene-variant annotation file ");
}

if ( $mta{'Gene-Appris->exit_status'} != 0 ) {
	usage("Gene-Appris annotation has non-zero exit status");
}else{
	$files_parallel{$mta{'Gene-Appris->var_output'}} = "e";
	$pr_names{$mta{'Gene-Appris->var_output'}} = "Gene-Appris";
}

my $gppris_vex =
  ( exists $mta{'Gene-Appris->exonic_output'} )
  ? $mta{'Gene-Appris->exonic_output'}
  : -1;    #  &usage2("Missing UCSC gene-exonic variant annotation file \n")  ;
if ( $gppris_vex == -1 ) {
	&usage2(
		"pipeline_meta: Missing Gene-Appris exonic variant annotation file ");
}else{
	$files_parallel{$mta{'Gene-Appris->exonic_output'}} = "ux";
	$pr_names{$mta{'Gene-Appris->exonic_output'}} = "Gene-Appris-exonic";
}



#gene_info->ucGene
#&usage("Please specify ucsc knowngene file \n")           unless ( -s $options{P} );
my $uc_knowngene =
  ( exists $mta{'gene_info->ucGene'} )
  ? $mta{'gene_info->ucGene'}
  : -1;    # &usage2("Please specify ucsc knowngene file \n")  ;
if ( $uc_knowngene == -1 ) {
	&usage("pipeline_meta : Please specify ucsc knowngene file ");
}

#CG->output
#&usage2("Missing cg annotation file \n")              unless ( -s $options{c} );
my $cg_anno =
  ( exists $mta{'CG->output'} )
  ? $mta{'CG->output'} 
  : -1;    # &usage2("Missing cg annotation file\n")  ;
if ( $cg_anno == -1 ) {
	&usage2("pipeline_meta:  Missing cg annotation file ");
}
if ( $mta{'CG->exit_status'} != 0 ) {
	usage("CG annotation has non-zero exit status");
}else{
	$files_parallel{$mta{'CG->output'}} = "f";
	$pr_names{$mta{'CG->output'}} = "CG";
}

#Welderly->output
#&usage2("Missing CG WLD annotation file \n")              unless ( -s $options{L} );
my $wld_anno =
  ( exists $mta{'Welderly->output'} )
  ? $mta{'Welderly->output'}
  : -1;    # &usage2("Missing CG WLD annotation file\n")  ;
if ( $wld_anno == -1 ) {
	&usage2("pipeline_meta:  Missing CG WLD annotation file  ");
}
if ( $mta{'Welderly->exit_status'} != 0 ) {
	usage("Welderly annotation has non-zero exit status");
}else{
	$files_parallel{$mta{'Welderly->output'}} = "f";
	$pr_names{$mta{'Welderly->output'}} = "Welderly";
}
#
#cg1KB->output
##&usage2("Missing CG 1KB annotation file \n")              unless ( -s $options{T} );
my $cg1KB_anno =
  ( exists $mta{'cg1KB->output'} )
  ? $mta{'cg1KB->output'}
  : -1;    #  &usage2("Missing CG 1KB annotation file\n")  ;
if ( $cg1KB_anno == -1 ) {
	&usage2("pipeline_meta:  Missing CG 1KB annotation file ");
}
if ( $mta{'cg1KB->exit_status'} != 0 ) {
	usage("cg1KB annotation has non-zero exit status");
}else{
	$files_parallel{$mta{'cg1KB->output'}} = "f";
	$pr_names{$mta{'cg1KB->output'}} = "cg1KB";
}

#dbsnp->output
#&usage2("Missing dbsnp annotation file \n")           unless ( -s $options{d} );
my $dbsnp_anno =
  ( exists $mta{'dbsnp->output'} )
  ? $mta{'dbsnp->output'}
  : -1;    # &usage2("Missing dbsnp annotation file\n")  ;
if ( $dbsnp_anno == -1 ) {
	&usage2("pipeline_meta:  Missing dbsnp annotation file  ");
}
if ( $mta{'dbsnp->exit_status'} != 0 ) {
	usage("dbsnp annotation has non-zero exit status");
}else{
	$files_parallel{$mta{'dbsnp->output'}} = "f";
	$pr_names{$mta{'dbsnp->output'}} = "dbsnp";
}

#dbsnp common->output
#&usage2("Missing dbsnpCommon annotation file \n")           unless ( -s $options{b} );
my $dbsnpC_anno =
  ( exists $mta{'dbsnpCommon->output'} )
  ? $mta{'dbsnpCommon->output'}
  : -1;    # &usage2("Missing dbsnpCommon annotation file\n")  ;
if ( $dbsnpC_anno == -1 ) {
	&usage2("pipeline_meta:  Missing dbsnpCommon annotation file ");
}
if ( $mta{'dbsnpCommon->exit_status'} != 0 ) {
	usage("dbsnp_common annotation has non-zero exit status");
}else{
	$files_parallel{$mta{'dbsnpCommon->output'}} = "f";
	$pr_names{$mta{'dbsnpCommon->output'}} = "dbsnpCommon";
}

#1000g->output
#&usage2("Missing 1000g annotation file \n")           unless ( -s $options{t} );
my $thg_anno =
  ( exists $mta{'1000g->output'} )
  ? $mta{'1000g->output'}
  : -1;    # &usage2("Missing 1000g annotation(ALL) file\n")  ;
if ( $thg_anno == -1 ) {
	&usage2("pipeline_meta: Missing 1000g annotation(ALL) file  ");
}
if ( $mta{'1000g->exit_status'} != 0 ) {
	usage("1000g annotation has non-zero exit status");
}else{
	$files_parallel{$mta{'1000g->output'}} = "f";
	$pr_names{$mta{'1000g->output'}} = "1000g";
}
#1000g amr->output
#&usage2("Missing 1000g AMR annotation file \n")           unless ( -s $options{A} );
my $thg_amr_anno =
  ( exists $mta{'1000g_AMR->output'} )
  ? $mta{'1000g_AMR->output'}
  : -1;    # &usage2("Missing 1000g AMR annotation file\n")  ;
if ( $thg_amr_anno == -1 ) {
	&usage2("pipeline_meta: Missing 1000g annotation(AMR) file ");
}else{
	$files_parallel{$mta{'1000g_AMR->output'}} = "f";
	$pr_names{$mta{'1000g_AMR->output'}} = "1000g_AMR";
}

if ( $mta{'1000g_AMR->exit_status'} != 0 ) {
	usage("1000g_AMR annotation has non-zero exit status");
}

#1000g eas->output
#&usage2("Missing 1000g ASN(or EAS) annotation file \n")           unless ( -s $options{B} );
my $thg_eas_anno =
  ( exists $mta{'1000g_EAS->output'} )
  ? $mta{'1000g_EAS->output'}
  : -1;    # &usage2("Missing 1000g ASN(or EAS) annotation file\n")  ;
if ( $thg_eas_anno == -1 ) {
	&usage2("pipeline_meta: Missing 1000g annotation(EAS) file ");
}else{
	$files_parallel{$mta{'1000g_EAS->output'}} = "f";
	$pr_names{$mta{'1000g_EAS->output'}} = "1000g_EAS";
}

if ( $mta{'1000g_EAS->exit_status'} != 0 ) {
	usage("1000g_EAS annotation has non-zero exit status");
}

#1000g eur->output
#&usage2("Missing 1000g EUR annotation file \n")           unless ( -s $options{C} );
my $thg_eur_anno =
  ( exists $mta{'1000g_EUR->output'} )
  ? $mta{'1000g_EUR->output'}
  : -1;    # &usage2("Missing 1000g EUR annotation file\n")  ;
if ( $thg_eur_anno == -1 ) {
	&usage2("pipeline_meta: Missing 1000g annotation(EUR) file  ");
}else{
	$files_parallel{$mta{'1000g_EUR->output'}} = "f";
	$pr_names{$mta{'1000g_EUR->output'}} = "1000g_EUR";
}
if ( $mta{'1000g_EUR->exit_status'} != 0 ) {
	usage("1000g_EUR annotation has non-zero exit status");
}

#1000g afr->output
#&usage2("Missing 1000g AFR annotation file \n")           unless ( -s $options{D} );
my $thg_afr_anno =
  ( exists $mta{'1000g_AFR->output'} )
  ? $mta{'1000g_AFR->output'}
  : -1;    #  &usage2("Missing 1000g AFR annotation file\n")  ;
if ( $thg_afr_anno == -1 ) {
	&usage2("pipeline_meta: Missing 1000g annotation(AFR) file  ");
}else{
	$files_parallel{$mta{'1000g_AFR->output'}} = "f";
	$pr_names{$mta{'1000g_AFR->output'}} = "1000g_AFR";
}
if ( $mta{'1000g_AFR->exit_status'} != 0 ) {
	usage("1000g afr annotation has non-zero exit status");
}

#1000g sas->output
#&usage2("Missing 1000g SAS annotation file \n")           unless ( -s $options{O} );
my $thg_sas_anno =
  ( exists $mta{'1000g_SAS->output'} )
  ? $mta{'1000g_SAS->output'}
  : -1;    # &usage2("Missing 1000g SAS annotation file\n")  ;
if ( $thg_sas_anno == -1 ) {
	&usage2("pipeline_meta: Missing 1000g annotation(SAS) file ");
}else{
	$files_parallel{$mta{'1000g_SAS->output'}} = "f";
	$pr_names{$mta{'1000g_SAS->output'}} = "1000g_SAS";
}
if ( $mta{'1000g_SAS->exit_status'} != 0 ) {
	usage("1000g sas annotation has non-zero exit status");
}

my $exac_anno =
  ( exists $mta{'exac->output'} )
  ? $mta{'exac->output'}
  : -1;    # &usage2("Missing ExAc annotation file \n")  ;
if ( $exac_anno == -1 ) {
	&usage2("pipeline_meta: Missing ExAc annotation file ");
}else{
	$files_parallel{$mta{'exac->output'}} = "f";
	$pr_names{$mta{'exac->output'}} = "exac";
}
if ( $mta{'exac->exit_status'} != 0 ) {
	usage("exac annotation has non-zero exit status");
}

#gnomAD annovar, added apr25,2017
#gnomad_exome211->output
my $gnADex_anno =
  ( exists $mta{'gnomad_exome211->output'} )
  ? $mta{'gnomad_exome211->output'}
  : -1;    # &usage2("Missing ExAc annotation file \n")  ;
if ( $gnADex_anno == -1 ) {
	&usage2("pipeline_meta: Missing gnomad_exome annotation file ");
}else{
	$files_parallel{$mta{'gnomad_exome211->output'}} = "f";
	$pr_names{$mta{'gnomad_exome211->output'}} = "gnomad_exome211";
}
if ( $mta{'gnomad_exome211->exit_status'} != 0 ) {
	usage("gnomad_exome211 annotation has non-zero exit status");
}

#my $gnADex_db =  ( exists $mta{'gnomad_exome->database'} ) ? $mta{'gnomad_exome->database'} : -1;
#if ( $gnADex_db == -1 ) {
#	&usage("pipeline_meta: Missing gnomeAD_exome(annovar generated) file \n");
#}

#my $gnADgn_anno211 =
#  ( exists $mta{'gnomad_genome211->output'} )
#  ? $mta{'gnomad_genome211->output'}
#  : -1;    # &usage2("Missing ExAc annotation file \n")  ;
#if ( $gnADgn_anno211 == -1 ) {
#	&usage2("pipeline_meta: Missing gnomad_genome211 annotation file ");
#}else{
#	$files_parallel{$mta{'gnomad_genome211->output'}} = "f";
#	$pr_names{$mta{'gnomad_genome211->output'}} = "gnomad_genome211";
#}
#if ( $mta{'gnomad_genome211->exit_status'} != 0 ) {
#	usage("gnomad_genome211 annotation has non-zero exit status\n");
#}
#my $gnADgn211_db =  ( exists $mta{'gnomad_genome211->database'} ) ? $mta{'gnomad_genome211->database'} : -1;
#if ( $gnADgn211_db == -1 ) {
#	&usage("pipeline_meta: Missing gnomeAD_genome211(annovar generated) file \n");
#}

#gerpwgs chrwise
foreach my $this_chr (@pipeline_valid_chrs){
	#gerpwgs:12->output
	$this_chr =~ s/chr//;
	my $out_key = "gerpwgs:" .  $this_chr; # . "->output";
	my $this_gerpwgs_anno =  ( exists $mta{$out_key . '->output'} )   ? $mta{$out_key. '->output'} : -1;
	if($this_gerpwgs_anno == -1 ){
		&usage2("pipeline_meta: Missing gerpwgs $this_chr annotation file \n");
	}else{
		$files_parallel{$mta{$out_key . '->output'}} = "f";
		$pr_names{$mta{$out_key . '->output'}} = $out_key;
	}
	if ( $mta{$out_key . '->exit_status'} != 0 ) {
		usage("gerpwgs $this_chr annotation has non-zero exit status\n");
	}
}



my $dbnsfp_db = "NA";
#dbnsfp chrwise
foreach my $this_chr (@pipeline_valid_chrs){
	#gerpwgs:12->output
	$this_chr =~ s/chr//;
	my $out_key = "dbnsfp:" .  $this_chr; # . "->output";
	my $this_dbnsfp_anno =  ( exists $mta{$out_key . '->output'} )   ? $mta{$out_key. '->output'} : -1;
	if($this_dbnsfp_anno == -1 ){
		&usage2("pipeline_meta: Missing dbnsfp $this_chr annotation file \n");
	}else{
		$files_parallel{$mta{$out_key . '->output'}} = "f";
		$pr_names{$mta{$out_key . '->output'}} = $out_key;
	}
	if ( $mta{$out_key . '->exit_status'} != 0 ) {
		usage("dbnsfp $this_chr annotation has non-zero exit status\n");
	}
	
	$dbnsfp_db =  ( exists $mta{$out_key . '->database'} ) ? $mta{$out_key . '->database'} : -1;
	if ( $dbnsfp_db == -1 ) {
		&usage("pipeline_meta: Missing dbnsfp database for $this_chr (annovar generated) file \n");
	}	
}


#my $gnADgn_db = "NA";
##gnad chrwise
#foreach my $this_chr (@pipeline_valid_chrs){
#	#gerpwgs:12->output
#	$this_chr =~ s/chr//;
#	my $out_key = "gnomad_genome:" .  $this_chr; # . "->output";
#	my $this_gnADgn_anno =  ( exists $mta{$out_key . '->output'} )   ? $mta{$out_key. '->output'} : -1;
#	if($this_gnADgn_anno == -1 ){
#		&usage2("pipeline_meta: Missing gnomAD_genome $this_chr annotation file \n");
#	}else{
#		$files_parallel{$mta{$out_key . '->output'}} = "f";
#		$pr_names{$mta{$out_key . '->output'}} = $out_key;
#	}
#	if ( $mta{$out_key . '->exit_status'} != 0 ) {
#		usage("gnomAD_genome $this_chr annotation has non-zero exit status\n");
#	}
#	
#	$gnADgn_db =  ( exists $mta{$out_key . '->database'} ) ? $mta{$out_key . '->database'} : -1;
#	if ( $gnADgn_db == -1 ) {
#		&usage("pipeline_meta: Missing gnomAD_genome database for $this_chr (annovar generated) file \n");
#	}	
#}

#gnad chrwise
foreach my $this_chr (@pipeline_valid_chrs){
	#gerpwgs:12->output
	$this_chr =~ s/chr//;
	my $out_key = "gnomad_genomeV30:" .  $this_chr; # . "->output";
	my $this_gnADgn_anno =  ( exists $mta{$out_key . '->output'} )   ? $mta{$out_key. '->output'} : -1;
	if($this_gnADgn_anno == -1  ){
		&usage2("pipeline_meta: Missing/empty gnomAD_genome $this_chr annotation file \n");
	}else{
		$files_parallel{$mta{$out_key . '->output'}} = "f";
		$pr_names{$mta{$out_key . '->output'}} = $out_key;
	}
	if ( $mta{$out_key . '->exit_status'} != 0 ) {
		usage("gnomAD_genome $this_chr annotation has non-zero exit status\n");
	}
	
	#$gnADgn_db =  ( exists $mta{$out_key . '->database'} ) ? $mta{$out_key . '->database'} : -1;
}


#gnad211 chrwise
foreach my $this_chr (@pipeline_valid_chrs){
	#gerpwgs:12->output
	$this_chr =~ s/chr//;
	my $out_key = "gnomad_genomeV211:" .  $this_chr; # . "->output";
	my $this_gnADgn211_anno =  ( exists $mta{$out_key . '->output'} )   ? $mta{$out_key. '->output'} : -1;
	if($this_gnADgn211_anno == -1  ){
		&usage2("pipeline_meta: Missing gnomAD211_genome $this_chr annotation file \n");
	}else{
		$files_parallel{$mta{$out_key . '->output'}} = "f";
		$pr_names{$mta{$out_key . '->output'}} = $out_key;
	}
	if ( $mta{$out_key . '->exit_status'} != 0 ) {
		usage("gnomAD_genome211 $this_chr annotation has non-zero exit status\n");
	}
	
	#$gnADgn_db =  ( exists $mta{$out_key . '->database'} ) ? $mta{$out_key . '->database'} : -1;
}



#my $gnADgn_db =  ( exists $mta{'gnomad_genome->database'} ) ? $mta{'gnomad_genome->database'} : -1;
#if ( $gnADgn_db == -1 ) {
#	&usage("pipeline_meta: Missing gnomeAD_genome(annovar generated) file \n");
#}

#PhastCon->output
#&usage2("Missing phastcons_placental annotation file \n")       unless ( -s $options{h} );
warn "**** Phastcon = $mta{'phastcon->output'} \n";
my $phastcon_anno =
  ( exists $mta{'phastcon->output'} )
  ? $mta{'phastcon->output'}
  : -1;    #  &usage2("Missing phastcons_placental annotation file \n")  ;
if ( $phastcon_anno == -1   ) {
	&usage2("pipeline_meta: Missing phastcons_placental annotation file ");
}else{
	$files_parallel{$mta{'phastcon->output'}} = "r";
	$pr_names{$mta{'phastcon->output'}} = "phastcon";
}
if ( $mta{'phastcon->exit_status'} != 0 ) {
	usage("phastcon annotation has non-zero exit status\n");
}


#phylop chrwise, added Nov 2,2018
foreach my $this_chr (@pipeline_valid_chrs){
	#gerpwgs:12->output
	$this_chr =~ s/chr//;
	my $out_key = "phylop:" .  $this_chr; # . "->output";
	my $this_phylop_anno =  ( exists $mta{$out_key . '->output'} )   ? $mta{$out_key. '->output'} : -1;
	if($this_phylop_anno == -1 ){
		&usage2("pipeline_meta: Missing Phylop $this_chr annotation file \n");
	}else{
		$files_parallel{$mta{$out_key . '->output'}} = "r";
		$pr_names{$mta{$out_key . '->output'}} = $out_key;
	}
	if ( $mta{$out_key . '->exit_status'} != 0 ) {
		usage("phylop $this_chr annotation has non-zero exit status\n");
	}
}



#phylop100 chrwise, added Nov 2,2018

#print scalar(@pipeline_valid_chrs) , " len \n";
foreach my $this_chr (@pipeline_valid_chrs){
	#gerpwgs:12->output
	$this_chr =~ s/chr//;
#	print "chr = $this_chr\n ";
	my $out_key = "phylop100:" .  $this_chr; # . "->output";
	my $this_phylop_anno =  ( exists $mta{$out_key . '->output'} )   ? $mta{$out_key. '->output'} : -1;
	if($this_phylop_anno == -1 ){
		&usage2("pipeline_meta: Missing Phylop100 $this_chr annotation file ");
	}else{
		$files_parallel{$mta{$out_key . '->output'}} = "r";
		$pr_names{$mta{$out_key . '->output'}} = $out_key;
	}
	if ( $mta{$out_key . '->exit_status'} != 0 ) {
		usage("phylop100 $this_chr annotation has non-zero exit status\n");
	}
}

#CG internal freq chrwise, added Mar 24,2020 (mssng only)

#print scalar(@pipeline_valid_chrs) , " len \n";

#pipeline->post_process
my $run_pp =
  ( exists $mta{'pipeline->post_process'} ) ? $mta{'pipeline->post_process'} : -1;
if ( $run_pp == -1 ) {
	&usage("pipeline_meta: Missing pipeline->post_process info ");
}

if($run_pp == 5){
	foreach my $this_chr (@pipeline_valid_chrs){
		$this_chr =~ s/chr//;
		my $out_key = "cgInt:" .  $this_chr; # . "->output";
		my $this_cgint_anno =  ( exists $mta{$out_key . '->output'} )   ? $mta{$out_key. '->output'} : -1;
		if($this_cgint_anno == -1 ){
			&usage2("pipeline_meta: Missing cgInt $this_chr annotation file ");
		}else{
			$files_parallel{$mta{$out_key . '->output'}} = "ifr";
			$pr_names{$mta{$out_key . '->output'}} = $out_key;
		}
		if ( $mta{$out_key . '->exit_status'} != 0 ) {
			usage("cgInt $this_chr annotation has non-zero exit status\n");
		}
	}
	
	foreach my $this_chr (@pipeline_valid_chrs){
		$this_chr =~ s/chr//;
		my $out_key = "illInt:" .  $this_chr; # . "->output";
		my $this_cgint_anno =  ( exists $mta{$out_key . '->output'} )   ? $mta{$out_key. '->output'} : -1;
		if($this_cgint_anno == -1 ){
			&usage2("pipeline_meta: Missing illInt $this_chr annotation file ");
		}else{
			$files_parallel{$mta{$out_key . '->output'}} = "ifr";
			$pr_names{$mta{$out_key . '->output'}} = $out_key;
		}
		if ( $mta{$out_key . '->exit_status'} != 0 ) {
			usage("illInt $this_chr annotation has non-zero exit status\n");
		}
	}
	
}

##dbnsfp->output  , added Jan 26, 2016
#my $dbnsfp_anno =
#  ( exists $mta{'dbnsfp->output'} )
#  ? $mta{'dbnsfp->output'}
#  : -1;    #  &usage2("Missing polyphen annotation file \n")  ;
#if ( $dbnsfp_anno == -1 ) {
#	&usage2("pipeline_meta: Missing dbnsfp annotation file ");
#}else{
#	$files_parallel{$mta{'dbnsfp->output'}} = "f";
#	$pr_names{$mta{'dbnsfp->output'}} = "dbnsfp";
#}
#if ( $mta{'dbnsfp->exit_status'} != 0 ) {
#	usage("dbnsfp annotation has non-zero exit status\n");
#}
#my $dbnsfp_db =  ( exists $mta{'dbnsfp->database'} ) ? $mta{'dbnsfp->database'} : -1;
#if ( $dbnsfp_db == -1 ) {
#	&usage("pipeline_meta: Missing dbnsfp database(annovar generated) file \n");
#}


##primateAI->output  , added Sep 10,2018
#my $prAI_anno =
#  ( exists $mta{'primateAI->output'} )
#  ? $mta{'primateAI->output'}
#  : -1;    #  &usage2("Missing polyphen annotation file \n")  ;
#if ( $prAI_anno == -1 ) {
#	&usage2("pipeline_meta: Missing  primateAI annotation file \n");
#}else{
#	$files_parallel{$mta{'primateAI->output'}} = "f";
#	$pr_names{$mta{'primateAI->output'}} = "primateAI";
#}
#if ( $mta{'primateAI->exit_status'} != 0 ) {
#	usage("primateAI annotation has non-zero exit status\n");
#}
#my $prAI_db =  ( exists $mta{'primateAI->database'} ) ? $mta{'primateAI->database'} : -1;
#if ( $prAI_db == -1 ) {
#	&usage("pipeline_meta: Missing $prAI_db database(?) file \n");
#}


#dbscsnv->output  , added Jan 26, 2016
my $dbscsnv_anno =
  ( exists $mta{'dbscsnv->output'} )
  ? $mta{'dbscsnv->output'}
  : -1;    #  &usage2("Missing polyphen annotation file \n")  ;
if ( $dbscsnv_anno == -1 ) {
	&usage2("pipeline_meta: Missing dbscsnv annotation file ");
}else{
	$files_parallel{$mta{'dbscsnv->output'}} = "f";
	$pr_names{$mta{'dbscsnv->output'}} = "dbscsnv";
}
if ( $mta{'dbscsnv->exit_status'} != 0 ) {
	usage("dbscsnv annotation has non-zero exit status\n");
}
my $dbscsnv_db =  ( exists $mta{'dbscsnv->database'} ) ? $mta{'dbscsnv->database'} : -1;
if ( $dbscsnv_db == -1 ) {
	&usage("pipeline_meta: Missing dbscsnv database(annovar generated) file \n");
}

#$annovar_filter_module{'spliceAIsnv'}{db} =  $filter_spliceAIsnv; 
my $spliceAIsnv_anno =
  ( exists $mta{'spliceAIsnv->output'} )
  ? $mta{'spliceAIsnv->output'}
  : -1;    #  &usage2("Missing polyphen annotation file \n")  ;
if ( $spliceAIsnv_anno == -1 ) {
	&usage2("pipeline_meta: Missing spliceAIsnv annotation file ");
}else{
	$files_parallel{$mta{'spliceAIsnv->output'}} = "f";
	$pr_names{$mta{'spliceAIsnv->output'}} = "spliceAIsnv";
}
if ( $mta{'spliceAIsnv->output'} != 0 ) {
	usage("spliceAIsnv annotation has non-zero exit status\n");
}
my $dbscsnv_db =  ( exists $mta{'dbscsnv->database'} ) ? $mta{'dbscsnv->database'} : -1;
if ( $dbscsnv_db == -1 ) {
	&usage("pipeline_meta: Missing dbscsnv database(annovar generated) file \n");
}

#pipeline->input ****** This need to be fixed
#&usage("Please specify ORIGINAL input  file \n")            unless ( -s $options{m} );
my $orgi_input =
  ( exists $mta{'pipeline->input_to_go'} )
  ? $mta{'pipeline->input_to_go'}
  : -1;    # &usage("Please specify ORIGINAL input  file \n")  ;
if ( $orgi_input == -1 ) {
	&usage("pipeline_meta: Please specify  input file  \n");
}else{
	$files_parallel{$orgi_input} = "i";
	$pr_names{$orgi_input} = "Input";
}

#dbsnpR_merge->output
#&usage2("Missing dnsnp-region annotation file \n")    unless ( -s $options{r} );

foreach my $this_chr (@pipeline_valid_chrs){
	#gerpwgs:12->output
	$this_chr =~ s/chr//;
	my $out_key = "dbsnpR:" .  $this_chr; # . "->output";
	my $this_phylop_anno =  ( exists $mta{$out_key . '->output'} )   ? $mta{$out_key. '->output'} : -1;
	if($this_phylop_anno == -1 ){
		&usage2("pipeline_meta: Missing dbsnpR $this_chr annotation file ");
	}else{
		$files_parallel{$mta{$out_key . '->output'}} = "r";
		$pr_names{$mta{$out_key . '->output'}} = $out_key;
	}
	if ( $mta{$out_key . '->exit_status'} != 0 ) {
		usage("dbsnpR $this_chr annotation has non-zero exit status\n");
	}
}

#foreach my $this_chr (@pipeline_valid_chrs){
#	#gerpwgs:12->output
#	$this_chr =~ s/chr//;
#	my $out_key = "spliceAIsnv:" .  $this_chr; # . "->output";
#	my $this_splAIsnv_anno =  ( exists $mta{$out_key . '->output'} )   ? $mta{$out_key. '->output'} : -1;
#	if($this_splAIsnv_anno == -1 ){
#		&usage2("pipeline_meta: Missing spliceAIsnv $this_chr annotation file ");
#	}else{
#		$files_parallel{$mta{$out_key . '->output'}} = "f";
#		$pr_names{$mta{$out_key . '->output'}} = $out_key;
#	}
#	if ( $mta{$out_key . '->exit_status'} != 0 ) {
#		usage("spliceAIsnv $this_chr annotation has non-zero exit status\n");
#	}
#}



#****

#my $dbsnpR_anno =
#  ( exists $mta{'dbsnpR_merge->output'} )
#  ? $mta{'dbsnpR_merge->output'}
#  : -1;    # &usage2("Missing dnsnp-region annotation file \n")  ;
#if ( $dbsnpR_anno == -1 ) {
#	&usage2("pipeline_meta: Missing dnsnp-region annotation file \n");
#}else{
#	$files_parallel{$mta{'dbsnpR_merge->output'}} = "r";
#	$pr_names{$mta{'dbsnpR_merge->output'}} = "dbsnpR";
#}
#if ( $mta{'dbsnpR_merge->exit_status'} != 0 ) {
#	usage("dbsnpR_merge annotation has non-zero exit status\n");
#}


#dbsnpW_merge->output
#&usage2("Missing dbsnp-window annotation file \n")    unless ( -s $options{w} );

foreach my $this_chr (@pipeline_valid_chrs){
	#gerpwgs:12->output
	$this_chr =~ s/chr//;
	my $out_key = "dbsnpW:" .  $this_chr; # . "->output";
	my $this_phylop_anno =  ( exists $mta{$out_key . '->output'} )   ? $mta{$out_key. '->output'} : -1;
	if($this_phylop_anno == -1 ){
		&usage2("pipeline_meta: Missing dbsnpW $this_chr annotation file \n");
	}else{
		$files_parallel{$mta{$out_key . '->output'}} = "r";
		$pr_names{$mta{$out_key . '->output'}} = $out_key;
	}
	if ( $mta{$out_key . '->exit_status'} != 0 ) {
		usage("dbsnpW $this_chr annotation has non-zero exit status\n");
	}
}


#gene_info->output
#&usage("Please specify entrez-gene map file \n")      unless ( -s $options{z} );
my $ent_gene_map =
  ( exists $mta{'gene_info->output'} )
  ? $mta{'gene_info->output'}
  : -1;    # &usage("Please specify entrez-gene map file \n")  ;
if ( $ent_gene_map == -1 ) {
	&usage("pipeline_meta: Please specify entrez-gene map file  \n");
}

#&usage("Please specify refflat.txt  file \n")         unless ( -s $options{l} );
my $gene_refFlat =
  ( exists $mta{'gene_info->refFlat'} )
  ? $mta{'gene_info->refFlat'}
  : -1;    # &usage("Please specify entrez-gene map file \n")  ;
if ( $gene_refFlat == -1 ) {
	&usage("pipeline_meta: Please specify refseq refFlat file  \n");
}

#Pfam->output
#&usage2("missing pfam annotation file \n")            unless ( -s $options{p} );
my $pfam_anno =
  ( exists $mta{'pfam->output'} )
  ? $mta{'pfam->output'}
  : -1;    # &usage2("missing pfam annotation file \n")  ;
if ( $pfam_anno == -1 ) {
	&usage2("pipeline_meta: Missing pfam annotation file ");
}else{
	$files_parallel{$mta{'pfam->output'}} = "r";
	$pr_names{$mta{'pfam->output'}} = "pfam";
}
if ( $mta{'pfam->exit_status'} != 0 ) {
	usage("pfam annotation has non-zero exit status\n");
}

my $gerp_elem_anno =
  ( exists $mta{'gerp_elem->output'} )
  ? $mta{'gerp_elem->output'}
  : -1;    # &usage2("missing pfam annotation file \n")  ;
if ( $gerp_elem_anno == -1 ) {
	&usage2("pipeline_meta: Missing gerp_elem annotation file ");
}else{
	$files_parallel{$mta{'gerp_elem->output'}} = "r";
	$pr_names{$mta{'gerp_elem->output'}} = "gerp_elem";
}
if ( $mta{'gerp_elem->exit_status'} != 0 ) {
	usage("gerp_elem annotation has non-zero exit status\n");
}


#Clinvar->output
#&usage2("Missing Clinvar  annotation file \n") unless ( -s $options{g} );
my $clinvar_anno =
  ( exists $mta{'clinvar->output'} ) ? $mta{'clinvar->output'} : -1;
if ( $clinvar_anno == -1 ) {
	&usage2("pipeline_meta: Missing Clinvar annotation file ");
}else{
	$files_parallel{$mta{'clinvar->output'}} = "f";
	$pr_names{$mta{'clinvar->output'}} = "clinvar";
}
if ( $mta{'clinvar->exit_status'} != 0 ) {
	usage("Clinvar annotation has non-zero exit status\n");
}

##added June 15,2017
#my $clinvar_annovar_anno =
#  ( exists $mta{'Clinvar_annovar->output'} ) ? $mta{'Clinvar_annovar->output'} : -1;
#if ( $clinvar_annovar_anno == -1 ) {
#	&usage2("pipeline_meta: Missing Clinvar(annovar) annotation file \n");
#}else{
#	$files_parallel{$mta{'Clinvar_annovar->output'}} = "f";
#	$pr_names{$mta{'Clinvar_annovar->output'}} = "Clinvar_annovar";
#}
#if ( $mta{'Clinvar_annovar->exit_status'} != 0 ) {
#	usage("Clinvar(annovar) annotation has non-zero exit status\n");
#}
#
#my $intervar_anno =
#  ( exists $mta{'intervar->output'} ) ? $mta{'intervar->output'} : -1;
#if ( $intervar_anno == -1 ) {
#	&usage2("pipeline_meta: Missing intervar annotation file \n");
#}else{
#	$files_parallel{$mta{'intervar->output'}} = "f";
#	$pr_names{$mta{'intervar->output'}} = "intervar";
#}
#if ( $mta{'intervar->exit_status'} != 0 ) {
#	usage("Intervar annotation has non-zero exit status\n");
#}


##added July 29, 2016 for clinvar summary
#my $clinvar_summary = ( exists $mta{'clinvar_summary->file'} ) ? $mta{'clinvar_summary->file'} : -1;
#if ( $clinvar_summary == -1 ) {
#	&usage2("pipeline_meta: Missing Clinvar summary file ");
#}

#CADD merge->output
#&usage2("Missing CADD  annotation file \n") unless ( -s $options{M} );
my $cadd_anno =
  ( exists $mta{'cadd->output'} ) ? $mta{'cadd->output'} : -1;
if ( $cadd_anno == -1 ) {
	&usage2("pipeline_meta: Missing CADD annotation file ");
}else{
	$files_parallel{$mta{'cadd->output'}} = "f";
	$pr_names{$mta{'cadd->output'}} = "cadd";
}
if ( $mta{'cadd->exit_status'} != 0 ) {
	usage("CADD annotation has non-zero exit status\n");
}

#spidex->output
#&usage2("Missing SPIDEX  annotation file \n") unless ( -s $options{S} );
my $spidex_anno =
  ( exists $mta{'spidex->output'} ) ? $mta{'spidex->output'} : -1;
if ( $spidex_anno == -1 ) {
	&usage2("pipeline_meta: Missing SPIDEX annotation file ");
}else{
	$files_parallel{$mta{'spidex->output'}} = "f";
	$pr_names{$mta{'spidex->output'}} = "SPIDEX";
}
if ( $mta{'spidex->exit_status'} != 0 ) {
	usage("spidex annotation has non-zero exit status\n");
}

#SegDup->output
#&usage2("Missing segdup annotation file \n")        unless ( -s $options{u} );
my $segdup_anno =
  ( exists $mta{'segdup->output'} ) ? $mta{'segdup->output'} : -1;
if ( $segdup_anno == -1 ) {
	&usage2("pipeline_meta: Missing SegDup annotation file ");
}else{
	$files_parallel{$mta{'segdup->output'}} = "r";
	$pr_names{$mta{'segdup->output'}} = "segdup";
}
if ( $mta{'segdup->exit_status'} != 0 ) {
	usage("SegDup annotation has non-zero exit status\n");
}

#Cosmic->output
#&usage2("Missing  COSMIC  annotation file \n")      unless ( -s $options{H} );
my $cosmic_anno =
  ( exists $mta{'cosmic->output'} ) ? $mta{'cosmic->output'} : -1;
if ( $cosmic_anno == -1 ) {
	&usage2("pipeline_meta: Missing cosmic annotation file ");
}else{
	$files_parallel{$mta{'cosmic->output'}} = "f";
	$pr_names{$mta{'cosmic->output'}} = "cosmic";
}
if ( $mta{'cosmic->exit_status'} != 0 ) {
	usage("cosmic annotation has non-zero exit status\n");
}

#HGMD pro->output
#&usage2("Missing  HGMD pro annotation file \n")      unless ( -s $options{k} );
my $hgmdpro_anno =
  ( exists $mta{'hgmd_pro->output'} ) ? $mta{'hgmd_pro->output'} : -1;
if ( $hgmdpro_anno == -1 ) {
	&usage2("pipeline_meta: Missing HGMD-professional annotation file ");
}else{
	$files_parallel{$mta{'hgmd_pro->output'}} = "f";
	$pr_names{$mta{'hgmd_pro->output'}} = "hgmd_pro";
}
if ( $mta{'hgmd_pro->exit_status'} != 0 ) {
	usage("HGMD pro annotation has non-zero exit status\n");
}


#gene_extra, added Sep25,2019
#&usage2("Missing  HGMD pro annotation file \n")      unless ( -s $options{k} );
#my $geneEx_out =
#  ( exists $mta{'gene_extra->output'} ) ? $mta{'gene_extra->output'} : -1;
#if ( $geneEx_out == -1 ) {
#	&usage2("pipeline_meta: Missing gene_extra output file ");
#}else{
##	$files_parallel{$mta{'hgmd_pro->output'}} = "f";
#	$pr_names{$mta{'gene_extra->output'}} = "gene_extra";
#}
#if ( $mta{'gene_extra->exit_status'} != 0 ) {
#	usage("gene_extra has non-zero exit status\n");
#}
#if (not exists $mta{'gene_extra->exit_status'}  ) {
#	usage("Critical pipeline error! could not find gene_extra modules... (number of cores > 16 ?) \n");
#}


#gene_extra, updated Jan10,2020
my @gene_extra_stats = ();
foreach my $this_chr (@pipeline_valid_chrsAll){
	my $out_key = "gene_extra:" .  $this_chr; # . "->output";
	#warn "***GeneEx: $out_key\n";
	my $this_geneEx_out =  ( exists $mta{$out_key . '->output'} )   ? $mta{$out_key. '->output'} : -1;
	if($this_geneEx_out == -1 ){
		&usage("pipeline_meta: Missing gene_extra $this_chr annotation file \n");
	}else{
		$files_parallel{$mta{$out_key . '->output'}} = "gex";
		$pr_names{$mta{$out_key . '->output'}} = $out_key;
	}
	if ( $mta{$out_key . '->exit_status'} != 0 ) {
		usage("gene_extra $this_chr annotation has non-zero exit status\n");
	}
	if (not exists $mta{$out_key . '->exit_status'}  ) {
		usage("Critical pipeline error! could not find gene_extra : $out_key module... (too many cores..?) \n");
	}
	if (exists $mta{$out_key . '->stats'}  ) {
		push(@gene_extra_stats , $mta{$out_key . '->stats'});
	}else{
		usage("Critical pipeline error! could not find gene_extra : $out_key stats file... (using old version of module script ? ) \n");
	}

	
}


#**********  Followings need to be added to %process_status
#&usage("Please specify CFG file \n")                unless ( -s $options{s} );
my $pipeline_config =
  ( exists $mta{'pipeline->config_file'} ) ? $mta{'pipeline->config_file'} : -1;
if ( $pipeline_config == -1 ) {
	&usage("pipeline_meta: Missing config file \n");
}

#&usage("Please specify morbidmap file \n")                unless ( -s $options{I} );
my $morbidmap_file =
  ( exists $mta{'Phenotype->morbidmap'} ) ? $mta{'Phenotype->morbidmap'} : -1;
if ( $morbidmap_file == -1 ) {
	&usage("pipeline_meta: OMIM morbidmap file is missing ");
}

#&usage("Please specify HPO file \n")                unless ( -s $options{J} );
my $hpo_file = ( exists $mta{'Phenotype->HPO'} ) ? $mta{'Phenotype->HPO'} : -1;
if ( $hpo_file == -1 ) {
	&usage("pipeline_meta: HPO file is missing ");
}

#&usage("Please specify MPO file \n")                unless ( -s $options{K} );
my $mpo_file = ( exists $mta{'Phenotype->MPO'} ) ? $mta{'Phenotype->MPO'} : -1;
if ( $mpo_file == -1 ) {
	&usage("pipeline_meta: MPO file is missing ");
}

#&usage("Please specify CGD file \n")                unless ( -s $options{i} );
my $cgd_file = ( exists $mta{'Phenotype->CGD'} ) ? $mta{'Phenotype->CGD'} : -1;
if ( $cgd_file == -1 ) {
	&usage("pipeline_meta: CGD file is missing ");
}

#&usage("Please specify GI file \n")                unless ( -s $options{U} );
my $gi_file = ( exists $mta{'Phenotype->GI'} ) ? $mta{'Phenotype->GI'} : -1;
if ( $gi_file == -1 ) {
	&usage2("pipeline_meta: GI file is missing ");
}

#&usage("Please specify HI file \n")                unless ( -s $options{V} );
my $hi_file = ( exists $mta{'Phenotype->HI'} ) ? $mta{'Phenotype->HI'} : -1;
if ( $hi_file == -1 ) {
	&usage2("pipeline_meta: HI file is missing ");
}

my $exmat_file = ( exists $mta{'Phenotype->ExAc_Matrics'} ) ? $mta{'Phenotype->ExAc_Matrics'} : -1;
if ( $exmat_file == -1 ) {
	&usage("pipeline_meta: ExAc_Matrics file is missing ");
}

#added July 26,2017
my $acmg_file = ( exists $mta{'Phenotype->ACMG'} ) ? $mta{'Phenotype->ACMG'} : -1;
if ( $acmg_file == -1 ) {
	&usage("pipeline_meta: ACMG file is missing ");
}

#*******************************************


my $rpt_anno =
  ( exists $mta{'repeat->output'} ) ? $mta{'repeat->output'} : -1;
if ( $rpt_anno == -1 ) {
	&usage2("pipeline_meta: Missing repeat annotation file ");
}else{
	$files_parallel{$mta{'repeat->output'}} = "r";
	$pr_names{$mta{'repeat->output'}} = "repeat";
}
if ( $mta{'repeat->exit_status'} != 0 ) {
	usage("repeat annotation has non-zero exit status\n");
}

##exon boundary , added June 12,2018
#my $exb_anno =  ( exists $mta{'exome_boundary_dist->output'} ) ? $mta{'exome_boundary_dist->output'} : -1;
#if ( $exb_anno == -1 ) {
#	&usage2("pipeline_meta: Missing exome_boundary_dist annotation file \n");
#}else{
#	$files_parallel{$mta{'exome_boundary_dist->output'}} = "b";
#	$pr_names{$mta{'exome_boundary_dist->output'}} = "exome_boundary_dist";
#}
#if ( $mta{'exome_boundary_dist->exit_status'} != 0 ) {
#	usage("exome_boundary_dist annotation has non-zero exit status\n");
#}



##Gene->var_output
#my $exb_refgene_var =
#  ( exists $mta{'exome_boundary_dist->var_output'} )
#  ? $mta{'exome_boundary_dist->var_output'}
#  : -1;    # &usage2("Missing gene-variant annotation file \n")  ;
#if ( $exb_refgene_var == -1 ) {
#	&usage2("pipeline_meta: Missing extended gene-variant annotation file  \n");
#}
#
#if ( $mta{'exome_boundary_dist->exit_status'} != 0 ) {
#	usage("Extended Gene annotation has non-zero exit status\n");
#}else{
#	$files_parallel{$mta{'exome_boundary_dist->var_output'}} = "b1";
#	$pr_names{$mta{'exome_boundary_dist->var_output'}} = "exome_boundary_dist_var";
#}
#
## unless ( -s $options{x} );
#my $exb_refgene_vex =
#  ( exists $mta{'exome_boundary_dist->exonic_output'} )
#  ? $mta{'exome_boundary_dist->exonic_output'}
#  : -1;    # &usage2("Missing gene-exonic variant annotation file \n")  ;
#if ( $exb_refgene_vex == -1 ) {
#	&usage2("pipeline_meta: Missing extended gene-exonic variant annotation file  \n");
#}else{
#	$files_parallel{$mta{'exome_boundary_dist->exonic_output'}} = "b2";
#	$pr_names{$mta{'exome_boundary_dist->exonic_output'}} = "exome_boundary_dist_var-exonic";
#}


#pipeline->run_type
my $run_type =
  ( exists $mta{'pipeline->run_type'} ) ? $mta{'pipeline->run_type'} : -1;
if ( $run_type == -1 ) {
	&usage("pipeline_meta: Missing run_type info \n");
}

#pipeline->input_count
my $input_var_count =
  ( exists $mta{'pipeline->input_count'} ) ? $mta{'pipeline->input_count'} : -1;
if ( $input_var_count == -1 ) {
	&usage("pipeline_meta: Missing input_count info \n");
}



#$process_status{'pipeline'}{impact}
my $calc_impact =
  ( exists $mta{'pipeline->impact'} ) ? $mta{'pipeline->impact'} : -1;
if ( $calc_impact == -1 ) {
	&usage("pipeline_meta: Missing pipeline->impact info ");
}

warn "post_process = $run_pp\n";
#extra->count=4999
my $extra_count =
  ( exists $mta{'extra->count'} ) ? $mta{'extra->count'} : -1;
if ( $extra_count == -1 ) {
	&usage("pipeline_meta: Missing extra->count info ");
}

#determine original run type ( added Aug 2,2019) 
if(exists $mta{'vcf_pp->pid'}){
	warn "**Detected vcf type input: \n";
	#vcf_pp->exit_status=0
	#vcf_pp->output
	if ( $mta{'vcf_pp->exit_status'} != 0 ) {
		&usage("Detected vcf type input, but vcf post-processing module (vcf_pp) has non-zero exit status\n");
	}else{
		$annovar_in = $mta{'vcf_pp->output'};
	}
}

#warn "Annovar-in = $annovar_in\n";

my $unit_stat = ($extra_count / 20);
my $invalid_file = $refgene_var;
$invalid_file =~ s/variant_function/invalid_input/;
my $invalid_out =
  $pipeline_tmp_dir . "/invalid_entries.tsv";
my $stat_file =
  $pipeline_tmp_dir . "/stats.meta";

#warn "$invalid_file  -- $invalid_out  \n";
#open( INVO, ">$invalid_out" ) or die "$! $invalid_out";
my $fh_invo = IO::File->new($invalid_out , q{>}) or die "$! $invalid_out\n";
#open( VST,  ">$stat_file" )   or die "$! $stat_file";
my $fh_vst = IO::File->new($stat_file , q{>}) or die "$! $stat_file\n";

my $fh_sam = IO::File->new($annovar_in , q{<}) or usage("annovar input , $annovar_in , is not readable-- $!");
  
my $fh_urfl = IO::File->new($uc_knowngene , q{<}) or usage("UCSC knowngene , $uc_knowngene , is not readable -- $!");

my $fh_gvar = IO::File->new($refgene_var , q{<}) or usage("gene-var annotation , $refgene_var , is not readable ( $! )");
  
#open( GEX, $refgene_vex )
#  or usage2(
#	"gene-exonic variant annotation file, $refgene_var , is not readable ( $! )"
#  );
my $fh_gex = IO::File->new($refgene_vex , q{<}) or usage("gene-exonic variant annotation file, $refgene_var , is not readable ( $! )");

#ucsc gene annotations
my $fh_uvar = IO::File->new($ucgene_var , q{<}) or usage("UCSC gene-var file , $ucgene_var , is not readable -- $!");
my $fh_uex = IO::File->new($ucgene_vex , q{<}) or usage("UCSC gene-exonic var annotation, $ucgene_vex , is not readable -- $!");


my $fh_gen = IO::File->new($ent_gene_map , q{<}) or usage("gene info file  , $ent_gene_map , is not readable ( $! )");
my $fh_rfl = IO::File->new($gene_refFlat , q{<}) or usage("gene refFlat file  , $gene_refFlat , is not readable ( $! )");

my $fh_cfg = IO::File->new($pipeline_config , q{<}) or usage("Config file , $pipeline_config , is not readable ( $! )");
 
my $fh_omi = IO::File->new($morbidmap_file, q{<}) or usage("Morbidmap file , $morbidmap_file , is not readable ( $! )");
my $fh_hpo = IO::File->new($hpo_file, q{<}) or usage("HPO file , $hpo_file , is not readable ( $! )");
my $fh_mpo = IO::File->new($mpo_file, q{<}) or usage("MPO file , $mpo_file , is not readable ( $! )");
my $fh_cgd = IO::File->new($cgd_file, q{<}) or usage("CGD file , $cgd_file , is not readable ( $! )");
my $fh_gi = IO::File->new($gi_file, q{<}) or usage2("GI file , $gi_file , is not readable ( $! )");
my $fh_hi = IO::File->new($hi_file, q{<}) or usage2("HI file , $hi_file , is not readable ( $! )");
my $fh_xmt = IO::File->new($exmat_file, q{<}) or usage("exac mat file , $exmat_file , is not readable ( $! )");
#my $fh_clvs = IO::File->new($clinvar_summary, q{<}) or usage2("clinvar summary , is not readable ( $! )");
my $fh_acmg = IO::File->new($acmg_file, q{<}) or usage("ACMG file , is not readable ( $! )");
#my $fh_gnEx = IO::File->new($geneEx_out, q{<}) or usage("$geneEx_out (gene_extra) , is not readable ( $! )");


#my $fh_hg38Ambg = IO::File->new($db_hg39Ambg, q{<}) or usage("hg38_ambiguous_sites , $db_hg39Ambg  is not readable ( $! )");

#**********************************************

 my $this_time = localtime;
 warn "\[$this_time\] Intial checking done \n";	

my %gvar = ();
my %gex  = ();

#ucsc genes
my %uvar = ();
my %uex  = ();
my %urfl = ();

my %cg      = ();
my %wld     = ();
my %ckb     = ();
my %cgf     = ();
my %snp     = ();
my %snpc    = ();
my %cln     = ();
my %thg     = ();
my %esp     = ();
my %phs     = ();
my %phsl    = ();
my %phy     = ();
my %phy_avg = ();

my %phypr     = ();
my %phypr_avg = ();

my %pol        = ();
my %sft        = ();
my %seg        = ();
my %mvar       = ();
my %dbr1       = ();
my %dbr2       = ();
my %dbrc       = ();
my %entz_id      = ();
my %uc_genes   = ();
my %gene_desc  = ();
my %gene_types = ();
my %gene_mim   = ();
my %mim_desc   = ();
my %gene_mpo   = ();
my %gene_hpo   = ();
my %gene_cgd   = ();
my %gene_gi    = ();
my %gene_hi    = ();
my %gene_exmat    = ();
my %cd_lens    = ();

my %uc_cd_lens        = ();
my %tr_counts         = ();
my %uc_tr_counts      = ();
my %pfam              = ();
my %dnf               = ();
my %dnr               = ();
my %dnw               = ();
#my %seg               = ();
my %pyr               = ();
my %hgm               = ();
my %cfg               = ();
my %nc1               = ();
my %nc2               = ();
my %nc3               = ();
my %nc4               = ();
my %nc5               = ();
my %nc6               = ();
my %nc7               = ();
my %nc4w              = ();
my %nc5w              = ();
my %nc6w              = ();
my %nc7w              = ();
my %reps              = ();
my %cadd              = ();
my %spdx              = ();
my %hgm2              = ();
my %nh_exac           = ();
my %refflats          = ();
my %uc_refflats       = ();
my %exon_dist         = ();
my %exon_dist_uc      = ();
#my %mta               = ();
my %transcripts_entrz = ();
my %dbnsfp_anno = ();
my %dbscsnv_anno = ();
my %db_giov = ();
my %clin_summary = ();
my %acmg = ();
my %hg38_ambg = ();

our ( %tam, %tsn, %theu, %thaf, %thss, %esaa, %esea, %maa, %cosmic, %no_of_tr )
  = ();

my %stats_final_merge = ();  
my $mvar_header = "";
warn "Program version : $0\n";
while ( my $line = <$fh_cfg> ) {
	chomp($line);
	my $key = ( split( "=", $line ) )[0];
	$cfg{$key} = ( split( "=", $line ) )[1];
}
close $fh_cfg;
$this_time = localtime;
warn "\[$this_time\] : Done populating CFG hash\n";

#refFlat ; This is for computing "cds affected"
#my $stat_count = 0;
#while ( my $line = <$fh_rfl> ) {
#	chomp($line);
#	my $key     = ( split( "\t", $line ) )[1];
#	my $rfl_chr = ( split( "\t", $line ) )[2];
#	$rfl_chr =~ s/chr//;
#	my $key_reff = $key . ":" . $rfl_chr;
#	$refflats{$key_reff} = $line;
#
#	# for transcript count
#	my $key_gene = ( split( "\t", $line ) )[0];
#	my $key_gene_rf = $rfl_chr . ":" . ( split( "\t", $line ) )[0];
#	
#	if ( exists $tr_counts{$key_gene_rf} ) {
#		$tr_counts{$key_gene_rf} =
#		  $tr_counts{$key_gene_rf} . ":" . ( split( "\t", $line ) )[1];
#	}
#	else {
#		$tr_counts{$key_gene_rf} = ( split( "\t", $line ) )[1];
#	}
#	$stat_count++;
#
#}
#close $fh_rfl;
#warn "Done populatiing RefFlat hash\n";
#$stats_final_merge{counts}{'refFlat'} = $stat_count;

#1       LOC102725121
#2       NR_148357
#3       chr1
#4       +
#5       11868
#6       14362
#7       14362
#8       14362
#9       3
#10      11868,12612,13220,
#11      12227,12721,14362,

my $stat_count = 0;
my %gene_exons = ();
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
	my $exon_ends = ( split( "\t", $line ) )[10];
	if($key !~ m/^NR_/){
		#print "$key -- not NR , $key_gene_rf\n";
		$gene_exons{$key_gene_rf}{'exstart'} = ($gene_exons{$key_gene_rf}{'exstart'}) ? $gene_exons{$key_gene_rf}{'exstart'} . $exon_starts : $exon_starts;
		$gene_exons{$key_gene_rf}{'exend'} = ($gene_exons{$key_gene_rf}{'exend'}) ? $gene_exons{$key_gene_rf}{'exend'} . $exon_ends : $exon_ends;
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
warn "Done populatiing RefFlat hash\n";
$stats_final_merge{counts}{'refFlat'} = $stat_count;


$stat_count = 0;

#gene_extra, added Sep25,2019

#1       1:1647814:1647814:T:C
#2       CDK11A
#3       29
#4       CDK11A:NM_024011:exon5:c.A459G:p.E153E,CDK11A:NM_001313982:exon5:c.A429G:p.E143E,CDK11A:NM_033529:exon5:c.A429G:p.E143E,CDK11A:NM_001313896:exon5:c.A429G:p.E143E
#5       CDK11A:NM_024011:exon5:c.A459G:p.E153E,CDK11A:NM_001313982:exon5:c.A429G:p.E143E,CDK11A:NM_033529:exon5:c.A429G:p.E143E,CDK11A:NM_001313896:exon5:c.A429G:p.E143E
#6       synonymous SNV
#7       synonymous SNV

#my %gene_extra = ();
#while ( my $line = <$fh_gnEx> ) {
#	chomp($line);
#	my ($ky, $this_typs, $this_gn, $this_dist, $rfq , $pr_rfq, $ef_cat, $pr_ef_cat, $rlxed, $ef_appris) = split("\t" , $line);
#	$gene_extra{$ky}{$this_gn}{$this_typs}{'dist'} = $this_dist;
#	$gene_extra{$ky}{$this_gn}{$this_typs}{'rfq'} = $rfq;
#	$gene_extra{$ky}{$this_gn}{$this_typs}{'pr_rfq'} = $pr_rfq;
#	$gene_extra{$ky}{$this_gn}{$this_typs}{'eff'} = $ef_cat;
#	$gene_extra{$ky}{$this_gn}{$this_typs}{'pr_eff'} = $pr_ef_cat;
#	$gene_extra{$ky}{$this_gn}{$this_typs}{'relaxed'} = $rlxed;
#	$gene_extra{$ky}{$this_gn}{$this_typs}{'eff_appris'} = $ef_appris;	
#	
#}
#close $fh_gnEx;
#$this_time = localtime;
#warn "\[$this_time\] : Done populating gene_extra hash\n";



#ucsc knowngene ; above for USCS genes
while ( my $line = <$fh_urfl> ) {
	chomp($line);
	my $key     = ( split( "\t", $line ) )[0];
	my $rfl_chr = ( split( "\t", $line ) )[1];
	$rfl_chr =~ s/chr//;
	my $key_reff = $key . ":" . $rfl_chr;
	$uc_refflats{$key_reff} = $line;
	$stat_count++;
}
close $fh_urfl;
warn "Done populatiing UCSC knowngene hash\n";
$stats_final_merge{counts}{'ucsc-knowngene'} = $stat_count;

$stat_count = 0;

#invalid file .. to be fixed
my %ani = ();
while ( my $line = <GINV> ) {
	chomp($line);
	if ( $run_type eq "v" ) {
		print INVO $line, "\n";
	}
	my $locus_id = extract_locusid($line);
	$ani{$locus_id} = 1;

}
close GINV;

#gene-id mapping , modified on Jan 2015 using refLink and process_reflink_v4

while ( my $line = <$fh_gen> ) {
	chomp($line);
	my @fields    = split( "\t", $line );
	my $gen_tr_id = $fields[0];
	my $gsymbol   = $fields[1];
	my $ent_id    = $fields[2];
#	print "Gene = $gsymbol\n";
	$entz_id{$gsymbol}               = $ent_id;
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
$stats_final_merge{counts}{'gene-id-map'} = $stat_count;

$stat_count = 0;

#Giovanna generated this file, populating gene-related info from here , updated July 27, 2016
#1       Hs.egID
#2       gene.symbol
#3       gene.alias
#4       OMIM.ID
#5       OMIM.name
#6       info
#7       gene.description
#8       gene.type

while ( my $line = <$fh_omi> ) {
	chomp($line);
	my @fields = split("\t" , $line);
	my $egID = $fields[0];
	$db_giov{$egID}{'OMIM.ID'} = $fields[3];
	$db_giov{$egID}{'OMIM.name'} = $fields[4];
	$db_giov{$egID}{'gene.description'} = $fields[6];
	$db_giov{$egID}{'gene.type'} = $fields[7];
	$stat_count++;
} 
close $fh_omi;
$stats_final_merge{counts}{'omim-ext'} = $stat_count;
$stat_count = 0;

#added clinvar summary July 29, 2016

##GeneID Gene_symbol     var_miss_pathogenic     var_miss_benign var_lof_pathogenic      var_lof_benign
#10      NAT2    0       0       0       0
#100     ADA     24      0       1       0

#while ( my $line = <$fh_clvs> ) {
#	chomp($line);
#	my ($clg_id, $clg_name, $miss_pathogenic , $miss_benign , $lof_pathogenic , $lof_benign ) = split("\t" , $line);
#	$clin_summary{$clg_id}{'miss_pathogenic'} = $miss_pathogenic;
#	$clin_summary{$clg_id}{'miss_benign'} = $miss_benign;
#	$clin_summary{$clg_id}{'lof_pathogenic'} = $lof_pathogenic;
#	$clin_summary{$clg_id}{'lof_benign'} = $lof_benign;
#	$stat_count++;
#} 
#close $fh_clvs;
#$stats_final_merge{counts}{'Clinvar-summary'} = $stat_count;

$stat_count = 0;


#ACMG added July 2017
while ( my $line = <$fh_acmg> ) {
	chomp($line);
	my ($acmg_gid, $acmg_symbol, $acmg_pheno) = split("\t" , $line);
	$acmg{$acmg_gid} = $acmg_pheno;
} 
close $fh_acmg;

#MPO
while ( my $line = <$fh_mpo> ) {
	chomp($line);
	my $key   = ( split( "\t", $line ) )[0];
	my $value = ( split( "\t", $line ) )[1];
	$gene_mpo{$key} =
	  ( exists $gene_mpo{$key} ) ? ( $gene_mpo{$key} . ";" . $value ) : $value;
	$stats_final_merge{counts}{'MPO'} = exists ($stats_final_merge{counts}{'MPO'}) ? ($stats_final_merge{counts}{'MPO'} + 1): 0 ; 
}

close $fh_mpo;

#HPO
while ( my $line = <$fh_hpo> ) {
	chomp($line);
	my $key   = ( split( "\t", $line ) )[0];
	my $value = ( split( "\t", $line ) )[1];
	$gene_hpo{$key} =
	  ( exists $gene_hpo{$key} ) ? ( $gene_hpo{$key} . ";" . $value ) : $value;
	$stats_final_merge{counts}{'HPO'} = exists ($stats_final_merge{counts}{'HPO'}) ? ($stats_final_merge{counts}{'HPO'} + 1): 0 ; 
}

close $fh_hpo;

#CGD
while ( my $line = <$fh_cgd> ) {
	chomp($line);
	my $key = ( split( "\t", $line ) )[0];
	my $value = ( split( "\t", $line ) )[1] . "@" . ( split( "\t", $line ) )[2];
	$gene_cgd{$key} =
	  ( exists $gene_cgd{$key} ) ? ( $gene_cgd{$key} . ";" . $value ) : $value;
}

close $fh_cgd;

#hg38 ambgeous sites:
#12      132241094       132241094       C       G
#12      132241118       132241118       A       C
#12      132241124       132241124       A       C

#while(my $line = <$fh_hg38Ambg>){
#	chomp($line);
#	my @fields = split("\t" , $line);
#	my $chrPre = ($fields[0] eq "MT") ? "M" : $fields[0];
#	$hg38_ambg{$fields[0] . ":" . $fields[1] . ":" . $fields[2]} = $fields[3];
#	$hg38_ambg{"chr" . $chrPre . ":" . $fields[1] . ":" . $fields[2]} = $fields[3];
#}
#close $fh_hg38Ambg;
#warn "Loaded ambiguous hg38 sites..\n";


##gi
#while ( my $line = <$fh_gi> ) {
#	chomp($line);
#	my $key   = ( split( "\t", $line ) )[3];
#	my $value = ( split( "\t", $line ) )[2];
#	if ( $key ne "NA" ) {
#		$gene_gi{$key} =
#		    ( exists $gene_gi{$key} )
#		  ? ( $gene_gi{$key} . ";" . $value )
#		  : $value;
#	}
#}
#close $fh_gi;
#
##hi
#while ( my $line = <$fh_hi> ) {
#	chomp($line);
#	my $key   = ( split( "\t", $line ) )[3];
#	my $value = ( split( "\t", $line ) )[1];
#	if ( $key ne "NA" ) {
#		$gene_hi{$key} =
#		    ( exists $gene_hi{$key} )
#		  ? ( $gene_hi{$key} . ";" . $value )
#		  : $value;
#	}
#}
#close $fh_hi;

#gnomAD constraint, updated Mar06,2019
#1       gene
#2       transcript
#3       canonical
#4       obs_lof
#5       exp_lof
#6       oe_lof
#7       oe_lof_lower
#8       oe_lof_upper
#9       obs_mis
#10      exp_mis
#11      oe_mis
#12      oe_mis_lower
#13      oe_mis_upper
#14      obs_syn
#15      exp_syn
#16      oe_syn
#17      oe_syn_lower
#18      oe_syn_upper
#19      lof_z
#20      mis_z
#21      syn_z
#22      pLI
#23      pRec
#24      pNull
#25      gene_issues
#26      range_lof_z
#27      range_mis_z
#28      range_syn_z
#29      range_pLI
#30      egID_byOfficialSymbol
#31      egID_bySynonym



my $xmat_header = <$fh_xmt>;
chomp($xmat_header);
my %xmt_hds = ();
my $xmat_idx = 0;
my @xma__headera =   split("\t", $xmat_header);

#gnomAD_oe_lof
#gnomAD_oe_lof_upper
#gnomAD_oe_mis
#gnomAD_oe_mis_upper
#gnomAD_pLI
#gnomAD_pRec

my @xma__required =   ("oe_lof" , "oe_lof_upper" , "oe_mis" , "mis_z",  "oe_mis_upper", "pLI" , "pRec"  );

foreach (@xma__headera){
	$xmt_hds{$_} = $xmat_idx++; 
}

while ( my $line = <$fh_xmt> ) {
	chomp($line);
	my @fields = split("\t" , $line);
	
#	my $key = (exists $xmt_hds{'egID_byOfficialSymbol'}) ? $fields[$xmt_hds{'egID_byOfficialSymbol'}]: die "Error parsing egID_byOfficialSymbol from gnomAD matrix line: ($line)\n";
#	if($key eq "NA"){ next; }
#	foreach my $xm_hd (@xma__headera){
#		$gene_exmat{$key}{$xm_hd} =  (exists $xmt_hds{$xm_hd}) ? $fields[$xmt_hds{$xm_hd}]: die "Error parsing $xm_hd from gnomAD matrix line: ($line)\n";
#	}

	my $key = (exists $xmt_hds{'consensus_egID'}) ? $fields[$xmt_hds{'consensus_egID'}]: die "Error parsing \'consensus_egID\' from gnomAD matrix line: ($line)\n";
	if($key eq "NA"){ next; }
	foreach my $xm_hd (@xma__headera){
		$gene_exmat{$key}{$xm_hd} =  (exists $xmt_hds{$xm_hd}) ? $fields[$xmt_hds{$xm_hd}]: die "Error parsing $xm_hd from gnomAD matrix line: ($line)\n";
	}
	
#	my $key = (exists $xmt_hds{'egID_byOfficialSymbol'}) ? $fields[$xmt_hds{'egID_byOfficialSymbol'}]: die "Error parsing egID_byOfficialSymbol from gnomAD matrix line: ($line)\n";
#	$gene_exmat{$key}{'mis_z'} = (exists $xmt_hds{'mis_z'}) ? $fields[$xmt_hds{'mis_z'}]: die "Error parsing mis_z from gnomAD matrix line: ($line)\n";
#	$gene_exmat{$key}{'lof_z'} = (exists $xmt_hds{'lof_z'}) ? $fields[$xmt_hds{'lof_z'}]: die "Error parsing lof_z from gnomAD matrix line: ($line)\n";
#	$gene_exmat{$key}{'pli'} = (exists $xmt_hds{'pLI'}) ? $fields[$xmt_hds{'pLI'}]: die "Error parsing pLI from gnomAD matrix line: ($line)\n";
}
close $fh_xmt;
warn "Done loading Phenotype data\n";

#
#while ( my $line = <$fh_xmt> ) {
#	chomp($line);
##	my ($key,$tr, $syn_z, $mis_z, $lof_z, $pli) = split("\t" , $line);
#	my ($tr, $syn_z, $mis_z, $lof_z, $pli, $key, $g_name) = split("\t" , $line);
#	#warn "$key,$tr, $syn_z, $mis_z, $lof_z, $pli\n";
#	$gene_exmat{$key}{'mis_z'} = $mis_z;
#	$gene_exmat{$key}{'lof_z'} = $lof_z;
#	$gene_exmat{$key}{'pli'} = $pli;
#}
#close $fh_xmt;
#warn "Done loading Phenotype data\n";
#


#dbsnfp header , added Jan 27,2015
#open(DBNSF, $dbnsfp_db) or usage("Cannot open $dbnsfp_db $!\n");
my $fh_dbnsf = IO::File->new($dbnsfp_db , q{<}) or usage("Cannot open dbnsfp_db  $dbnsfp_db $!\n");
my $dbnsfp_db_header = <$fh_dbnsf>;
chomp($dbnsfp_db_header);
close $fh_dbnsf;
if($dbnsfp_db_header !~ /^#/){
	usage("Critical: Cannot extract the header from $dbnsfp_db_header ...exiting.. \n")
}

#dbscsnv header , added Feb 09,2015
#open(DBSC, $dbscsnv_db) or usage("Cannot open $dbscsnv_db $!\n");
my $fh_dbsc = IO::File->new($dbscsnv_db , q{<}) or usage("Cannot open dbscsnv_db  $dbscsnv_db $!\n");
my $dbscsnv_db_header = <$fh_dbsc>;
chomp($dbscsnv_db_header);
close $fh_dbsc;
if($dbscsnv_db_header !~ /^#/){
	usage("Critical: Cannot extract the header from $dbscsnv_db_header ...exiting.. \n")
}



#open(GADG, $gnADgn_db) or usage("Cannot open $gnADgn_db $!\n");
#my $fh_gadg = IO::File->new($gnADgn_db , q{<}) or usage("Cannot open $gnADgn_db $!\n");
#
#my $gnADgn_db_header = <$fh_gadg>;
#chomp($gnADgn_db_header);
#warn "gnADgn_db_header == $gnADgn_db_header\n";
#my @gnADgn_hds = ();
#my $gnADgn_index = 0;
#foreach (split("\t" , $gnADgn_db_header)){
#	if($gnADgn_index++ > 4){
#		push(@gnADgn_hds , $_);
#	}
#}
#close $fh_gadg;
#if($gnADgn_db_header !~ /^#/){
#	usage("Critical: Cannot extract the header from $gnADgn_db_header ...exiting.. \n")
#}
#my @gnADgn_db_hds = map { "gnomAD_genome30_"  . $_ } split("\t" , $gnADgn_db_header );

if ( $orgi_input eq 'stdin' ) {
	*MVAR = *STDIN;
	warn "** input is from stdin ..\n";
}
elsif ( $orgi_input =~ m/.+\.gz$/ ) {
	warn "** input is gz ..\n";
	open( MVAR, "gunzip -c $orgi_input |" )
	  or usage("Error reading orginal input , $orgi_input ( $! )");
}
elsif ( $orgi_input =~ m/.+\.bz2$/ ) {
	open( MVAR, "bzip2 -c -d  $orgi_input |" )
	  or usage("Error reading orginal input , $orgi_input ( $! )");
}
else {
	open( MVAR, $orgi_input )
	  or usage("Error reading orginal input , $orgi_input ( $! )");
}


my @vcf_header = ();
$mvar_header = "NA";
my $ovcf_key_index = -1;
while ( my $line = <MVAR> ) {
	chomp($line);
	if ( $line =~ m/^#.+/ || $line =~ m/^>.+/ || $line eq "" ) {
		if ( $line =~ m/^>.+/ ) {
			$mvar_header = $line;
			
			$mvar_header =~ s/>//;
			print $fh_invo $mvar_header, "\n";
		}
		push( @vcf_header, $line );
	}else{
		last;
	}
}
close MVAR;


foreach (@vcf_header){
	if($_ =~ m/^#CHROM/){
		my $chrom_line = $_;
		my $idx = 0;
		foreach my $cvf (split("\t" , $chrom_line)){
			if($cvf eq "OVCF_KEY"){
				$ovcf_key_index = $idx;
				
			}
			$idx++;
		} 
	}
}

#Non-parallel loading , Feb 05, 2016
my %retrieved_hashes = ();
my $total_2load = scalar (keys %files_parallel);

my $load_count = 0;
foreach my $m (keys %files_parallel){
	$load_count++;
	my $idt = $pr_names{$m};
   	my $start_time = localtime;
	my $sub_type = $files_parallel{$m};
	next if($sub_type eq "e") ;
	warn "\[$start_time\] ** Loading $idt ($load_count\/$total_2load)\n";
	if($sub_type eq "r" || $sub_type eq "f"){
		$retrieved_hashes{$idt} =  load_hash($m);
	}elsif($sub_type eq "g"){
		warn "Skipped\n";
		#$retrieved_hashes{$idt} = load_genev_hash($m);
	}elsif($sub_type eq "x" || $sub_type eq "ux"){
		#$retrieved_hashes{$idt} = load_genex_hash_ext($m, $sub_type);
		warn "Skipped\n";
	}elsif($sub_type eq "i"){
		$retrieved_hashes{$idt} = load_input_hash($m, $run_type);
	}elsif($sub_type eq "w"){
		$retrieved_hashes{$idt} = load_hashW($m);
	}elsif($sub_type eq "b"){
		$retrieved_hashes{$idt} = load_hashBedDistance($m);
	}elsif($sub_type eq "gex"){
		$retrieved_hashes{$idt} = load_geneExtra($m);
	}elsif($sub_type eq "v"){
		$retrieved_hashes{$idt} = load_vcfanno($m);
	}elsif($sub_type eq "vt"){
		$retrieved_hashes{$idt} = load_vcfannoText($m);
	}elsif($sub_type eq "ifr"){
		$retrieved_hashes{$idt} = load_cgInt($m);
	}
   	my $end_time = localtime;
   	warn  "\[$end_time\] Done $idt \n";  	
}


my $start_time = localtime;
warn "\[$start_time\] *** Generating Table.... \n";
my $text_has_header = 1;

my @spx_hds = ();
foreach (split(":" , $spidex_header)){
	push(@spx_hds , "spx_" . $_);
}
my @spx_hds2go = ();

my @spliceAIsnv_hds = ();
foreach (split(":" , $spliceAIsnv_header)){
	push(@spliceAIsnv_hds , "spliceAI_" . $_);
}

#Output configuration , added June 09,2017
my @cols_show = ();
my %cols_config = ();
#open(OCFG, $ouput_config_file) or die "$! $ouput_config_file\n";
my $fh_ocfg = IO::File->new($ouput_config_file , q{<}) or usage("Cannot open $ouput_config_file $!\n");

while(my $ocfg_line = <$fh_ocfg>){
	chomp($ocfg_line);
	my $col_name = (split("\t" , $ocfg_line))[0];
	$cols_config{$col_name} = 1;
	if($ocfg_line =~ /^#/){}else{
		push(@cols_show , (split("\t" , $ocfg_line))[0] );
		$cols_config{$col_name} = 0;
	}
}
close $fh_ocfg;
my @exac_cols = ("ExAC_Freq","ExAC_AFR","ExAC_AMR","ExAC_EAS","ExAC_FIN","ExAC_NFE","ExAC_OTH","ExAC_SAS");

my $fh_outvcf = undef;
if($vcfout == 1){
	my $out_vcfname = $orgvcf_in_basename . "_final.vcf";
	open($fh_outvcf , ">$out_vcfname");
}


#my $fh_org = undef;
#if ($orgvcf_in =~ /\.gz$/) {
#	#warn "gzipped input file\n";
#	open($fh_org, "gunzip -c $orgvcf_in |") || die "can’t open pipe to $orgvcf_in";
#}elsif($orgvcf_in =~ /\.vcf$/){
#	#warn "$orgvcf_in normal input file\n";
#	$fh_org = IO::File->new($orgvcf_in , q{<}) or usage("$! cannnot load input vcf $orgvcf_in\n");
#}else{
#	die "Cannot determine input file type (vcf or gz ?)..";
#}	

my @headers = ();
#while(my $line = <$fh_org>){
#	chomp($line);
#	if($line =~ m/^#/){
#		push(@headers , $line);
#	}else{
#		last;
#	}
#}

#my $ret_header = add2header(\@headers , "Gene", "A" , "Something");
#$ret_header = add2header($ret_header , "typeseq_pr", "A" , "typesee");
#$ret_header = add2header($ret_header , "GENE_Source", "A" , "ucsc,refseq");
#$ret_header = add2header($ret_header , "effect", "A" , "eff");
#$ret_header = add2header($ret_header , "Clinvar_alleleID", "A" , "clinvar");
#$ret_header = add2header($ret_header , "Clinvar_SIG", "A" , "clinvar");
#$ret_header = add2header($ret_header , "dbsnp", "A" , "dbsnp v150");
#$ret_header = add2header($ret_header , "SpliceJuncDist", "A" , "splice dist");
#
##dbsnpCommon
#$ret_header = add2header($ret_header , "dbsnpCommon", "A" , "dbsnpCommon v150");
#$ret_header = add2header($ret_header , "dbsnpRegion", "A" , "dbsnpRegion v150");
#$ret_header = add2header($ret_header , "dbsnpWindow", "A" , "dbsnpWindow padded 7bp,  v150");
#$ret_header = add2header($ret_header , "SIFT_score", "A" , "sift score from dbsnfp");
#$ret_header = add2header($ret_header , "polyphen_score", "A" , "polyphen_score from dbsnfp");
#$ret_header = add2header($ret_header , "PROVEAN_score", "A" , "PROVEAN_score from dbsnfp");
#$ret_header = add2header($ret_header , "MA_score", "A" , "MutationAssessor_score from dbsnfp");
#$ret_header = add2header($ret_header , "MT_score", "A" , "MutationTaster_score from dbsnfp");
##108     CADD_Raw
##109     CADD_phred
#$ret_header = add2header($ret_header , "CADD_Raw", "A" , "CADD_Raw");
#$ret_header = add2header($ret_header , "CADD_Phred", "A" , "CADD_phred");
#$ret_header = add2header($ret_header , "Cosmic", "A" , "Cosmic info");
##Phylop_100
#$ret_header = add2header($ret_header , "Phylop_100", "A" , "Phylop_100");
#$ret_header = add2header($ret_header , "Phylop_20", "A" , "Phylop_20");
#$ret_header = add2header($ret_header , "pfam_ucsc", "A" , "pfam");
#$ret_header = add2header($ret_header , "SegDup", "A" , "segmental dup");
##repeat
#$ret_header = add2header($ret_header , "Repeat", "A" , "rmsk");
##dbscsnv
#
#$ret_header = add2header($ret_header , "dbscSNV_ADA_SCORE", "A" , "dbscSNV_ADA_SCORE");
#$ret_header = add2header($ret_header , "dbscSNV_RF_SCORE", "A" , "dbscSNV_RF_SCORE");

my $ret_header  = add2headerFull(\@headers , \@cols_show);

if($vcfout == 1){
	print $fh_outvcf join("\n" , @{$ret_header}) , "\n";
}

#warn "*** RUN TYPE = $run_type\n";

my $text_header = "-1";
if ( $run_pp == 5 ) {
	print "#" , join("\t" , @cols_show) ;
	$text_header = <$fh_sam>;
}elsif ( $run_type eq "v" ) {
	my $v_header = join( "\n", @vcf_header );
	print $v_header , "\t" , join("\t" , @cols_show) , "\n"; 
}elsif ( $run_type eq "t" ) {
	$text_header = <$fh_sam>;
	if ( $text_header =~ m/^#/ ) {
		chomp($text_header);
#		$text_header =~ s/^\s+//;    #remove leading spaces
#		$text_header =~ s/\s+$//;    #remove trailing spaces
		$text_header =~ s/^\s*(\S*(?:\s+\S+)*)\s*$/$1/;
		#@array = @array[ 5 .. $#array ];
		my @text_headers = split("\t" , $text_header);
		my @text_headers2go = @text_headers[5 .. $#text_headers];
		  #print "chromosome\tbegin\tend\tannovar_chr\tannovar_start\tannovar_end\tref_allele\talt_allele\tsequenceOverlap\trefseq-id\taa_flag\teffect\tleftD\trightD\tgene-symbol\tentrez-id\tgene_desc\tgene_type\tomim_id\tomim_Phenotype\tMPO\tHPO\tCGD_disease\tCGD_inheritance\tGene_GI\tGene_HI\tdbsnp\tdbsnp_common\tdbsnp_region\tdbsnp_common_region\tdbsnp_wind\tClinvar_SIG\tClinvar_CLNDBN\tClinvar_CLNACC\t1000g_all\t1000g_eur\t1000g_amr\t1000g_eas\t1000g_afr\t1000g_sas\tNHLBI_all\tNHLBI_aa\tNHLBI_eu\tExAC_Freq\tExAC_AFR\tExAC_AMR\tExAC_EAS\tExAC_FIN\tExAC_NFE\tExAC_OTH\tExAC_SAS\tcg\tcg_filtered\tcgW597_AllFreq\tcgW597_CalledFreq\tcgW597_11s\tcgW597_Hs\tcgW597_Ls\tcg1KG436_AllFreq\tcg1KG436_CalledFreq\tcg1KG436_11s\tcg1KG436_Hs\tcg1KG436_Ls\tsift_score\tpolyphen_score\tma_score\tmt_score\tcosmic\tphylopPMam\tphylopPMam_avg\tphylopVert100\tphylopVert100_avg\tphastCons_placental\tpfam_annovar\tper_cds_affected\tper_transcripts_affected\tSegDup\tHGMD_Accession\tHGMD_type\tHGMD_tag\tHGMD_Disease\tHGMD_PubmedId\tRepeat\tCADD_Raw\tCADD_phred\n";
		print join("\t" , @text_headers2go), "\t",
#"annovar_chr\tannovar_start\tannovar_end\tref_allele\talt_allele\ttypeseq\ttypeseq_priority\trefseq_id\taa_flag\teffect\teffect_priority\tleftD\trightD\tgene_symbol\tentrez_id\tgene_desc\tgene_type\tomim_id\tomim_phenotype\tMPO\tHPO\tCGD_disease\tCGD_inheritance\tExAc_mis_z\tExAc_lof_z\tExAc_pLI\tClinvar_mis_pathogenic_count\tClinvar_mis_benign_count\tClinvar_lof_pathogenic_count\tClinvar_lof_benign_count\tClinvar_SIG\tClinvar_CLNREF\tClinvar_CLNACC\tHGMD_Accession\tHGMD_type\tHGMD_tag\tHGMD_Disease\tHGMD_PubmedId\tcosmic\t1000g_all\t1000g_eur\t1000g_amr\t1000g_eas\t1000g_afr\t1000g_sas\tNHLBI_all\tNHLBI_aa\tNHLBI_eu\tExAC_Freq\tExAC_AFR\tExAC_AMR\tExAC_EAS\tExAC_FIN\tExAC_NFE\tExAC_OTH\tExAC_SAS\tcg\tcg_filtered\tcgW597_AllFreq\tcgW597_CalledFreq\tcgW597_11s\tcgW597_Hs\tcgW597_Ls\tcg1KG436_AllFreq\tcg1KG436_CalledFreq\tcg1KG436_11s\tcg1KG436_Hs\tcg1KG436_Ls\t" , join("\t" , @gnADex_hds) , "\t" , join("\t" , @gnADgn_hds) , "\tdbsnp\tdbsnp_common\tdbsnp_region\tdbsnp_common_region\tdbsnp_wind\tsift_score\tPROVEAN_score\tpolyphen_score\tma_score\tmt_score\tCADD_Raw\tCADD_phred\tphylopPMam\tphylopPMam_avg\tphylopVert100\tphylopVert100_avg\tphastCons_placental\tpfam_annovar\tper_cds_affected\tper_transcripts_affected\t" , join("\t" , @spx_hds) , "\t" ,"dbscSNV_ADA_SCORE\tdbscSNV_RF_SCORE\tSegDup\tRepeat\n";
join("\t" , @cols_show) ;
	}
	else {
		$text_header     = "";
		$text_has_header = 0;
		seek $fh_sam , 0, 0;
		print
#"annovar_chr\tannovar_start\tannovar_end\tref_allele\talt_allele\ttypeseq\ttypeseq_priority\trefseq_id\taa_flag\teffect\teffect_priority\tleftD\trightD\tgene_symbol\tentrez_id\tgene_desc\tgene_type\tomim_id\tomim_phenotype\tMPO\tHPO\tCGD_disease\tCGD_inheritance\tExAc_mis_z\tExAc_lof_z\tExAc_pLI\tClinvar_mis_pathogenic_count\tClinvar_mis_benign_count\tClinvar_lof_pathogenic_count\tClinvar_lof_benign_count\tClinvar_SIG\tClinvar_CLNREF\tClinvar_CLNACC\tHGMD_Accession\tHGMD_type\tHGMD_tag\tHGMD_Disease\tHGMD_PubmedId\tcosmic\t1000g_all\t1000g_eur\t1000g_amr\t1000g_eas\t1000g_afr\t1000g_sas\tNHLBI_all\tNHLBI_aa\tNHLBI_eu\tExAC_Freq\tExAC_AFR\tExAC_AMR\tExAC_EAS\tExAC_FIN\tExAC_NFE\tExAC_OTH\tExAC_SAS\tcg\tcg_filtered\tcgW597_AllFreq\tcgW597_CalledFreq\tcgW597_11s\tcgW597_Hs\tcgW597_Ls\tcg1KG436_AllFreq\tcg1KG436_CalledFreq\tcg1KG436_11s\tcg1KG436_Hs\tcg1KG436_Ls\t" , join("\t" , @gnADex_hds) , "\t" , join("\t" , @gnADgn_hds), "\tdbsnp\tdbsnp_common\tdbsnp_region\tdbsnp_common_region\tdbsnp_wind\tsift_score\tPROVEAN_score\tpolyphen_score\tma_score\tmt_score\tCADD_Raw\tCADD_phred\tphylopPMam\tphylopPMam_avg\tphylopVert100\tphylopVert100_avg\tphastCons_placental\tpfam_annovar\tper_cds_affected\tper_transcripts_affected\t" , join("\t" , @spx_hds) , "\t" ,"dbscSNV_ADA_SCORE\tdbscSNV_RF_SCORE\tSegDup\tRepeat\n";
join("\t" , @cols_show) ;
	}
}
else {
	print
#"chromosome\tbegin\tend\tannovar_chr\tannovar_start\tannovar_end\tref_allele\talt_allele\ttypeseq\ttypeseq_priority\trefseq-id\taa_flag\teffect\teffect_priority\tleftD\trightD\tgene-symbol\tentrez-id\tgene_desc\tgene_type\tomim_id\tomim_Phenotype\tMPO\tHPO\tCGD_disease\tCGD_inheritance\tExAc_mis_z\tExAc_lof_z\tExAc_pLI\tClinvar_missense_pathogenic\tClinvar_missense_benign\tClinvar_lof_pathogenic\tClinvar_lof_benign\tdbsnp\tdbsnp_common\tdbsnp_region\tdbsnp_common_region\tdbsnp_wind\tClinvar_SIG\tClinvar_CLNDBN\tClinvar_CLNACC\t1000g_all\t1000g_eur\t1000g_amr\t1000g_eas\t1000g_afr\t1000g_sas\tNHLBI_all\tNHLBI_aa\tNHLBI_eu\tExAC_Freq\tExAC_AFR\tExAC_AMR\tExAC_EAS\tExAC_FIN\tExAC_NFE\tExAC_OTH\tExAC_SAS\tcg\tcg_filtered\tcgW597_AllFreq\tcgW597_CalledFreq\tcgW597_11s\tcgW597_Hs\tcgW597_Ls\tcg1KG436_AllFreq\tcg1KG436_CalledFreq\tcg1KG436_11s\tcg1KG436_Hs\tcg1KG436_Ls\tsift_score\tPROVEAN_score\tpolyphen_score\tma_score\tmt_score\tcosmic\tphylopPMam\tphylopPMam_avg\tphylopVert100\tphylopVert100_avg\tphastCons_placental\tpfam_annovar\tper_cds_affected\tper_transcripts_affected\tSegDup\tHGMD_Accession\tHGMD_type\tHGMD_tag\tHGMD_Disease\tHGMD_PubmedId\tRepeat\tCADD_Raw\tCADD_phred\t" , join("\t" , @spx_hds) , "\tdbscSNV_ADA_SCORE\tdbscSNV_RF_SCORE\n";
"chromosome\tbegin\tend\t" , join("\t" , @cols_show) ; 
}


if($calc_impact == 1) {
	print "\tImpact\tComment\n";
}else{
	print "\tComment\n";
}

my $tot_var                 = 0;
my $tot_var_with_genes                 = 0;
my $multi_gene              = 0;
my $multi_typeseq			= 0;
my $no_entrez_id_genes      = 0;
my %no_entrez_id_genes_list = ();
my %type_warnings           = ();
$stat_count = 0;
my $multi_typeseq_per_gene = 0;
my $iter = 0;

my %fh_sam_header = ();
my $fh_sam_hd_idx = 0;
foreach (split("\t" , $text_header)){
	$fh_sam_header{$_} = $fh_sam_hd_idx++;
}


while ( my $line = <$fh_sam> ) {
 	#added Dec 2014 to correctly handle multigene
	my %multigene_tr_affected       = ();
	my %multigene_total_transcripts = ();
	my %final_out = ();
	my %impact_col = ();
	my @freq_max = ();
	my @x1000g_freq_max = ();
	my @ExAC_freq_max = ();
	my @cg_freq_max = ();
	my @gnomAD_exome_freq_max = ();
	my @gnomAD_genome_freq_max = (0);	

	$tot_var++;
	chomp($line);
	
	#chr1     13813   13813   T       G       chr1    13813   .       T       G       121.8   VQSRTrancheSNP99.95to100.00     AC=2;AF=1;AN=2;DP=4;ExcessHet=3.0103;FS=0;MLEAC=1;MLEAF=0.5;MQ=23;QD=25.36;SOR=2.833;VQSLOD=-36.01;culprit=MQ;OVCF_KEY=chr1:13813:T:G;MULTI_ALLELIC=0   OGT:OZYG:GT:AD:DP:GQ:PGT:PID:PL 1/1:hom-alt:1/1:0,3:3:9:1|1:13813_T_G:135,9,0
	my @fields_an      = split( "\t", $line );
	my $fields_size_an = scalar @fields_an;
	my @vcf_fields     = ();
	my @text_fields    = ();
	my @text_fields2go    = ();
	my $locus_id       = "";
	my $ovcf_value = undef;
	
	my @fields = ();
	my $fields_size = 0;
	
	if($vcfout == 1){
		@vcf_fields = @fields_an[ 5 .. ( @fields_an - 1 ) ];
		$ovcf_value = get_ovcfkey($vcf_fields[7]);
	}
	if ( $run_type eq "m" ) {
		$locus_id = extract_locusid($line);
		$locus_id =~ s/^\s*(\S*(?:\s+\S+)*)\s*$/$1/;
		my $mvar_val       = $mvar{$locus_id};
		my $mvar_val_count = $mvar_val =~ tr/\t/\t/;
	
		@fields      = split( "\t", $mvar_val, $mvar_val_count );
		$fields_size = scalar @fields;
		
	}
	elsif ( $run_type eq "v" || $run_type eq "v2" ) {
		# @vcf_fields = @fields_an[ 8 .. ( @fields_an - 1 ) ];
		@vcf_fields = @fields_an[ 5 .. ( @fields_an - 1 ) ];   # this has to adjust for vcf4old (5 .. ), and for vcf4 (8 ..) , October 2, 2017
		
	}elsif ( $run_type eq "t" && $text_has_header == 1 ) {
		@text_fields2go = @fields_an[5 .. $#fields_an]; 
		#@vcf_fields = @fields_an[ 5 .. ( @fields_an - 1 ) ];
		#warn "TEXT_line = " ,join("\t" , @text_fields2go) , "\n";
	}
	
#	$locus_id =~ s/^\s+//;    #remove leading spaces
#	$locus_id =~ s/\s+$//;    #remove trailing spaces
	
	my $key =
	    $fields_an[0] . ":"
	  . $fields_an[1] . ":"
	  . $fields_an[2] . ":"
	  . $fields_an[3] . ":"
	  . $fields_an[4];

	my $region_key = $fields_an[0] . ":" . $fields_an[1] . ":" . $fields_an[2];
	my $is_hg38_ambg = 0;
	my @anno_comments = ();
	if(exists $hg38_ambg{$region_key}){
		$is_hg38_ambg = 1;
		push(@anno_comments, $hg38_ambg{$region_key});
	}
	
	
	if($run_pp == 5){  
		$final_out{'Annovar_key'} = $key;
		#extract mssngkey
		$final_out{'MSSNG_key'} = extract_mssngkey($fields_an[12]);
		$final_out{'id'} = 		 $final_out{'MSSNG_key'};

		#reference_name
		#start
		#end
		#reference_bases
		#alternate_bases
		#chr1-952179-952180-A-C
		my @mssng_coords = split('-' , $final_out{'MSSNG_key'});
		if(scalar(@mssng_coords) < 4){
			usage("Critical: MSSNG_key " . $final_out{'MSSNG_key'} . " has wrong format!\n");
		}		
		$final_out{'reference_name'} = $mssng_coords[0];
		$final_out{'start'} =   (exists $fh_sam_header{'POS'}) ?  ($fields_an[$fh_sam_header{'POS'}] -1) : usage("Critical! Cannot extract POS from header!");   #$mssng_coords[1];
		$final_out{'end'} =   (exists $fh_sam_header{'ann_end'}) ?  $fields_an[$fh_sam_header{'ann_end'}] : usage("Critical! Cannot extract ann_end from header!");  # $mssng_coords[2];
		$final_out{'reference_bases'} =  (exists $fh_sam_header{'REF'}) ?  $fields_an[$fh_sam_header{'REF'}] : usage("Critical! Cannot extract REF from header!"); # $mssng_coords[3];
		$final_out{'alternate_bases'} =  (exists $fh_sam_header{'ALT'}) ?  $fields_an[$fh_sam_header{'ALT'}] : usage("Critical! Cannot extract ALT from header!");  # $mssng_coords[4];
	}  
	
	#key for region annotation ( phylop , dbsnp region etc.)
	#my $key_dbr2 = $fields_an[0] . ":" . $fields_an[1] . ":" . $fields_an[2];
	#$key_dbr2 =~ s/chr//i;
	#my $key_ncw = $fields_an[0] . ":" . ($fields_an[1] - 7) . ":" . ($fields_an[2] + 7);

	# other annotations
	my $cg_val = 0;
	my $cgf_val = 0;
	
	if ( exists $retrieved_hashes{'CG'}->{$key} ) {
		my $cgboth = $retrieved_hashes{'CG'}->{$key};
		$cg_val = (split(":", $cgboth))[0];
		$cg_val =~ s/U=//;
		$cgf_val = (split(":", $cgboth))[1];
		$cgf_val =~ s/F=//;
	}
	$final_out{'cg'} = $cg_val;
	$final_out{'cg_filtered'} = $cgf_val;
#    push(@freq_max , $cg_val);
#    push(@freq_max , $cgf_val);
    push(@cg_freq_max , $cg_val);
    push(@cg_freq_max , $cgf_val);
    
    
	my $wld_val = "NA";
	my $w597_AllFreq = 0; my $w597_CalledFreq = 0 ; my $w597_11s = 0 ; my $w597_Hs = 0 ; my $w597_Ls = 0;
	if ( exists  $retrieved_hashes{'Welderly'}->{$key} ) {
		$wld_val = $retrieved_hashes{'Welderly'}->{$key};
		if($wld_val =~ m/$pat_weld/){
			$w597_AllFreq = $1;
			$w597_CalledFreq = $2;
			$w597_11s = $3;
			$w597_Hs  = $4;
			$w597_Ls = $5;
		} 
	}
	$final_out{'cgW597_AllFreq'} = $w597_AllFreq;
	$final_out{'cgW597_CalledFreq'} = $w597_CalledFreq;
	$final_out{'cgW597_11s'} = $w597_11s; 
	$final_out{'cgW597_Hs'} = $w597_Hs;
	$final_out{'cgW597_Ls'} = $w597_Ls;
   # push(@freq_max , $w597_AllFreq);
   # push(@freq_max , $w597_CalledFreq);
    push(@cg_freq_max , $w597_AllFreq);
    push(@cg_freq_max , $w597_CalledFreq);
	

	my $ckb_val = "NA";
	#my $cg1kb_AllFreq = 0; my $cg1kb_CalledFreq = 0; my $cg1kb_11s = 0 ; my $cg1kb_Hs = 0 ; my $cg1kb_Ls = 0;
	if ( exists $retrieved_hashes{'cg1KB'}->{$key} ) {
		$ckb_val = $retrieved_hashes{'cg1KB'}->{$key};
		if($ckb_val =~ m/$pat_cg2/){
#			$cg1kb_AllFreq = $1;
#			$cg1kb_CalledFreq = $2;
#			$cg1kb_11s = $3;
#			$cg1kb_Hs  = $4;
#			$cg1kb_Ls = $5;
			$final_out{'cg1KG436_AllFreq'} = $1;
			$final_out{'cg1KG436_CalledFreq'} = $2;
			$final_out{'cg1KG436_11s'} = $3;
			$final_out{'cg1KG436_Hs'} = $4;
			$final_out{'cg1KG436_Ls'} = $5;
		}		
	}else{
			$final_out{'cg1KG436_AllFreq'} = 0;
			$final_out{'cg1KG436_CalledFreq'} = 0;
			$final_out{'cg1KG436_11s'} = 0;
			$final_out{'cg1KG436_Hs'} = 0;
			$final_out{'cg1KG436_Ls'} = 0;
	}
	
    #push(@freq_max , $final_out{'cg1KG436_AllFreq'});
    #push(@freq_max , $final_out{'cg1KG436_CalledFreq'});
    push(@cg_freq_max , $final_out{'cg1KG436_AllFreq'});
    push(@cg_freq_max , $final_out{'cg1KG436_CalledFreq'});
	
	$final_out{'NHLBI_all'} = ( exists  $retrieved_hashes{'NHLBI'}->{$key} ) ? $retrieved_hashes{'NHLBI'}->{$key} : 0 ;
	
	$final_out{'dbsnp'} = ( exists  $retrieved_hashes{'dbsnp'}->{$key} ) ? $retrieved_hashes{'dbsnp'}->{$key} : "NA" ;
	$impact_col{'dbsnp'} = $final_out{'dbsnp'};
	$stats_final_merge{counts}{'dbsnp'} = ($final_out{'dbsnp'} ne "NA") ? ($stats_final_merge{counts}{'dbsnp'} + 1) : $stats_final_merge{counts}{'dbsnp'};

	$final_out{'dbsnp_common'} = ( exists  $retrieved_hashes{'dbsnpCommon'}->{$key} ) ? $retrieved_hashes{'dbsnpCommon'}->{$key} : "NA" ;
	$impact_col{'dbsnpc'} = $final_out{'dbsnp_common'};
	$stats_final_merge{counts}{'dbsnp_common'} = ($final_out{'dbsnp_common'} ne "NA") ? ($stats_final_merge{counts}{'dbsnp_common'} + 1) : $stats_final_merge{counts}{'dbsnp_common'};
	
	
	my $cln_val = 'NA,NA,NA,NA,NA,NA';
	if ( exists $retrieved_hashes{'clinvar'}->{$key} ) {
		$cln_val = $retrieved_hashes{'clinvar'}->{$key};
		#warn "Clin = $key  $cln_val\n";
	}

	#Updated July 29 with  in-house clinvar
	#Clinvar_SIG
	#Clinvar_CLNREF
	#Clinvar_CLNACC
	
	#Since Aug 2017
	#RCV000410214@Likely pathogenic@Likely pathogenic@criteria provided, single submitter@20301788,20631546
	
	#my ($Clinvar_3 , $Clinvar_1 , $Clinvar_1ord , $Clinvar_2, $Clinvar_reviewStat) = split('@', $cln_val);
	#my ($Clinvar_acc , $Clinvar_sig , $Clinvar_sig_ord ,  $Clinvar_reviewStat, $Clinvar_ref) = split('@', $cln_val);
	#514941,Myasthenic_syndrome\x2c_congenital\x2c_8,MedGen:C3808739\x2cOMIM:615120,criteria_provided\x2c_single_submitter,Uncertain_significance 
	
	#updated on Jan21,2019
	 #AF_ESP=0.68430;AF_EXAC=0.59339;AF_TGP=0.56130;
	 #ALLELEID=249669;CLNDISDB=MedGen:CN169374|MedGen:CN517202;
	 #CLNDN=not_specified|not_provided;
	 #CLNHGVS=NC_000001.11:g.2025598T>C;
	 #CLNREVSTAT=criteria_provided,_multiple_submitters,_no_conflicts;
	 #CLNSIG=Benign;
	 #CLNVC=single_nucleotide_variant;CLNVCSO=SO:0001483;GENEINFO=GABRD:2563;MC=SO:0001819|synonymous_variant;ORIGIN=1;RS=2229110
	 
	#my ($Clinvar_acc , $Clinvar_pheno , $Clinvar_ref, $Clinvar_reviewStat , $Clinvar_sig) = split(',', $cln_val);
	 #Benign,MedGen:CN517202,not_provided,657839,criteria_provided\x2c_single_submitter,0
	my ($Clinvar_sig , $Clinvar_ref , $Clinvar_pheno, $Clinvar_acc ,  $Clinvar_reviewStat , $clinvar_sigsim) = split(',' , $cln_val);
	#my ($Clinvar_acc , $Clinvar_pheno , $Clinvar_ref, $Clinvar_reviewStat , $Clinvar_sig, $clinvar_sigsim) = parse_clinvar($cln_val);
	$final_out{'Clinvar_SIG'} = ($Clinvar_sig ne "") ?  $Clinvar_sig : "NA";
	$final_out{'Clinvar_CLNREF'} = ($Clinvar_ref ne "") ?   $Clinvar_ref : "NA";;
	$final_out{'Clinvar_AlleleID'} = ($Clinvar_acc ne "") ? $Clinvar_acc : "NA";
	$final_out{'Clinvar_ReviewStatus'} = ($Clinvar_reviewStat ne "") ?  $Clinvar_reviewStat : "NA";
	$final_out{'Clinvar_Disease'} =  ($Clinvar_pheno ne "") ?  $Clinvar_pheno : "NA";
	$final_out{'Clinvar_SIG_Simple'} =   ($clinvar_sigsim ne "") ? $clinvar_sigsim : "NA";
	
	$stats_final_merge{counts}{'Clinvar_SIG'} = ($final_out{'Clinvar_SIG'} ne "NA") ? ($stats_final_merge{counts}{'Clinvar_SIG'} + 1) : $stats_final_merge{counts}{'Clinvar_SIG'};
	$stats_final_merge{'Clinvar_stats'}{Pathogenic} = ($final_out{'Clinvar_SIG_Simple'} == 1) ? ($stats_final_merge{'Clinvar_stats'}{Pathogenic} + 1) : $stats_final_merge{'Clinvar_stats'}{Pathogenic};
		
	my $dnf_val = "NA";
	if ( exists $dnf{$key} ) {
		$dnf_val = $dnf{$key};
	}
	$final_out{'A1000g_all'} = ( exists  $retrieved_hashes{'1000g'}->{$key} ) ? $retrieved_hashes{'1000g'}->{$key} : 0;	
	$stats_final_merge{counts}{'1000g_all'} = ($final_out{'1000g_all'} != 0) ? ($stats_final_merge{counts}{'1000g_all'} + 1) : $stats_final_merge{counts}{'1000g_all'};
	$final_out{'A1000g_amr'} = ( exists  $retrieved_hashes{'1000g_AMR'}->{$key} ) ? $retrieved_hashes{'1000g_AMR'}->{$key} : 0;	
	$final_out{'A1000g_eas'} = ( exists  $retrieved_hashes{'1000g_EAS'}->{$key} ) ? $retrieved_hashes{'1000g_EAS'}->{$key} : 0;	
	$final_out{'A1000g_sas'} = ( exists  $retrieved_hashes{'1000g_SAS'}->{$key} ) ? $retrieved_hashes{'1000g_SAS'}->{$key} : 0;	
	$final_out{'A1000g_eur'} = ( exists  $retrieved_hashes{'1000g_EUR'}->{$key} ) ? $retrieved_hashes{'1000g_EUR'}->{$key} : 0;	
	$final_out{'A1000g_afr'} = ( exists  $retrieved_hashes{'1000g_AFR'}->{$key} ) ? $retrieved_hashes{'1000g_AFR'}->{$key} : 0;	
	$final_out{'NHLBI_aa'} = ( exists  $retrieved_hashes{'NHLBI aa'}->{$key} ) ? $retrieved_hashes{'NHLBI aa'}->{$key} : 0;	
	$final_out{'NHLBI_ea'} = ( exists  $retrieved_hashes{'NHLBI ea'}->{$key} ) ? $retrieved_hashes{'NHLBI ea'}->{$key} : 0;	
	#$final_out{'intervar'} = ( exists  $retrieved_hashes{'intervar'}->{$key} ) ? $retrieved_hashes{'intervar'}->{$key} : "NA";	
	#$final_out{'Clinvar_annovar'} = ( exists  $retrieved_hashes{'Clinvar_annovar'}->{$key} ) ? $retrieved_hashes{'Clinvar_annovar'}->{$key} : "NA";	
	#$final_out{'gerp_wgs'} = ( exists  $retrieved_hashes{'gerp_wgs'}->{$key} ) ? $retrieved_hashes{'gerp_wgs'}->{$key} : "NA";	
#	$final_out{'gerp_elem'} = ( exists  $retrieved_hashes{'gerp_elem'}->{$key} ) ? $retrieved_hashes{'gerp_elem'}->{$key} : "NA";

#    push(@freq_max , $final_out{'A1000g_all'});
#    push(@freq_max , $final_out{'A1000g_amr'});
#    push(@freq_max , $final_out{'A1000g_eas'});
#    push(@freq_max , $final_out{'A1000g_sas'});
#    push(@freq_max , $final_out{'A1000g_eur'});
#    push(@freq_max , $final_out{'A1000g_afr'});
    push(@x1000g_freq_max , $final_out{'A1000g_all'});
    push(@x1000g_freq_max , $final_out{'A1000g_amr'});
    push(@x1000g_freq_max , $final_out{'A1000g_eas'});
    push(@x1000g_freq_max , $final_out{'A1000g_sas'});
    push(@x1000g_freq_max , $final_out{'A1000g_eur'});
    push(@x1000g_freq_max , $final_out{'A1000g_afr'});

     
#    #gerpwgs parallel
    my $gerp_wgs_parallel = "NA";
	my $this_chr = $fields_an[0];
	$this_chr =~ s/chr//;
	my $out_key = "gerpwgs:" .  $this_chr;
	#warn "####gerpwgs: $out_key\n";
	if($retrieved_hashes{$out_key}->{$key}){
		$gerp_wgs_parallel = $retrieved_hashes{$out_key}->{$key};
	}
	$final_out{'gerp_wgs'} = $gerp_wgs_parallel;
	$stats_final_merge{counts}{'gerp_wgs'} = ($final_out{'gerp_wgs'} ne "NA") ? ($stats_final_merge{counts}{'gerp_wgs'} + 1) : $stats_final_merge{counts}{'gerp_wgs'};

	#added on Nov11, 2014 ; modified June 2017
	my $exac_val = "0,0,0,0,0,0,0,0";
	#my $each_exact_val = "0\t0\t0\t0\t0\t0\t0\t0";
	if ( exists $retrieved_hashes{'exac'}->{$key} ) {
		$exac_val = $retrieved_hashes{'exac'}->{$key};
	}
		
	my $exc_col_count = 0;
	foreach my $each_exact (split("," , $exac_val)){
		my $exac_val_go = ($each_exact eq ".") ? 0: $each_exact;
		#push(@freq_max , $exac_val_go);
		#push(@ExAC_freq_max, $exac_val_go);
		$final_out{$exac_cols[$exc_col_count++]} = $exac_val_go;
	}
	
	my $phs_val = "NA";
	if ( exists $retrieved_hashes{'phastcon'}->{$key} ) {
		my $phs_val0 = $retrieved_hashes{'phastcon'}->{$key};
		#Score=467;Name=lod=120 , parsing only score , added April 04, 2016
		if($phs_val0 =~ m/Score=(\d+).+/){
			$phs_val = $1;
			
		}else{
			die "Critical error! unexpected format for phastcons: $phs_val0 , example for expected format : Score=467;Name=lod=120\n "
		}
	}
    $final_out{'phastCons_placental'} = $phs_val;
    $impact_col{'phastCons'} = $phs_val;
    $stats_final_merge{counts}{'phastCons_placental'} = ($final_out{'phastCons_placental'} ne "NA") ? ($stats_final_merge{counts}{'phastCons_placental'} + 1) : $stats_final_merge{counts}{'phastCons_placental'};
	$final_out{'SegDup'} = ( exists  $retrieved_hashes{'segdup'}->{$key} ) ? $retrieved_hashes{'segdup'}->{$key} : "NA";	
    $stats_final_merge{counts}{'segdup'} = ($final_out{'SegDup'} ne "NA") ? ($stats_final_merge{counts}{'segdup'} + 1) : $stats_final_merge{counts}{'segdup'};

##********** fixed Dec 22, 2015
#	my $phy_val = "NA";
#	my $phy_avg_val = "NA";
#	if ( exists $retrieved_hashes{'phylop'}->{$key} ) {
#		$phy_val = $retrieved_hashes{'phylop'}->{$key};
#		$phy_avg_val  = get_phylop_avarage($phy_val);
#	}
#	$final_out{'phylopPMam'} = $phy_val;
#	$stats_final_merge{counts}{'phylopPMam'} = ($final_out{'phylopPMam'} ne "NA") ? ($stats_final_merge{counts}{'phylopPMam'} + 1) : $stats_final_merge{counts}{'phylopPMam'};
#	$final_out{'phylopPMam_avg'} = $phy_avg_val;
#
#	my $phypr_val = "NA";
#	my $phypr_avg_val = "NA";
#	if ( exists $retrieved_hashes{'phylop100w'}->{$key} ) {
#		$phypr_val = $retrieved_hashes{'phylop100w'}->{$key};
#		$phypr_avg_val  = get_phylop_avarage($phypr_val);
#	}
#	$final_out{'phylopVert100'} = $phypr_val;
#	$stats_final_merge{counts}{'phylopVert100'} = ($final_out{'phylopVert100'} ne "NA") ? ($stats_final_merge{counts}{'phylopVert100'} + 1) : $stats_final_merge{counts}{'phylopVert100'};
#	$final_out{'phylopVert100_avg'} = $phypr_avg_val;
#
##*******************************************

#********** updated Nov 5, 2018
	my $phy_val = "NA";
	my $phy_avg_val = "NA";
	my $this_chr = $fields_an[0];
	$this_chr =~ s/chr//;
	my $out_key = "phylop:" .  $this_chr;
	#warn "####phylop: $out_key\n";
	if($retrieved_hashes{$out_key}->{$key}){
		$phy_val = $retrieved_hashes{$out_key}->{$key};
		$phy_avg_val  = get_phylop_avarage($phy_val);
	}
	$final_out{'phylopMam'} = $phy_val;

	$final_out{'phylopMam_avg'} = $phy_avg_val;
	$impact_col{'phymam'} = $phy_avg_val;
	$stats_final_merge{counts}{'phylopMam'} = ($final_out{'phylopMam'} ne "NA") ? ($stats_final_merge{counts}{'phylopMam'} + 1) : $stats_final_merge{counts}{'phylopMam'};

	my $phypr_val = "NA";
	my $phypr_avg_val = "NA";
	#my $this_chr = $fields_an[0];
	#$this_chr =~ s/chr//;
	my $out_key = "phylop100:" .  $this_chr;
	#warn "####gerpwgs: $out_key\n";
	if($retrieved_hashes{$out_key}->{$key}){
		$phypr_val = $retrieved_hashes{$out_key}->{$key};
		$phypr_avg_val  = get_phylop_avarage($phypr_val);
	}
	$final_out{'phylopVert100'} = $phypr_val;
	$final_out{'phylopVert100_avg'} = $phypr_avg_val;
	$impact_col{'phyvert'} = $phypr_avg_val;
	
	$stats_final_merge{counts}{'phylopVert100'} = ($final_out{'phylopVert100'} ne "NA") ? ($stats_final_merge{counts}{'phylopVert100'} + 1) : $stats_final_merge{counts}{'phylopVert100'};


#*********************************

#cg internal

#$final_out{'MSSNG_key'}
if($run_pp == 5){
#	warn "*** for int\n";
	if(not exists $final_out{'MSSNG_key'}){
		usage("MSSNG_key not found in final_out\n");
	}
	my $out_key = "cgInt:" .  $this_chr;
	if(exists $retrieved_hashes{$out_key}->{$final_out{'MSSNG_key'}}){
		my @this_cgintfrqs = split('\t' , $retrieved_hashes{$out_key}->{$final_out{'MSSNG_key'}});
		my $ix = 0;	
		foreach my $this_cgint (@cgInt_header){
			$final_out{$this_cgint} = $this_cgintfrqs[$ix++];
		}
	}else{
		foreach my $this_cgint (@cgInt_header){
			$final_out{$this_cgint} = 0;
		}
		
	}
	
	my $out_key = "illInt:" .  $this_chr;
	if(exists $retrieved_hashes{$out_key}->{$final_out{'MSSNG_key'}}){
		my @this_illintfrqs = split('\t' , $retrieved_hashes{$out_key}->{$final_out{'MSSNG_key'}});
		my $ix = 0;	
		foreach my $this_illint (@illInt_header){
			$final_out{$this_illint} = $this_illintfrqs[$ix++];
		}
	}else{
		foreach my $this_illint (@illInt_header){
			$final_out{$this_illint} = 0;
		}
		
	}	
	
	 
}

# added Jan 27, 2016    extra PROVEAN_score   PROVEAN_pred
    my $pol_val = "NA";
	my $dbnsfp_val = "NA";
	my $ma_val = "NA";
	my $mt_val = "NA";
	my $sft_val = "NA";
	my $sft_cat = "NA";
	my $prov_val = "NA";
	my $prov_cat = "NA";
		



	#dbnsfp parallel
	my $this_chr = $fields_an[0];
	$this_chr =~ s/chr//;
	my $out_key = "dbnsfp:" .  $this_chr;
	if($retrieved_hashes{$out_key}->{$key}){
		$dbnsfp_val = $retrieved_hashes{$out_key}->{$key};
		my @dbnsfp_vals = split("," , $dbnsfp_val);
		my @dbnsfp_db_hds = split("\t" , $dbnsfp_db_header );
	    my $dbnfsp_val_index = 0;
        while( my( $hd_index, $hd_value ) = each @dbnsfp_db_hds ) {
           if($hd_index > 4) {	 # assume db header is #Chr    Start   End     Ref     Alt     SIFT_score      SIFT_pred   ...
      	     $dbnsfp_anno{ $hd_value }  = $dbnsfp_vals[$dbnfsp_val_index++];
            }
        }

        $pol_val = $dbnsfp_anno{'Polyphen2_HVAR_score'};
        $pol_val = ($pol_val eq '.') ? "NA" : $pol_val;
        
        $ma_val = $dbnsfp_anno{'MutationAssessor_score'};
        $ma_val = ($ma_val eq '.') ? "NA" : $ma_val;
        $mt_val = $dbnsfp_anno{'MutationTaster_score'};
        $mt_val = ($mt_val eq '.') ? "NA" : $mt_val;
        $sft_val = $dbnsfp_anno{'SIFT_score'};
        $sft_val = ($sft_val eq '.') ? "NA" : $sft_val;

        $sft_cat = $dbnsfp_anno{'SIFT_pred'};
        $sft_cat = ($sft_cat eq '.') ? "NA" : $sft_cat;
        
        $prov_val = $dbnsfp_anno{'PROVEAN_score'};
        $prov_val = ($prov_val eq '.') ? "NA" : $prov_val;
        $prov_cat = $dbnsfp_anno{'PROVEAN_pred'};
        $prov_cat = ($prov_cat eq '.') ? "NA" : $prov_cat;
        
		
	}	
	

#	if ( exists $retrieved_hashes{'dbnsfp'}->{$key}  ) {
#		$dbnsfp_val = $retrieved_hashes{'dbnsfp'}->{$key};
#		my @dbnsfp_vals = split("," , $dbnsfp_val);
#	    #parsing dbnsfp
#	    my @dbnsfp_db_hds = split("\t" , $dbnsfp_db_header );
#	    my $dbnfsp_val_index = 0;
#        while( my( $hd_index, $hd_value ) = each @dbnsfp_db_hds ) {
#           if($hd_index > 4) {	 # assume db header is #Chr    Start   End     Ref     Alt     SIFT_score      SIFT_pred   ...
#      	     $dbnsfp_anno{ $hd_value }  = $dbnsfp_vals[$dbnfsp_val_index++];
#            }
#        }
#        $pol_val = $dbnsfp_anno{'Polyphen2_HVAR_score'};
#        $pol_val = ($pol_val eq '.') ? "NA" : $pol_val;
#        
#        $ma_val = $dbnsfp_anno{'MutationAssessor_score'};
#        $ma_val = ($ma_val eq '.') ? "NA" : $ma_val;
#        $mt_val = $dbnsfp_anno{'MutationTaster_score'};
#        $mt_val = ($mt_val eq '.') ? "NA" : $mt_val;
#        $sft_val = $dbnsfp_anno{'SIFT_score'};
#        $sft_val = ($sft_val eq '.') ? "NA" : $sft_val;
#
#        $sft_cat = $dbnsfp_anno{'SIFT_pred'};
#        $sft_cat = ($sft_cat eq '.') ? "NA" : $sft_cat;
#        
#        $prov_val = $dbnsfp_anno{'PROVEAN_score'};
#        $prov_val = ($prov_val eq '.') ? "NA" : $prov_val;
#        $prov_cat = $dbnsfp_anno{'PROVEAN_pred'};
#        $prov_cat = ($prov_cat eq '.') ? "NA" : $prov_cat;
#        
# 	}
 	$final_out{'polyphen_score'} = $pol_val;
 	$impact_col{'poly'} = $pol_val;
    $final_out{'ma_score'} = $ma_val;
    $impact_col{'ma'} = $ma_val;
    
    $final_out{'mt_score'} = $mt_val;
    $impact_col{'mt'} = $mt_val;
     
    $final_out{'sift_score'} = $sft_val;
    $impact_col{'sft'} = $sft_val;
    
    $final_out{'PROVEAN_score'} = $prov_val;
    $stats_final_merge{counts}{'polyphen_score'} = ($final_out{'polyphen_score'} ne "NA") ? ($stats_final_merge{counts}{'polyphen_score'} + 1) : $stats_final_merge{counts}{'polyphen_score'};
    $stats_final_merge{counts}{'sift_score'} = ($final_out{'sift_score'} ne "NA") ? ($stats_final_merge{counts}{'sift_score'} + 1) : $stats_final_merge{counts}{'sift_score'};
 	
 	my $primateAI_val = "NA";
 	if ( exists $retrieved_hashes{'primateAI'}->{$key}  ) {
 		$primateAI_val = $retrieved_hashes{'primateAI'}->{$key};
 		
 	}
 	$final_out{'primateAI'} = $primateAI_val;

	#added Feb 09, 2016    dbscsnv:  dbscSNV_ADA_SCORE       dbscSNV_RF_SCORE
    my $dbscSNV_ADA_SCORE = "NA";
	my $dbscsnv_val = "NA";
	my $dbscSNV_RF_SCORE = "NA";
	if ( exists $retrieved_hashes{'dbscsnv'}->{$key}  ) {
		$dbscsnv_val = $retrieved_hashes{'dbscsnv'}->{$key};
		my @dbscsnv_vals = split("," , $dbscsnv_val);
	    #parsing dbnsfp
	    my @dbscsnv_db_hds = split("\t" , $dbscsnv_db_header );
	    my $dbscsnv_val_index = 0;
        while( my( $hd_index, $hd_value ) = each @dbscsnv_db_hds ) {
           if($hd_index > 4) {	 # assume db header is #Chr    Start   End     Ref     Alt     score      ..   ...
      	     $dbscsnv_anno{ $hd_value }  = $dbscsnv_vals[$dbscsnv_val_index++];
            }
        }
        $dbscSNV_ADA_SCORE = $dbscsnv_anno{'dbscSNV_ADA_SCORE'};
        $dbscSNV_RF_SCORE = $dbscsnv_anno{'dbscSNV_RF_SCORE'};
 	}
    $final_out{'dbscSNV_ADA_SCORE'} = $dbscSNV_ADA_SCORE;
    $final_out{'dbscSNV_RF_SCORE'} = $dbscSNV_RF_SCORE;
    $impact_col{'dbscAda'} = $dbscSNV_ADA_SCORE;
    $impact_col{'dbscRf'} = $dbscSNV_RF_SCORE;
    
    $stats_final_merge{counts}{'dbscSNV'} = ($final_out{'dbscSNV_ADA_SCORE'} ne "NA") ? ($stats_final_merge{counts}{'dbscSNV'} + 1) : $stats_final_merge{counts}{'dbscSNV'};
	
	#gnomeAD annotation
	#gnomAD_genome_ASJ
	#initialize
	#my @gnADex_db_hds = split("\t" , $gnADex_db_header );
	my %gnADex_value = ();
	foreach (@gnADex211_reqs){
		$gnADex_value{$_} = 0;
		$final_out{$_} = 0;
	}
	#my @gnADgn_db_hds = split("\t" , $gnADgn_db_header );
	my %gnADgn_value = ();
	foreach (@gnADgn30_reqs){
		$gnADgn_value{$_} = 0;
		$final_out{$_} = 0;
	}
	
	
	my %gnADgn211_value = ();
	foreach (@gnADgn211_reqs){
		$gnADgn211_value{$_} = 0;
		$final_out{$_} = 0;
	}

	if($gnAD_fix == 1){
		foreach (keys %gnADex_reqs){
			$gnADex_value{$gnADex_reqs{$_}} = 0;
			$final_out{$gnADex_reqs{$_}} = 0;
		}
		foreach (keys %gnADgn_reqs){
			$gnADgn_value{$gnADgn_reqs{$_}} = 0;
			$final_out{$gnADgn_reqs{$_}} = 0;
		}
		
	}

	$final_out{'gnomAD_genome30_FILTER'} = "NA";
	$final_out{'gnomAD_genome211_FILTER'} = "NA";
	$final_out{'gnomAD_exome211_FILTER'} = "NA";
	
	#my @gnADex_db_hds = split("\t" , $gnADex_db_header );
	#warn "** run_pp = $run_pp\n";
	if ( exists $retrieved_hashes{'gnomad_exome211'}->{$key}  ) {
		my $gnADex_val_info = $retrieved_hashes{'gnomad_exome211'}->{$key};
		#$final_out{'gnomAD_exome211_FILTER'} = $retrieved_hashes{'gnomadExV211'}->{$key}{FILTER};

		my @gnADex_vals = split(";" , $gnADex_val_info);
		foreach (@gnADex_vals){
			my ($this_gnADex_tag , $this_gnADex_val0) = split('=' , $_);
			my $this_gnADex_val = ($this_gnADex_val0 ne '.') ? $this_gnADex_val0 : 0;
			#for mssng, only selected columns
			
			if($gnAD_fix == 1){
				
			}else{
		        $final_out{'gnomAD_exome211_' . $this_gnADex_tag} = $this_gnADex_val;
		        if($this_gnADex_tag =~ m/faf95/){
					push(@gnomAD_exome_freq_max, $this_gnADex_val);
					push(@freq_max , $this_gnADex_val)
				}
		        
			}
		}
 	}else{
			push(@freq_max , 0);
			push(@gnomAD_exome_freq_max, 0);
	}


#**************** genomeAD_211*************************

	my $this_chr = $fields_an[0];
	$this_chr =~ s/chr//;
	my $out_key = "gnomad_genomeV211:" .  $this_chr;
	if(exists $retrieved_hashes{$out_key}->{$key}){
		my $gnADgn_val_info = $retrieved_hashes{$out_key}->{$key};
		#my $gnADgn_val_filter = $retrieved_hashes{$out_key}->{$key}{FILTER};
		#$final_out{'gnomAD_genome211_FILTER'} = $retrieved_hashes{$out_key}->{$key}{FILTER};
		
		my @gnADgn_vals = split(";" , $gnADgn_val_info);
#		my $gnomAA_index = 5;
		foreach (@gnADgn_vals){
			my ($this_gnADgn_tag , $this_gnADgn_val0) = split('=' , $_);
			my $this_gnADgn_val = ($this_gnADgn_val0 ne '.') ? $this_gnADgn_val0 : 0;
			if($gnAD_fix == 1){
#				my $this_gnHd = $gnADgn_db_hds[$gnomAA_index];
#				if(exists $gnADgn_reqs{$this_gnHd}){
#					my $this_gnHdOld = $gnADgn_reqs{$this_gnHd};
#					$final_out{$this_gnHdOld} = $this_gnADgn_val ;
					#push(@freq_max , 0);
					#push(@gnomAD_genome_freq_max, 0);
#					$stats_final_merge{counts}{$this_gnHdOld} = ($final_out{$this_gnHdOld} != 0) ? ($stats_final_merge{counts}{$this_gnHdOld} + 1) : $stats_final_merge{counts}{$this_gnHdOld};
#					
#				}
#				$gnomAA_index++;

			}else{
				$final_out{'gnomAD_genome211_' . $this_gnADgn_tag} = $this_gnADgn_val;
				#push(@freq_max , 0);
				if($this_gnADgn_tag =~ m/faf95/){
					#push(@gnomAD_genome_freq_max, $this_gnADgn_val);
					#push(@freq_max , $this_gnADgn_val);
				}
				#$stats_final_merge{counts}{$this_gnAD_col} = ($final_out{$this_gnAD_col} != 0) ? ($stats_final_merge{counts}{$this_gnAD_col} + 1) : $stats_final_merge{counts}{$this_gnAD_col};
			}
			
		}
		
		
	}else{
			#push(@freq_max , 0);
			push(@gnomAD_genome_freq_max, 0);
 	}


#***********************************************************************	

	#gnomAD genome parallel
#	$dbval{$key}{INFO} = $fields[8];
#	$dbval{$key}{FILTER} = $fields[7];
		
	my $this_chr = $fields_an[0];
	$this_chr =~ s/chr//;
	my $out_key = "gnomad_genomeV30:" .  $this_chr;
	if(exists $retrieved_hashes{$out_key}->{$key}){
		my $gnADgn_val_info = $retrieved_hashes{$out_key}->{$key};
		#my $gnADgn_val_filter = $retrieved_hashes{$out_key}->{$key}{FILTER};
		#$final_out{'gnomAD_genome30_FILTER'} = $retrieved_hashes{$out_key}->{$key}{FILTER};
		
		my @gnADgn_vals = split(";" , $gnADgn_val_info);
#		my $gnomAA_index = 5;
		foreach (@gnADgn_vals){
			my ($this_gnADgn_tag , $this_gnADgn_val0) = split('=' , $_);
			my $this_gnADgn_val = ($this_gnADgn_val0 ne '.') ? $this_gnADgn_val0 : 0;
			if($gnAD_fix == 1){
#				my $this_gnHd = $gnADgn_db_hds[$gnomAA_index];
#				if(exists $gnADgn_reqs{$this_gnHd}){
#					my $this_gnHdOld = $gnADgn_reqs{$this_gnHd};
#					$final_out{$this_gnHdOld} = $this_gnADgn_val ;
					#push(@freq_max , 0);
					#push(@gnomAD_genome_freq_max, 0);
#					$stats_final_merge{counts}{$this_gnHdOld} = ($final_out{$this_gnHdOld} != 0) ? ($stats_final_merge{counts}{$this_gnHdOld} + 1) : $stats_final_merge{counts}{$this_gnHdOld};
#					
#				}
#				$gnomAA_index++;

			}else{
				$final_out{'gnomAD_genome30_' . $this_gnADgn_tag} = $this_gnADgn_val;
				#push(@freq_max , 0);
				if($this_gnADgn_tag =~ m/faf95/){
					push(@gnomAD_genome_freq_max, $this_gnADgn_val) unless ($this_gnADgn_tag eq 'faf95_adj') ;
					push(@freq_max , $this_gnADgn_val) unless ($this_gnADgn_tag eq 'faf95_adj');
				}
				#$stats_final_merge{counts}{$this_gnAD_col} = ($final_out{$this_gnAD_col} != 0) ? ($stats_final_merge{counts}{$this_gnAD_col} + 1) : $stats_final_merge{counts}{$this_gnAD_col};
			}
			
		}
		
		
	}else{
			push(@freq_max , 0);
			push(@gnomAD_genome_freq_max, 0);
 	}
	
	################
	
#	if ( exists $retrieved_hashes{'gnomad_genome'}->{$key}  ) {
#		my $gnADgn_val = $retrieved_hashes{'gnomad_genome'}->{$key};
#		my @gnADgn_vals = split("," , $gnADgn_val);
#		my $gnomAA_index = 5;
#		foreach (@gnADgn_vals){
#		#	$gnADgn_value{$gnADgn_db_hds[$gnomAA_index++]} = $_;
#			my $this_gnADgn_val = ($_ ne '.') ? $_ : 0;
#			if($run_pp == 5){
#				my $this_gnHd = $gnADgn_db_hds[$gnomAA_index];
#				if(exists $gnADgn_reqs{$this_gnHd}){
#					my $this_gnHdOld = $gnADgn_reqs{$this_gnHd};
#					$final_out{$this_gnHdOld} = $this_gnADgn_val ;
#					push(@freq_max , $this_gnADgn_val) unless ($this_gnHdOld eq "gnomAD_genome_ASJ");
#					push(@gnomAD_genome_freq_max, $this_gnADgn_val) unless ($this_gnHdOld eq "gnomAD_genome_ASJ");
#					$stats_final_merge{counts}{$this_gnHdOld} = ($final_out{$this_gnHdOld} != 0) ? ($stats_final_merge{counts}{$this_gnHdOld} + 1) : $stats_final_merge{counts}{$this_gnHdOld};
#				}
#				$gnomAA_index++; 
#			}else{			
#				$final_out{$gnADgn_db_hds[$gnomAA_index++]} = $this_gnADgn_val; # fixed for . , sep 28,2017
#				push(@freq_max , $this_gnADgn_val);
#				push(@gnomAD_genome_freq_max, $this_gnADgn_val);
#				$stats_final_merge{counts}{$gnADgn_db_hds[$gnomAA_index-1]} = ($final_out{$gnADgn_db_hds[$gnomAA_index-1]} != 0) ? ($stats_final_merge{counts}{$gnADgn_db_hds[$gnomAA_index-1]} + 1) : $stats_final_merge{counts}{$gnADgn_db_hds[$gnomAA_index-1]};
#			}
#		}
# 	}else{
#			push(@freq_max , 0);
#			push(@gnomAD_genome_freq_max, 0);
# 	}
#	my @gnADgn_outs = ();
#	foreach (@gnADgn_hds){
#		push(@gnADgn_outs , $gnADgn_value{$_});
#	}

	
	#adding max frequecies
	
	
	$final_out{'freq_max'} = max(@freq_max);
	$final_out{'A1000g_freq_max'} = max(@x1000g_freq_max);
	$final_out{'ExAC_freq_max'} = max(@ExAC_freq_max);
	$final_out{'cg_freq_max'} = max(@cg_freq_max);
	$final_out{'gnomAD_exome_freq_max'} = max(@gnomAD_exome_freq_max);
	
	$final_out{'gnomAD_genome_freq_max'} = max(@gnomAD_genome_freq_max);

	$final_out{'cosmic'} = ( exists  $retrieved_hashes{'cosmic'}->{$key} ) ? $retrieved_hashes{'cosmic'}->{$key} : "NA";	
	$stats_final_merge{counts}{'cosmic'} = ($final_out{'cosmic'} ne "NA") ? ($stats_final_merge{counts}{'cosmic'} + 1) : $stats_final_merge{counts}{'cosmic'};
	
#********** has to do different way
#	my $dbr1_val = "NA";
#	if ( exists $retrieved_hashes{'dbsnpR'}->{$key} ) {
#		$dbr1_val = $retrieved_hashes{'dbsnpR'}->{$key};
#		$dbr1_val =~ s/Name=//;
#	}
#	$final_out{'dbsnp_region'} = $dbr1_val;
#	my $dbrc_val = "NA";
#	
##	if ( exists $retrieved_hashes{'dbsnpc region'}->{$key} ) {
##		$dbrc_val = $retrieved_hashes{'dbsnpc region'}->{$key};
##		$dbrc_val =~ s/Name=//;
##	}
##	$final_out{'dbsnp_common_region'} = $dbrc_val;
#	
#	my $dbr2_val = "NA";
#	if ( exists $retrieved_hashes{'dbsnpW'}->{$key} ) {
#		$dbr2_val = $retrieved_hashes{'dbsnpW'}->{$key};
#		$dbr2_val =~ s/Name=//;
#	}
#	$final_out{'dbsnp_wind'} = $dbr2_val;


# updated on Nov 5,2018

    my $dbr1_val = "NA";
	$this_chr = $fields_an[0];
	$this_chr =~ s/chr//;
	my $out_key = "dbsnpR:" .  $this_chr;
	#warn "####gerpwgs: $out_key\n";
	if($retrieved_hashes{$out_key}->{$key}){
		$dbr1_val = $retrieved_hashes{$out_key}->{$key};
		$dbr1_val =~ s/Name=//;
	}
	$final_out{'dbsnp_region'} = $dbr1_val;
	$impact_col{'dbsnpr'} = $dbr1_val;
	$stats_final_merge{counts}{'dbsnp_region'} = ($final_out{'dbsnp_region'} ne "NA") ? ($stats_final_merge{counts}{'dbsnp_region'} + 1) : $stats_final_merge{counts}{'dbsnp_region'};

    my $dbr2_val = "NA";
	#my $this_chr = $fields_an[0];
	#$this_chr =~ s/chr//;
	my $out_key = "dbsnpW:" .  $this_chr;
	#warn "####gerpwgs: $out_key\n";
	if($retrieved_hashes{$out_key}->{$key}){
		$dbr2_val = $retrieved_hashes{$out_key}->{$key};
		$dbr2_val =~ s/Name=//;
	}
	$final_out{'dbsnp_wind'} = $dbr2_val;
	$stats_final_merge{counts}{'dbsnp_wind'} = ($final_out{'dbsnp_wind'} ne "NA") ? ($stats_final_merge{counts}{'dbsnp_wind'} + 1) : $stats_final_merge{counts}{'dbsnp_wind'};

#

	my $pfam_annovar = "NA";
	#pfam is only for exionic-type variants , fixed Feb 13, 2016
	if ( exists $retrieved_hashes{'pfam'}->{$key} ) {
		$pfam_annovar = $retrieved_hashes{'pfam'}->{$key};
		$pfam_annovar =~ s/Name=//;
	}
	$final_out{'pfam_annovar'} = $pfam_annovar;

	#$final_out{'gerp_elem'} = ( exists  $retrieved_hashes{'gerp_elem'}->{$key} ) ? $retrieved_hashes{'gerp_elem'}->{$key} : "NA";

	my $gerp_elem_annovar = "NA";
	if ( exists $retrieved_hashes{'gerp_elem'}->{$key} ) {
		$gerp_elem_annovar = $retrieved_hashes{'gerp_elem'}->{$key};
		$gerp_elem_annovar =~ s/Name=//;
	}
	$final_out{'gerp_elem'} = $gerp_elem_annovar;


 #*************************

	my $hgm_value = ( exists $hgm{$key} ) ? $hgm{$key} : "NA";
	my $hgm2_value = ( exists $retrieved_hashes{'hgmd_pro'}->{$key} ) ? $retrieved_hashes{'hgmd_pro'}->{$key} : "NA\@NA\@NA\@NA\@NA\@NA";
	$final_out{'HGMD_type'} = (split("\@", $hgm2_value ))[0];
	#my $hgmd2_class = (split("\@", $hgm2_value))[1];
	$final_out{'HGMD_tag'} = (split("\@", $hgm2_value))[1];
	#my $hgmd2_acc = (split("\@", $hgm2_value))[2];
	$final_out{'HGMD_Accession'} = (split("\@", $hgm2_value))[2];
	$final_out{'HGMD_Disease'}  = (split("\@", $hgm2_value))[3];
	$final_out{'HGMD_PubmedId'} = (split("\@", $hgm2_value))[4];
	$final_out{'HGMD_RANKSCORE'} = (split("\@", $hgm2_value))[5];
	
	my $cadd_value = ( exists $retrieved_hashes{'cadd'}->{$key} ) ? $retrieved_hashes{'cadd'}->{$key} : "NA,NA";  # since v2.6
	($final_out{'CADD_Raw'}, $final_out{'CADD_phred'}) = split("," ,  $cadd_value);
	$impact_col{'cadd'} = $final_out{'CADD_phred'};
	$stats_final_merge{counts}{'CADD_phred'} = ($final_out{'CADD_phred'} ne "NA") ? ($stats_final_merge{counts}{'CADD_phred'} + 1) : $stats_final_merge{counts}{'CADD_phred'};

	#added July 20,2016
	my $spidex_value = ( exists $retrieved_hashes{'SPIDEX'}->{$key} ) ? $retrieved_hashes{'SPIDEX'}->{$key} : join(":" , @spx_valna);  # since Jul2016
	my @spidex_values = split(":" , $spidex_value);  # since Jul2016
	my %spidex2go = ();
	my $spx_index = 0;
	foreach my $this_spx_hd  (@spx_hds){
		#$spidex2go{$this_spx_hd} = $spidex_values[$spx_index++];
		$final_out{$this_spx_hd} = $spidex_values[$spx_index++];
	}
	$impact_col{'spx'} = $final_out{'spx_dpsi'};
	

#	my $spliceAIsnv_value = ( exists $retrieved_hashes{'spliceAIsnv'}->{$key} ) ? $retrieved_hashes{'spliceAIsnv'}->{$key} : join("|" , @spliceAIsnv_valna);  #
#	my @spliceAIsnv_values = split(';' , $spliceAIsnv_value);
#	my $spliceAI_index = 0;
#	foreach my $this_spliceAIhd (@spliceAIsnv_hds){
#		$final_out{$this_spliceAIhd} = $spliceAIsnv_values[$spliceAI_index++];
#	}

	#SpliceAI=G|OR4F5|0.01|0.08|0.00|0.00|-10|26|-28|-25
	#ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL
	if(exists $retrieved_hashes{'spliceAIsnv'}->{$key}){
		my $spliceAIsnv_value =  $retrieved_hashes{'spliceAIsnv'}->{$key};
		my ($spliceAI_1, $spliceAI_2) = split('=' , $spliceAIsnv_value);
		my $spliceAI_index = 0;
		foreach my $this_spliceAI_val (split('\|' , $spliceAI_2)){
			$final_out{$spliceAIsnv_hds[$spliceAI_index++]} = $this_spliceAI_val;
		}
	}else{
		my $spliceAI_index = 0;
		foreach my $this_spliceAI_val (@spliceAIsnv_valna){
			$final_out{$spliceAIsnv_hds[$spliceAI_index++]} = $this_spliceAI_val;
		}
	}	
	
	
	
	
	#SYMBOL=ARHGEF16;STRAND=+;TYPE=E;DIST=-51;DS_AG=0.0014;DS_AL=0.0000;DS_DG=0.4581;DS_DL=0.0000;DP_AG=-47;DP_AL=-2;DP_DG=-3;DP_DL=14
	
	 
#	my @spidex_values = split(":" , $spidex_value);  # since Jul2016
#	my %spidex2go = ();
#	my $spx_index = 0;
#	foreach my $this_spx_hd  (@spx_hds){
#		#$spidex2go{$this_spx_hd} = $spidex_values[$spx_index++];
#		$final_out{$this_spx_hd} = $spidex_values[$spx_index++];
#	}


	my $rep_value = ( exists $retrieved_hashes{'repeat'}->{$key} ) ? $retrieved_hashes{'repeat'}->{$key} : "NA";
	$rep_value =~ s/Name=//;
	$final_out{'Repeat'} = $rep_value;
	$stats_final_merge{counts}{'Repeat'} = ($final_out{'Repeat'} ne "NA") ? ($stats_final_merge{counts}{'Repeat'} + 1) : $stats_final_merge{counts}{'Repeat'};
	
#	my $exb_distance = ( exists $retrieved_hashes{'exome_boundary_dist'}->{$key} ) ? $retrieved_hashes{'exome_boundary_dist'}->{$key} : "NA";
#	$final_out{'distance_spliceJunction'} = $exb_distance;
	if($run_pp != 5){
		$final_out{'annovar_chr'} = $fields_an[0];
		$final_out{'annovar_start'} = $fields_an[1];
		$final_out{'annovar_end'} = $fields_an[2];
		$final_out{'ref_allele'} = $fields_an[3];
		$final_out{'alt_allele'} = $fields_an[4];
	}	
	
 #***********************************
	
	my $seq_overlap_type = "NA";
	my $gene_symbol      = "N_A";
	my $ex_rfq = "NA";
	my $ex_eff = "NA";
	my $eff_category = "NA";
	my $eff_category_pri = "NA";
	my $refSeq_id        = "NA";
	my $refSeq_id_sorted = "NA";
	my $cds_prop                = "NA";
	my $each_entz_id = "NA";
	my $per_transcript_affected = "NA";
	my $aa_flag = "NA";
	my ($exon_dist_valueL,$exon_dist_valueR)  = ("NA", "NA");
	my $gene_iter = 0;
	my %for_gene_iter = ();
	my $typeseq_iter = 0;
	my $exb_distance = "NA";
	
	#main iteration, for each seqOverlap and each genes
	#Gene-Select
	
	my $rt_gene_h = $retrieved_hashes{'Gene-var'};
	my $rt_genex_h = $retrieved_hashes{'Gene-exonic'};
	my $rt_geneucx_h = $retrieved_hashes{'Gene-Rlx-exonic'};
	my $rt_geneapr_h = $retrieved_hashes{'Gene-Appris'};
	my $rt_geneaprx_h = $retrieved_hashes{'Gene-Appris-exonic'};
	my $rt_genesel_h = $retrieved_hashes{'Gene-Select'};
	my $rt_geneselx_h = $retrieved_hashes{'ene-Select-exonic'};
	
	my $this_geneExtra_mod = "gene_extra:" . $fields_an[0];
	my $rt_geneExtraNew = $retrieved_hashes{$this_geneExtra_mod};
	
	##gencode-Gene-var
	
	#print "RT VAL = $rt_gene_h $key\n"; 
	
	my %var_outs = ();
	my %var_seqOverlaps = ();
	my %var_effects = ();
	my %primary_var_seqOverlaps = ();
	my %refseqid_for_effect = ();
	my $pri_refseq = "NA";
	my $pri_refSeq_id_sorted = "NA";
	
	my $res2go_left = -1;
	my $res2go_right = -1;
	my $primary_seqOverlap = "NA";
	my %this_vcfout = ();
	my @vcfout_field = ();
	my @vcfout_field_beg = ();
	my $vcfout_info_field = ".";
	if($vcfout == 1){
		@vcfout_field_beg = @vcf_fields[0..6];
		$vcfout_info_field = $vcf_fields[7];
	}
	
	if ( exists $rt_geneExtraNew->{$key} ) {

		my @gvar_kys = keys %{$rt_geneExtraNew->{$key}};
		
		foreach my $this_gn (sort {$b cmp $a} @gvar_kys){			 # each genes
		 	my @res_field = ();
			if($run_pp == 5){
			  #no prefixes
			  
 			}elsif($run_type eq "v"){
    				push(@res_field, @vcf_fields);
 			}elsif($run_type eq "t"){
    				push(@res_field, @text_fields2go);
   			}else{
				my $line_info =
	    			$fields_an[5] . "\t"
	  				. $fields_an[6] . "\t"
	  				. $fields_an[7];
   				
     			push(@res_field, $line_info);
   			}

#		$gva{$key}->{$gn}->{'cds_prop' } = $fields[8];
#		$gva{$key}->{$gn}->{'tr_prop'} = $fields[9];
#		$gva{$key}->{$gn}->{'typeseq_select'} = $fields[17];
#		$gva{$key}->{$gn}->{'effect_select'} = $fields[18];
   			
   			$final_out{'typeseq'} = $rt_geneExtraNew->{$key}->{$this_gn}->{'typeseq'};
   			$final_out{'typeseq_priority'} = $rt_geneExtraNew->{$key}->{$this_gn}->{'typeseq_pr'};
   			$impact_col{'typeseq'} = $rt_geneExtraNew->{$key}->{$this_gn}->{'typeseq_pr'};
			$final_out{'refseq_id'} = $rt_geneExtraNew->{$key}->{$this_gn}->{'refseq-id'};
			$final_out{'aa_flag'} = $rt_geneExtraNew->{$key}->{$this_gn}->{'aa_flag'};
			$final_out{'effect'} = $rt_geneExtraNew->{$key}->{$this_gn}->{'effect'};
			$final_out{'effect_priority'} = $rt_geneExtraNew->{$key}->{$this_gn}->{'effect_pr'};
			$final_out{'effect_appris'} = $rt_geneExtraNew->{$key}->{$this_gn}->{'effect_appris'};
			
			#get_uniqValueComma
			$final_out{'typeseq_RefseqSelect'} = get_uniqValueComma($rt_geneExtraNew->{$key}->{$this_gn}->{'typeseq_select'});
			$final_out{'effect_RefseqSelect'} = get_uniqValueComma($rt_geneExtraNew->{$key}->{$this_gn}->{'effect_select'});
			
			#splice_dist
			$final_out{'distance_spliceJunction'} = $rt_geneExtraNew->{$key}->{$this_gn}->{'splice_dist'};

			$impact_col{'eff'} = $rt_geneExtraNew->{$key}->{$this_gn}->{'effect_pr'};

			$final_out{'gene_symbol'} = $rt_geneExtraNew->{$key}->{$this_gn}->{'gene'};
			$each_entz_id = $rt_geneExtraNew->{$key}->{$this_gn}->{'entrez_id'};
			$final_out{'entrez_id'} = $each_entz_id;
			$final_out{'gene_desc'} = $rt_geneExtraNew->{$key}->{$this_gn}->{'gene_desc'};
			$final_out{'gene_type'} = $rt_geneExtraNew->{$key}->{$this_gn}->{'gene_type'};
			$final_out{'omim_id'} = $rt_geneExtraNew->{$key}->{$this_gn}->{'omim_id'};
			$final_out{'omim_phenotype'} = $rt_geneExtraNew->{$key}->{$this_gn}->{'omim_phenotype'};
		    my $current_mpo = ($each_entz_id eq "" || $gene_mpo{$each_entz_id} eq "") ? "NA" : $gene_mpo{$each_entz_id};
        	my $current_hpo = ($each_entz_id eq "" || $gene_hpo{$each_entz_id} eq "" ) ? "NA" : $gene_hpo{$each_entz_id};
        	my $current_cgd = ($each_entz_id eq "" || $gene_cgd{$each_entz_id} eq "") ? "NA\@NA" : $gene_cgd{$each_entz_id};
        	my ($current_cgd_dis, $current_cgd_inh) = split('@', $current_cgd);
			
			$final_out{'MPO'} = $current_mpo;
			$final_out{'HPO'} = $current_hpo;
			$final_out{'CGD_disease'} = $current_cgd_dis;
			$final_out{'CGD_inheritance'} = $current_cgd_inh;
			#ACMG_disease
			$final_out{'ACMG_disease'} = ($each_entz_id eq "" || $acmg{$each_entz_id} eq "" ) ? "NA" : $acmg{$each_entz_id};
			
			#cds_prop
			$final_out{'per_cds_affected'} = $rt_geneExtraNew->{$key}->{$this_gn}->{'cds_prop'};
			$final_out{'per_transcripts_affected'} = $rt_geneExtraNew->{$key}->{$this_gn}->{'tr_prop'};
			

			foreach my $xm_hd (@xma__required){
				my $this_xmt_key = "gnomAD_" . $xm_hd;
				$final_out{$this_xmt_key} = ($each_entz_id eq "" || $gene_exmat{$each_entz_id}{$xm_hd} eq "" ) ? "NA" : $gene_exmat{$each_entz_id}{$xm_hd};;
			}

			#calculate impact
			my $s_bit = 0;
			my $s_bit_str = "NA";
			if($run_type eq "t"){
	
				my $this_thr =  threshold($sft_val, $pol_val, $ma_val , $mt_val, $phy_avg_val, $phypr_avg_val, $impact_col{'cadd'});
				$impact_col{'thr'} = $this_thr;
				my $this_effs = get_effects(\%impact_col);
				foreach  (@{ $this_effs } ){
					my $s1 = (exists $EFFECTS{$_}) ? (2 ** ($EFFECTS{$_} + bitwise_impact_of_effect($_, \%impact_col)) ) : 0;
					$s_bit += $s1;
				}
				$s_bit_str = decode_bit($s_bit);
			}

				

			
		 	foreach my $cshw (@cols_show){
		 		if(exists $final_out{$cshw} &&  $final_out{$cshw} ne ""){
		 			my $final_out_val = $final_out{$cshw};
		 			push(@res_field , $final_out_val);
					if( $final_out_val eq "NA" || $final_out_val eq "0"){}else{
						$validated{$cshw} = 1;
					}
					

		 			
		 		}else{
					die "**Critical error! Column $cshw has not been populated: ( variant $key )\n";
				}
				
				
		 	} #my $cshw (@cols_show)
			push(@res_field ,  $s_bit  );
			my $anno_comments2go = (scalar (@anno_comments) == 0) ? "NA" : join("," , @anno_comments);
			push(@res_field , $anno_comments2go);
		 	
			print  join("\t" , @res_field) , "\n";
			$tot_var_with_genes++;
		 	
		}   #each genes
	}

###############################################################

	
	#warn "Refseq after loop: $refSeq_id\n";
#	my %this_vcfout = ();
#	my @vcfout_field = ();
#	my @vcfout_field_beg = ();
#	my $vcfout_info_field = ".";

	 
	if(scalar(keys %for_gene_iter)  > 1){
		$multi_gene++;
	}
	
	if($vcfout == 1){
		foreach my $cshw (@cols_show){
			if(exists $this_vcfout{$cshw}){
				my $vcf_val2go = get_uniqValue($this_vcfout{$cshw});
				push(@vcfout_field , $cshw . "=" . $vcf_val2go );
			}
		}
	}
	
	#print OVCF $ovcf_value , "\t" , join(";" , @vcfout_field ) , "\n";
	if($vcfout == 1){ 
		print $fh_outvcf join("\t" , @vcfout_field_beg) , "\t" , ($vcfout_info_field . ";" . join(";" , @vcfout_field)) , "\t" , $vcf_fields[8] , "\t" , $vcf_fields[9] , "\n";
	}
	$stats_final_merge{counts}{'multi-gene-overlap'} = $multi_gene ;
	
	if($typeseq_iter > 1){
		$multi_typeseq++;
	}
	$stats_final_merge{counts}{'multi-sequenceOverlap'} = $multi_typeseq ;
	$stats_final_merge{counts}{'multi-sequenceOverlapPerGene'} = $multi_typeseq_per_gene ;
	
	if($iter >= $unit_stat){
		my $pg_string = print_progress($extra_count, $tot_var, $iter);
		$iter = 0;
		warn $pg_string , "\n";
	}
	$iter++;

}
close SAM;
#close $fh_outvcf;
my $gene_stat_file = $pipeline_work_dir . "/" . $annovar_in_basename . "_gene_var_stats.tsv";
print_gene_stats($gene_stat_file);

#close WAR;

my $var_count_diff = $input_var_count - $tot_var;
if($var_count_diff != 0){
	usage("Critical: inconsistent input/output counts: \(input = $input_var_count vs output =  $tot_var \)");
}

print $fh_vst "Total : $tot_var \n";
print $fh_vst "Diff : $var_count_diff\n";
print $fh_vst "Total output rows : $tot_var_with_genes\n";

my $end_time = localtime;
my $done_time = $end_time - $start_time;
warn "\[$end_time\] Done ( " . sprintf("%.2f" ,  $$done_time/60)  .  " min ) , Total $tot_var entries processed. $tot_var_with_genes rows output.\n";



#print $fh_vst "MultiGene : $multi_gene \n";
#print $fh_vst "MultiSeqOverlap : $multi_typeseq \n";
#print $fh_vst "\nExonic : " , $full_summary{exonic}{'count'} , "\n";
##print VST "\nexonic;splicing : " , $full_summary{'exonic;splicing'}{'count'} , "\n";
#print $fh_vst "Intronic : " , $full_summary{intronic}{'count'} , "\n";
#print $fh_vst "UTR3 : " , $full_summary{UTR3}{'count'} , "\n";
#print $fh_vst "UTR5 : " , $full_summary{UTR5}{'count'} , "\n";
#print $fh_vst "Splicing : " , $full_summary{splicing}{'count'} , "\n";
#print $fh_vst "Upstream : " , $full_summary{upstream}{'count'} , "\n";
#print $fh_vst "Downstream : " , $full_summary{downstream}{'count'} , "\n";
#print $fh_vst "Intergenic : " , $full_summary{intergenic}{'count'} , "\n";
##refseq_UKNOWN
#print $fh_vst "refseq_Unknown_effect : " , $full_summary{refseq_UKNOWN_effect}{'count'} , "\n\n";

#Generating gene stats, Jan 14,2020
warn "\nGenerating Gene-level stats...\n";



my %gene_extra_mat = ();
warn join("\n" , @gene_extra_stats) , "\n";
foreach my $this_st (@gene_extra_stats){
	my $fh_gst = IO::File->new($this_st , q{<}) or die "$! $this_st\n";
	while(my $line = <$fh_gst>){
		chomp($line);
		my ($item, $gst) = split('=' , $line);
		$item =~ s/^\s+//; #remove leading spaces
		$item =~ s/\s+$//; #remove trailing spaces
		$gst =~ s/^\s+//; #remove leading spaces
		$gst =~ s/\s+$//; #remove trailing spaces
		my ($item_type, $item_name) = split('->' , $item);
		$gene_extra_mat{$item_type}->{$item_name} = (exists $gene_extra_mat{$item_type}->{$item_name}) ? ($gene_extra_mat{$item_type}->{$item_name} + $gst) : $gst; 
	}
}

foreach my $ck (keys %{$gene_extra_mat{'General'}}  ){
	print $fh_vst $ck , " : " , $gene_extra_mat{'General'}{$ck} , "\n";
}
foreach my $ck (keys %{$gene_extra_mat{'Typeseq'}}  ){
	print $fh_vst $ck , " : " , $gene_extra_mat{'Typeseq'}{$ck} , "\n";
}
foreach my $ck (keys %{$gene_extra_mat{'Effect'}}  ){
	print $fh_vst $ck , " : " , $gene_extra_mat{'Effect'}{$ck} , "\n";
}


#writing new stat file, oct16,2017
print $fh_mrg_mta Dumper \%stats_final_merge , "\n";
#print $fh_mrg_mta Dumper \%gene_stats , "\n";
close $fh_mrg_mta;

#%{$rt_gene_h->{$key}}
foreach (sort keys %{$full_summary{'exonic'}}){
	print $fh_vst $_ , ":\t" , $full_summary{'exonic'}{$_} , "\n" unless ($_ eq "count");
}

#print VST "Gene without entrez-id :  $no_entrez_id_genes\n";
if($no_entrez_id_genes != 0	){
print $fh_vst
"\nThese genes has no mapping to entrez-id ( $no_entrez_id_genes variants ) :-> ",
  join( ",", ( keys %no_entrez_id_genes_list ) ), "\n";
}

print $fh_vst "\nDetailed Stats:\n";
print $fh_vst "*****************\n";
print $fh_vst Dumper \%stats_final_merge , "\n";

print $fh_vst "\n>>Output columns validation(based on full set, vcf specific columns are excluded): \n";
foreach (@cols_show){
	if(not exists $validated{$_}){
		print $fh_vst "Column $_ has only NAs or 0s\n";
	}
}

close $fh_vst;



#### *** subroutines ****************



sub load_input_hash{
	my ($mm, $rn_ty) = @_;
	my %dbval = ();
	$dbval{'in_type'} = $rn_ty;
	if($rn_ty eq "m"){
	#open(IX, $mm) or warn "$mm $!";
	my $fh_ix = IO::File->new($mm , q{<}) or warn "$mm  $!\n";
	while(my $line = <$fh_ix>){
		chomp($line);
		if ( $line =~ m/^#.+/ || $line =~ m/^>.+/ || $line eq "" ) {
			
		}else{
			my $locus_id = ( split( "\t", $line ) )[0];
#			$locus_id =~ s/^\s+//;    #remove leading spaces		$gva{$key}->{$gn}->{'cds_prop' } = $fields[8];

#			$locus_id =~ s/\s+$//;    #remove trailing spaces
			$locus_id =~ s/^\s*(\S*(?:\s+\S+)*)\s*$/$1/;
			$dbval{$locus_id} = $line;
			
		}
	}
	close $fh_ix;
	}
	return \%dbval;
}


sub load_hash{
	my ($mm) = @_;
	my %dbval = ();
	#open(X, $mm) or warn "File(load_hash): $mm $!";
	#warn "**Load hash: $mm\n";
	my $fh_x = IO::File->new($mm , q{<}) or warn "$mm  $!\n";
	if(defined $fh_x){
		while(my $line = <$fh_x>){
			chomp($line);
			my @fields = split("\t", $line);
			my $key = $fields[2] . ":" . $fields[3] . ":" . $fields[4] . ":" . $fields[5] . ":" . $fields[6];
			$dbval{$key} = $fields[1];
		}
		close $fh_x;
	}
	return \%dbval;
}

#1       chr1
#2       10146
#3       rs779258992
#4       AC
#5       A
#6       5648728.00
#7       PASS
#8       AC=9915;AN=15668;AF=6.32818e-01;AC_female=5302;AN_female=8372;AF_female=6.33301e-01;nhomalt_female=1365;AF_ami=5.57143e-01;AF_oth=6.29771e-01;AC_male=4613;AN_male=7296;AF_male=6.32264e-01;nhomalt_male=1199;AF_afr=6.29990e-01;AF_sas=6.65698e-01;AF_raw=6.10203e-01;AF_asj=6.45923e-01;nhomalt=2564;AF_fin=5.20833e-01;AF_amr=6.38009e-01;AF_eas=5.98802e-01;faf95_afr=6.09760e-01;faf95_sas=5.95029e-01;faf95_adj=6.22401e-01;faf95_amr=6.02362e-01;faf95_nfe=6.26793e-01;faf95_eas=5.30899e-01

sub load_vcfanno{
	my ($mm) = @_;
	my %dbval = ();
	#open(X, $mm) or warn "File(load_hash): $mm $!";
	
	my $fh_x = IO::File->new($mm , q{<}) or warn "$mm  $!\n";
	if(defined $fh_x){
		while(my $line = <$fh_x>){
			chomp($line);
			my @fields = split("\t", $line);
			my $key = $fields[0] . ":" . $fields[1] . ":" . $fields[3] . ":" . $fields[4];
			$dbval{$key}{INFO} = $fields[7];
			$dbval{$key}{FILTER} = $fields[6];
		}
		close $fh_x;
	}
	return \%dbval;
}

sub load_cgInt{
	my ($mm) = @_;
	my %dbval = ();
	#open(X, $mm) or warn "File(load_hash): $mm $!";
	
	my $fh_x = IO::File->new($mm , q{<}) or warn "$mm  $!\n";
	if(defined $fh_x){
		while(my $line = <$fh_x>){
			chomp($line);
			my @fields = split("\t", $line);
			my @val = @fields[2 .. $#fields];
			my $key = $fields[0];
			$dbval{$key} = join("\t" , @val );
		}
		close $fh_x;
	}
	return \%dbval;
}

sub load_vcfannoText{
	my ($mm) = @_;
	my %dbval = ();
	#open(X, $mm) or warn "File(load_hash): $mm $!";
	
	my $fh_x = IO::File->new($mm , q{<}) or warn "$mm  $!\n";
	if(defined $fh_x){
		while(my $line = <$fh_x>){
			chomp($line);
			my @fields = split("\t", $line);
			my $key = $fields[0];
			$dbval{$key}{INFO} = $fields[8];
			$dbval{$key}{FILTER} = $fields[7];
		}
		close $fh_x;
	}
	return \%dbval;
}


sub load_hashW{
	my ($mm) = @_;
	my %dbval = ();
#	open(X, $mm) or warn "File: $mm $!";
	my $fh_x = IO::File->new($mm , q{<}) or warn "$mm  $!\n";
	while(my $line = <$fh_x>){
		chomp($line);
		my @fields = split("\t", $line);
		my $key = $fields[2] . ":" . $fields[3]  . ":" . $fields[4]; #":" . $fields[5] . ":" . $fields[6];
		$dbval{$key} = $fields[1];
	}
	close $fh_x;
	return \%dbval;
}

sub load_genev_hash{
	my ($mm) = @_;
	my %gva = ();
	#open(GV, $mm) or warn "$mm $!";
	my $fh_gv = IO::File->new($mm , q{<}) or warn "$mm  $!\n";
	while(my $line = <$fh_gv>){
		chomp($line);
		my @fields = split("\t", $line);
		my $key = $fields[2] . ":" . $fields[3] . ":" . $fields[4] . ":" . $fields[5] . ":" . $fields[6];
		my $gn = Parse::Annovar::get_genename( $fields[0], $fields[1] );
#		warn "Gene = $gn\n";
		$gva{$key}->{ $fields[0] } = $gn;
		#$gva{$key}->{'refseqid'} = $fields[1] ;
		$gva{$key}->{'refseqid'} = (exists $gva{$key}->{'refseqid'} ) ?  ($gva{$key}->{'refseqid'} . ";" . $fields[1]) : $fields[1] ;
	}
	close $fh_gv;
	return \%gva;
}

sub load_genex_hash{
	my ($mm) = @_;
	my %gx = ();
#	open(GX, $mm) or warn "$mm $!";
	my $fh_gx = IO::File->new($mm , q{<}) or warn "$mm  $!\n";
	while(my $line = <$fh_gx>){
		chomp($line);
		my @fields = split("\t", $line);
		my $key = $fields[3] . ":" . $fields[4] . ":" . $fields[5] . ":" . $fields[6] . ":" . $fields[7];
		foreach my $each_gex (split("," , $fields[2] )){
			my $this_gex_gene = (split(":", $each_gex))[0];
			$gx{$key}->{$this_gex_gene}->{rfq} = (exists $gx{$key}->{$this_gex_gene}->{rfq}) ?  $gx{$key}->{$this_gex_gene}->{rfq} . "," .   $each_gex : $each_gex; 
			$gx{$key}->{$this_gex_gene}->{eff} = (exists $gx{$key}->{$this_gex_gene}->{eff}) ?  $gx{$key}->{$this_gex_gene}->{eff} . "," . $fields[1] :   $fields[1];
		}
	}
	close $fh_gx;
	return \%gx;
}

sub load_genex_hash_ext{
	my ($mm, $st) = @_;
	my %gx = ();
	#open(GX, $mm) or warn "$mm $!";
	my $fh_gx = IO::File->new($mm , q{<}) or warn "$mm  $!\n";
	while(my $line = <$fh_gx>){
		chomp($line);
		my @fields = split("\t", $line);
		my $key = $fields[3] . ":" . $fields[4] . ":" . $fields[5] . ":" . $fields[6] . ":" . $fields[7];
		foreach my $each_gex (split("," , $fields[2] )){
			my $this_gex_gene = (split(":", $each_gex))[0];
			$gx{$key}->{$this_gex_gene}->{rfq} = (exists $gx{$key}->{$this_gex_gene}->{rfq}) ?  $gx{$key}->{$this_gex_gene}->{rfq} . "," .   $each_gex : $each_gex; 
			$gx{$key}->{$this_gex_gene}->{eff} = (exists $gx{$key}->{$this_gex_gene}->{eff}) ?  $gx{$key}->{$this_gex_gene}->{eff} . "," . $fields[1] :   $fields[1];
			$gx{$key . ":" . $st}->{$this_gex_gene}->{$fields[1]} = (exists $gx{$key . ":" . $st}->{$this_gex_gene}->{$fields[1]}) ? $gx{$key . ":" . $st}->{$this_gex_gene}->{$fields[1]} . "," .  $each_gex : $each_gex ;
		}
	}
	close $fh_gx;
	return \%gx;
}


sub load_hashBedDistance{
	my ($mm) = @_;
	my %dbval = ();
	
	#open(X, $mm) or warn "File(load_hash): $mm $!";
	my $fh_x = IO::File->new($mm , q{<}) or warn "$mm  $!\n";
	warn "********In load hash bed $mm\n";

#chr1    10002   10003   A       C|NA|NA|chr1    11868   11869   NR_148357:LOC102725121:ExonStart_1      100     +|1866
#chr1    10014   10015   A       C|NA|NA|chr1    11868   11869   NR_148357:LOC102725121:ExonStart_1      100     +|1854
#chr1    10018   10019   T       A|NA|NA|chr1    11868   11869   NR_148357:LOC102725121:ExonStart_1      100     +|1850

#chr1    84491   84492   C       G|chr1  70007   70008   NM_001005484:OR4F5:ExonEnd_1    100     +|-14484|chr1   134772  134773  NR_039983:LOC729737:ExonStart_1 100     -|50281
#chr1    91535   91536   G       T|chr1  70007   70008   NM_001005484:OR4F5:ExonEnd_1    100     +|-21528|chr1   134772  134773  NR_039983:LOC729737:ExonStart_1 100     -|43237
#chr1    91580   91581   G       A|chr1  70007   70008   NM_001005484:OR4F5:ExonEnd_1    100     +|-21573|chr1   134772  134773  NR_039983:LOC729737:ExonStart_1 100     -|43192
#chr1    98682   98683   G       A|chr1  70007   70008   NM_001005484:OR4F5:ExonEnd_1    100     +|-28675|chr1   134772  134773  NR_039983:LOC729737:ExonStart_1 100     -|36090

	
	while(my $line = <$fh_x>){
		chomp($line);
		#print "Bed -- $line\n";
		my @each_sec = split('\|', $line);
		my $var_info = shift @each_sec;		
		my @fieldsA = split("\t" , $var_info);
		my $key = $fieldsA[0] . ":" . ($fieldsA[1] + 1) . ":" . $fieldsA[2] . ":" . $fieldsA[3] . ":" . $fieldsA[4];
		
		#chr1    14907   14907
		while(@each_sec){
			my $this_sec = shift @each_sec;
			if($this_sec eq "NA"){
				$this_sec = shift @each_sec;
				next;
			}
			my @fieldsB = split("\t" , $this_sec);
			my $geneSec = $fieldsB[3];
			if($geneSec =~ m/^NR_/){
				$this_sec = shift @each_sec;
				next;
			}
			my $this_gn = (split(":" , $geneSec))[1];
			my $this_dist = shift @each_sec;

	#		$bed_dist{$key}{$this_gn} = (exists $bed_dist{$key}{$this_gn}) ?  $bed_dist{$key}{$this_gn} . "," .  $this_dist : $this_dist ;
			#print "$key , $this_gn -- $this_dist \n";
			$dbval{$key}{$this_gn} =  (exists $dbval{$key}{$this_gn} ) ? $dbval{$key}{$this_gn} . "," . $this_dist : $this_dist ;
		}
	}
	close $fh_x;
	return \%dbval;
	
}
sub load_hashBedDistance0{
	my ($mm) = @_;
	my %dbval = ();
	#open(X, $mm) or warn "File(load_hash): $mm $!";
	my $fh_x = IO::File->new($mm , q{<}) or warn "$mm  $!\n";
	warn "obs********In load hash bed\n";

#15      38973673        38973674        G       T|15    38988798        38988799        NM_207444:C15orf53:ExonStart_1  100     +|15125
#15      71491881        71491882        A       G|15    71481397        71481398        NR_120351:THSD4-AS1:ExonEnd_3   100     -|-10484
#15      71491884        71491888        TGTG    TAGT|15 71481397        71481398        NR_120351:THSD4-AS1:ExonEnd_3   100     -|-10487
	
	while(my $line = <$fh_x>){
		chomp($line);
		my @each_sec = split('\|', $line);
		#warn $each_sec[0] , "\n";
		my @fieldsA  = split("\t" , $each_sec[0]);
		my @fieldsB  = split("\t" , $each_sec[1]);
		my $this_spliceDist = $each_sec[2];
		my $this_strand = $fieldsB[5];
		my $this_spliceDist2go = ($this_spliceDist eq "-") ? (-1 * $this_spliceDist ) : $this_spliceDist;
		
		my $key = $fieldsA[0] . ":" . ($fieldsA[1] + 1) . ":" . $fieldsA[2] . ":" . $fieldsA[3] . ":" . $fieldsA[4];
		#warn "key = $key\n";
		$dbval{$key} = $this_spliceDist2go;
	}
	close $fh_x;
	return \%dbval;
}

#added Jan 10,2020
#1       chr1:69270:69270:A:G 
#2       OR4F5 
#3       79501
#4       OR4F5:NM_001005484:exon1:c.A180G:p.S60S 
#5       exonic
#6       exonic
#7       synonymous SNV
#8       synonymous SNV
#9       0.804
#10      100
#11      179



sub load_geneExtra{
	my ($mm) = @_;
	my %gva = ();
	#open(GV, $mm) or warn "$mm $!";
	my $fh_gv = IO::File->new($mm , q{<}) or warn "$mm  $!\n";
	while(my $line = <$fh_gv>){
		chomp($line);
		my @fields = split("\t", $line);
		my $key = $fields[0];   
		my $gn = $fields[1]; 
#		warn "Gene = $gn\n";
		$gva{$key}->{$gn}->{'gene'} = $fields[1];
		$gva{$key}->{$gn}->{'entrez_id'} = $fields[2];
		$gva{$key}->{$gn}->{'refseq-id'} = $fields[3];
		$gva{$key}->{$gn}->{'typeseq'} = $fields[4];
		$gva{$key}->{$gn}->{'typeseq_pr'} = $fields[5];
		$gva{$key}->{$gn}->{'effect'} = $fields[6];
		$gva{$key}->{$gn}->{'effect_pr'} = $fields[7];
		$gva{$key}->{$gn}->{'cds_prop' } = $fields[8];
		$gva{$key}->{$gn}->{'tr_prop'} = $fields[9];
		$gva{$key}->{$gn}->{'splice_dist'} = $fields[10];
		$gva{$key}->{$gn}->{'aa_flag'} = $fields[11];
		$gva{$key}->{$gn}->{'gene_desc'} = $fields[12];
		$gva{$key}->{$gn}->{'gene_type'} = $fields[13];
		$gva{$key}->{$gn}->{'omim_id'} = $fields[14];
		$gva{$key}->{$gn}->{'omim_phenotype'} = $fields[15];
		$gva{$key}->{$gn}->{'effect_appris'} = $fields[16];
		$gva{$key}->{$gn}->{'typeseq_select'} = $fields[17];
		$gva{$key}->{$gn}->{'effect_select'} = $fields[18];
		
	}
	close $fh_gv;
	return \%gva;
}

sub getGene {

#   to parse thsese entries: PLEKHN1:NM_001160184:exon12:c.C1253T:p.S418L,PLEKHN1:NM_032129:exon13:c.C1358T:p.S453L,
	my ($r)            = @_;
	my %result         = ();
	my %genes          = ();
	my %nm_ids         = ();
	my @am_changes     = ();
	my @trs            = split( ",", $r );
	my $this_gene      = "";
	my $this_aa_change = "";
	foreach (@trs) {

		if ( ( split( ":", $_ ) )[0] ne $this_gene ) {

			#	push(@genes, (split(":",$_))[0]);
			$genes{ ( split( ":", $_ ) )[0] } = 1;
			$nm_ids{ ( split( ":", $_ ) )[1] . ":" . ( split( ":", $_ ) )[0] } =
			  1;
		}
		$this_gene = ( split( ":", $_ ) )[0];
		my $aa = ( split( ":", $_ ) )[4];
		if ( $aa =~ m/p?(.)\d+(.)/ ) {
			if ( $this_aa_change ne "$1$2" ) {
				push( @am_changes, $aa );
			}
			$this_aa_change = "$1$2";
		}
	}
	my $res     = join( ",", keys %genes );
	my $res2    = join( ",", @am_changes );
	my $res_trs = join( ",", keys %nm_ids );
	my $flag    = 0;
	if ( ( scalar @am_changes ) > 1 ) {
		$flag = 1;
	}
	$result{"genes"}        = $res;
	$result{"aa_flag"}      = $flag;
	$result{"isoform_gene"} = $res_trs;
	return %result;
}




# output for help and errors
sub usage {
	my $error = shift;
	print $wfh "Critical: $error\n";
	warn "Error: $error\n\n" if ($error);
	exit(-1);
}

sub usage2 {
	my $error = shift;

	warn <<__EOT;



  ----
__EOT
	my $print_dt = localtime;
	print $wfh "$error\n";
	warn "**** $print_dt ::  Final merge Warning: $error ****\n\n" if ($error);

	# exit(-1);
}

sub usage0 {
	my $error = shift;

	print <<__EOT;
		
__EOT

	warn "Error: $error\n\n" if ($error);
	exit(-1);
}

sub extract_locusid {

	my ($l) = @_;
	my $lc_id = 0;
	my $info = ( split( "\t", $l ) )[12];
	if ( $info =~ m/.+LO=(\d+)/ ) {
		$lc_id = $1;
	}

	return $lc_id;
}

sub check_type {
	my ( $this_data, $expected ) = @_;
	my $sub_type = $expected;
	if ( looks_like_number($this_data) ) {
		$sub_type = 1;
	}
	elsif ( $this_data eq "NA" ) {
		$sub_type = 1;
	}
	else {
		$sub_type = 0;
	}

	my $ret_value = ( $sub_type == $expected ) ? 1 : 0;
	return $ret_value;
}

sub extract_gvarGenes {
	my ($rfinfo) = @_;

	#print "sub : $rfinfo\n";
	my @out_gene = ();
	my @spl = split( /\(.+?\)/, $rfinfo );
	foreach my $eg (@spl) {
		if ( $eg =~ m/^,/ ) {
			$eg =~ s/,//g;
		}
		push( @out_gene, $eg );
	}

	return join( ",", @out_gene );
}

sub compute_exon_min_dist {
	my ( $v_chr, $v_st, $v_end, $v_refseqid ) = @_;
	my $leftd    = undef;
	my $rightd   = undef;
	my @leftdis  = ();
	my @rightdis = ();
	foreach my $rf ( split( ",", $v_refseqid ) ) {
		my $rf_trid = ( split( ":", $rf ) )[1];
		my $rf_ex   = ( split( ":", $rf ) )[2];
		$rf_ex =~ s/exon//;
		my $rf_exons     = $refflats{ $rf_trid . ":" . $v_chr };
		my $rf_ex_strand = ( split( "\t", $rf_exons ) )[3];
		my $rf_ex_sts    = ( split( "\t", $rf_exons ) )[9];
		my $rf_ex_ends   = ( split( "\t", $rf_exons ) )[10];
		my @rf_ex_stsA   = split( ",", $rf_ex_sts );
		my @rf_ex_endsA  = split( ",", $rf_ex_ends );

		#depends on strand
		my $this_rf_ex_st  = $rf_ex_stsA[ $rf_ex - 1 ];
		my $this_rf_ex_end = $rf_ex_endsA[ $rf_ex - 1 ];

		if ( $rf_ex_strand eq "-" ) {
			my @rev_rf_ex_stsA  = reverse @rf_ex_stsA;
			my @rev_rf_ex_endsA = reverse @rf_ex_endsA;
			$this_rf_ex_st  = $rev_rf_ex_stsA[ $rf_ex - 1 ];
			$this_rf_ex_end = $rev_rf_ex_endsA[ $rf_ex - 1 ];
		}

		my $left = abs( $v_st - $this_rf_ex_st );
		push( @leftdis, $left );
		my $right = abs( $this_rf_ex_end - $v_end );
		push( @rightdis, $right );
		$leftd  = min(@leftdis);
		$rightd = min(@rightdis);
	}
	return ( $leftd, $rightd );

}

###not done yet
sub compute_exon_min_dist_uc {
	my ( $v_chr, $v_st, $v_end, $v_refseqid ) = @_;
	my $leftd    = undef;
	my $rightd   = undef;
	my @leftdis  = ();
	my @rightdis = ();
	foreach my $rf ( split( ",", $v_refseqid ) ) {
		my $rf_trid = ( split( ":", $rf ) )[1];
		my $rf_ex   = ( split( ":", $rf ) )[2];
		$rf_ex =~ s/exon//;
		my $rf_exons     = $uc_refflats{ $rf_trid . ":" . $v_chr };
		my $rf_ex_strand = ( split( "\t", $rf_exons ) )[2];
		my $rf_ex_sts    = ( split( "\t", $rf_exons ) )[8];
		my $rf_ex_ends   = ( split( "\t", $rf_exons ) )[9];
		my @rf_ex_stsA   = split( ",", $rf_ex_sts );
		my @rf_ex_endsA  = split( ",", $rf_ex_ends );

		#depends on strand
		my $this_rf_ex_st  = $rf_ex_stsA[ $rf_ex - 1 ];
		my $this_rf_ex_end = $rf_ex_endsA[ $rf_ex - 1 ];

		if ( $rf_ex_strand eq "-" ) {
			my @rev_rf_ex_stsA  = reverse @rf_ex_stsA;
			my @rev_rf_ex_endsA = reverse @rf_ex_endsA;
			$this_rf_ex_st  = $rev_rf_ex_stsA[ $rf_ex - 1 ];
			$this_rf_ex_end = $rev_rf_ex_endsA[ $rf_ex - 1 ];
		}

		my $left = abs( $v_st - $this_rf_ex_st );
		push( @leftdis, $left );
		my $right = abs( $this_rf_ex_end - $v_end );
		push( @rightdis, $right );
		$leftd  = min(@leftdis);
		$rightd = min(@rightdis);
	}
	return ( $leftd, $rightd );

}




#added on Dec 2014, to sort refseq id
sub sort_refseq_id {

#rid = PDE4DIP:NM_001002811:exon8:c.T2047C:p.L683L,PDE4DIP:NM_001002812:exon12:c.T1558C:p.L520L,PDE4DIP:NM_014644:exon12:c.T1558C:p.L520L,PDE4DIP:NM_001198832:exon15:c.T1756C:p.L586L,PDE4DIP:NM_001198834:exon12:c.T1558C:p.L520L,
	my ($rid)    = @_;
	my %tmp_rids = ();
	my @out_rid  = ();
	foreach my $this_rid ( split( ",", $rid ) ) {
		my $rid_tr = ( split( ":", $this_rid ) )[1];
		$tmp_rids{$this_rid} = $cd_lens{$rid_tr};
	}

	foreach my $id ( sort { $tmp_rids{$b} <=> $tmp_rids{$a} } keys %tmp_rids ) {
		push( @out_rid, $id );
	}

	return join( ",", @out_rid );

}

sub sort_uc_refseq_id {

#rid = PDE4DIP:NM_001002811:exon8:c.T2047C:p.L683L,PDE4DIP:NM_001002812:exon12:c.T1558C:p.L520L,PDE4DIP:NM_014644:exon12:c.T1558C:p.L520L,PDE4DIP:NM_001198832:exon15:c.T1756C:p.L586L,PDE4DIP:NM_001198834:exon12:c.T1558C:p.L520L,
	my ($rid)    = @_;
	my %tmp_rids = ();
	my @out_rid  = ();
	foreach my $this_rid ( split( ",", $rid ) ) {
		my $rid_tr = ( split( ":", $this_rid ) )[1];
		$tmp_rids{$this_rid} = $uc_cd_lens{$rid_tr};
	}

	foreach my $id ( sort { $tmp_rids{$b} <=> $tmp_rids{$a} } keys %tmp_rids ) {
		push( @out_rid, $id );
	}

	return join( ",", @out_rid );

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

sub validate_gene {

	my ($this_gene_name) = @_;
	my $valid_gene = 0;
	if ( exists $tr_counts{$this_gene_name} ) {
		$valid_gene = 1;
	}
	return $valid_gene;
}

sub compute_cds_prop {
	my ($rid) = @_;

	#my %cd_lens = %{$cdl};
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


sub print_gene_stats{
	my ($gene_stat_file) = @_;
	#open(GNST, ">$gene_stat_file") or die "Error $gene_stat_file $!";
	my $fh_gnst = IO::File->new($gene_stat_file , q{>}) or die "$! $gene_stat_file\n";
	
	print $fh_gnst "gene\tTotal\t" , join "\t" , ( sort keys %eff_categories) , "\n";
	foreach my $gs (keys %entz_id){
		my @out_line = ();
		my $gs_total = (exists $gene_stats{$gs}{'total'}) ? $gene_stats{$gs}{'total'} : 0;
		if($gs_total != 0){
		push(@out_line , $gs_total);
		foreach (sort keys %eff_categories){
			my $gs_val = (exists $gene_stats{$gs}{$_}) ? $gene_stats{$gs}{$_} : 0;
			push(@out_line , $gs_val);
		}
		
		print $fh_gnst $gs , "\t" , join("\t", @out_line) , "\n";
		}
	}
	close $fh_gnst;
	return 1;
}

sub check_precedence{
	my ($cur_ovr, $ovrs) = @_;
	my  %prec = (
        'splicing' => 2,
        'exonic' => 1,
        'ncRNA_splicing' => 4,
        'ncRNA_exonic' => 3,
        'UTR3' => 5,
        'UTR5' => 5,
        'intronic' => 6,
        'ncRNA_intronic' => 7,
        'upstream' => 8,
        'downstream' => 8,
        'intergenic' => 9,
    );
	my $ret = 1;
	my $this_prec = (exists $prec{$cur_ovr} ) ? ($prec{$cur_ovr}) : 1100;  
	foreach my $ov (split(";" , $ovrs )){
		my $pr = exists ($prec{$ov}) ? ($prec{$ov}) :1100;
		if($pr < $this_prec){
			$ret = 0;
		}
	}
	warn "ret = $ret\n";
	return $ret;	
}

sub check_precedence_ext{
	my ($cur_ovr, $ovrs, $povr) = @_;
	my %primary_ovrls = ();
	my  %prec = (
        'splicing' => 1,
        'exonic' => 1,
        'ncRNA_splicing' => 2,
        'ncRNA_exonic' => 2,
        'UTR3' => 5,
        'UTR5' => 5,
        'intronic' => 6,
        'ncRNA_intronic' => 7,
        'upstream' => 8,
        'downstream' => 8,
        'intergenic' => 9,
    );
	my $ret = 1;
	my $this_prec = (exists $prec{$cur_ovr} ) ? ($prec{$cur_ovr}) : 1100;  
	foreach my $ov (split(";" , $ovrs )){
		my $pr = exists ($prec{$ov}) ? ($prec{$ov}) :1100;
		if($pr < $this_prec){
			$ret = 0;
		}
	}
	return $ret;	
}

sub get_priortized_typeseq{
	my ($ovrs) = @_;
#	warn "Ins sub: $ovrs\n";
	my  %prec = (
        'exonic' => 1,
        'splicing' => 1,
        'ncRNA_splicing' => 2,        
        'ncRNA_exonic' => 2,
        'UTR3' => 3,
        'UTR5' => 3,
        'intronic' => 4,
        'ncRNA_intronic' => 5,
        'upstream' => 6,
        'downstream' => 6,
        'intergenic' => 7,
    );
    
    my @ov = split(";" , $ovrs);
    my $first = pop @ov;
    my $pr_first = $prec{$first};
    my @out = ($first);
    foreach (@ov){
    	my $this_pr = $prec{$_};
    	if($this_pr == $pr_first ){
    		push(@out , $_);
    	}
    	if($this_pr < $pr_first ){
    		@out = ($_);
    		$pr_first = $this_pr;
    	}
    }
    return join(";" , @out);
    
}


sub get_priortized_effect{
	my ($ovrs) = @_;
#	warn "Ins sub: $ovrs\n";
	my  %prec = (
        'frameshift insertion' => 1,
        'frameshift deletion' => 2,
        'frameshift block substitution' => 3,        
        'stopgain' => 4,
        'stoploss' => 5,
        'nonframeshift insertion' => 6,
        'nonframeshift deletion' => 7,
        'nonframeshift block substitution' => 8,
        'nonsynonymous SNV' => 9,
        'synonymous SNV' => 10,
        'UNKNOWN' => 11,
        'unknown' => 11,
        'unknown*' => 11,
    );
    
    my @ov = split("," , $ovrs);
    my $first = pop @ov;
    my $pr_first = $prec{$first};
    my @out = ($first);
    foreach (@ov){
    	my $this_pr = $prec{$_};
    	if($this_pr == $pr_first ){
    		push(@out , $_);
    	}
    	if($this_pr < $pr_first ){
    		@out = ($_);
    		$pr_first = $this_pr;
    	}
    }
    return join(";" , @out);
    
}



sub get_phylop_avarage{
	my ($phy) = @_;
	my @phys = split(";" , $phy); 
	return sum(@phys)/scalar(@phys);
}

sub get_lowestAbs{
	my ($list) = @_;
	my @tmplist = split("," , $list);
	my $lowest = shift(@tmplist);
	my $lowest_abs = abs($lowest);
	
	while(@tmplist){
		my $next_list = shift (@tmplist);
		if(abs($next_list) < $lowest_abs ){
			$lowest = $next_list;
			$lowest_abs = abs($next_list);
		}
	}
	
	return $lowest; 
	 
}

#my $exb_distance0 = getExb_distance($key, $each_gene);
sub getExb_distanceLow{
	my ($sub_key, $sub_gene) = @_;
	my $sub_dist = "NA";
	my @key_fields = split(":" , $sub_key);
	my $sub_chr = $key_fields[0];
	$sub_chr =~ s/chr//;
	my $gene_key = $sub_chr . ":" . $sub_gene;
	#print "Gene-key = $gene_key\n";
	my @dists = ();
	if(exists $gene_exons{$gene_key}{'exstart'} && exists $gene_exons{$gene_key}{'exend'}){
		#my @exsts = split("," , $gene_exons{$gene_key}{'exstart'});
		#my @exeds = split("," , $gene_exons{$gene_key}{'exend'});
		foreach ( split("," , $gene_exons{$gene_key}{'exstart'}) ){
			push(@dists, abs($key_fields[1] - ($_ + 1) ) ); 
		}
		foreach ( split("," , $gene_exons{$gene_key}{'exend'}) ){
			push(@dists, abs($key_fields[1] - $_) ); 
		}
		$sub_dist = min(@dists);
	}
	
	return $sub_dist;
}

#trying to make faster
sub getExb_distanceLow2{
	my ($sub_key, $sub_gene) = @_;
	my $d_limit = 1000000;
	my $sub_dist = "NA";
	my @key_fields = split(":" , $sub_key);
	my $sub_chr = $key_fields[0];
	$sub_chr =~ s/chr//;
	my $gene_key = $sub_chr . ":" . $sub_gene;
	#print "Gene-key = $gene_key\n";
	my @dists = ();
	if(exists $gene_exons{$gene_key}{'exstart'} && exists $gene_exons{$gene_key}{'exend'}){
		#my @exsts = split("," , $gene_exons{$gene_key}{'exstart'});
		#my @exeds = split("," , $gene_exons{$gene_key}{'exend'});
		foreach ( split("," , $gene_exons{$gene_key}{'exstart'}) ){
			#push(@dists, abs($key_fields[1] - ($_ + 1) ) );
			my $tmp_dist = abs($key_fields[1] - ($_ + 1));
			if($tmp_dist < $d_limit){
				$sub_dist = $tmp_dist;
				$d_limit = $tmp_dist;
			} 
		}
		foreach ( split("," , $gene_exons{$gene_key}{'exend'}) ){
#			push(@dists, abs($key_fields[1] - $_) ); 
			my $tmp_dist = abs($key_fields[1] - $_);
				if($tmp_dist < $d_limit){
					$sub_dist = $tmp_dist;
					$d_limit = $tmp_dist;
				} 
			
		}
		#$sub_dist = min(@dists);
	}
	
	return $sub_dist;
}

sub getExb_distance{
	my ($sub_key, $sub_gene) = @_;
	my $sub_dist = "NA";
	my @key_fields = split(":" , $sub_key);
	my $sub_chr = $key_fields[0];
	$sub_chr =~ s/chr//;
	my $gene_key = $sub_chr . ":" . $sub_gene;
	#print "Gene-key = $gene_key\n";
	my @dists = ();
	if(exists $gene_exons{$gene_key}{'exstart'} && exists $gene_exons{$gene_key}{'exend'}){
		my @exsts = split("," , $gene_exons{$gene_key}{'exstart'});
		my @exeds = split("," , $gene_exons{$gene_key}{'exend'});
		foreach (@exsts){
			push(@dists, ($key_fields[1] - ($_ + 1) ) ); 
		}
		foreach (@exeds){
			push(@dists, ($key_fields[1] - $_) ); 
		}
		
		$sub_dist = join("," , @dists);
	}
	
	return $sub_dist;
}

###INFO=<ID=SNVHPOL,Number=1,Type=Integer,Description="SNV contextual homopolymer length">

sub add2header{
	my ($hdr, $nm,  $type, $desc) = @_;
	my @hd2 = @{$hdr};
	#splice @dwarfs, 3, 0, 'SnowWhite';
	my $info2go  = qq(##INFO=<ID=)  . $nm . qq(,Number=) . $type . qq(,Type=String,Description=") . $desc . qq(">);
	splice @hd2 , 5, 0, $info2go; 
	return \@hd2;
}

sub add2headerFull{
	my ($hdr, $nms) = @_;
	my @hd2 = @{$hdr};
	my @nm2 = @{$nms};
	
	my @newinfo = ();
	#splice @dwarfs, 3, 0, 'SnowWhite';
	foreach (@nm2){
		my $info2go  = qq(##INFO=<ID=)  . $_ . qq(,Number=) . "A" . qq(,Type=String,Description=") . $_ . qq(">);
		push(@newinfo , $info2go);
	}
	splice @hd2 , 5, 0, @newinfo; 
	return \@hd2;
}

sub get_ovcfkey{
	my ($inf) = @_;
	my $ret = "-1";
	my ($this_info) = @_;
	if( $inf =~ m/OVCF_KEY=(.+?);/ || $this_info =~ m/OVCF_KEY=(.+)/ ){
		$ret = $1;
	} 
	return $ret;	
	
}

sub format_forVcf{
	my ($fld) = @_;
	$fld =~ s/,/\\x2c/g;
	$fld =~ s/;/\\x3b/g;
	$fld =~ s/=/\\x3d/g;
	return $fld;
	
}

sub get_uniqValue{
	my ($v) = @_;
	my @vs = split('\|' , $v);
	my @unique_vs = uniq @vs;
	return join('|' , @unique_vs);
	
}

sub get_uniqValueComma{
	my ($v) = @_;
	my @vs = split(',' , $v);
	my @unique_vs = uniq @vs;
	return join(',' , @unique_vs);
	
}

#my ($Clinvar_acc , $Clinvar_pheno , $Clinvar_ref, $Clinvar_reviewStat , $Clinvar_sig) = parse_clinvar($cln_val);
#ALLELEID=556509;
#CLNDISDB=MedGen:C4015293,OMIM:616126,Orphanet:ORPHA319563;
#CLNDN=Immunodeficiency_38_with_basal_ganglia_calcification;
#CLNHGVS=NC_000001.11:g.1014255G>A;
#CLNREVSTAT=criteria_provided,_single_submitter;
#CLNSIG=Uncertain_significance;
#CLNVC=single_nucleotide_variant;CLNVCSO=SO:0001483;GENEINFO=ISG15:9636;MC=SO:0001583|missense_variant;ORIGIN=1
#AF_ESP=0.01761;AF_EXAC=0.00493;AF_TGP=0.01418;ALLELEID=249265;CLNDISDB=MedGen:CN169374;CLNDN=not_specified;CLNHGVS=NC_000001.11:g.1022188A>G;CLNREVSTAT=criteria_provided,_multiple_submitters,_no_conflicts;CLNSIG=Benign;CLNVC=single_nucleotide_variant;CLNVCSO=SO:0001483;GENEINFO=AGRN:375790;MC=SO:0001627|intron_variant;ORIGIN=1;RS=115704555;ClinSigSimple=0

#return ($Clinvar_acc , $Clinvar_pheno , $Clinvar_ref, $Clinvar_reviewStat , $Clinvar_sig)
sub parse_clinvar{
	my ($cln_v) = @_;
	my $out_id = ($cln_v =~ m/ALLELEID=(.+?);/) ? $1: "NA";
	my $out_db = ($cln_v =~ m/CLNDISDB=(.+?);/) ? $1: "NA";
	my $out_dis = ($cln_v =~ m/CLNDN=(.+?);/) ? $1: "NA";
	my $out_rev = ($cln_v =~ m/CLNREVSTAT=(.+?);/) ? $1: "NA";
	my $out_sig = ($cln_v =~ m/CLNSIG=(.+?);/) ? $1: "NA";
	my $out_sig_sim = ($cln_v =~ m/ClinSigSimple=(.+)/) ? $1: "NA";
	if($out_sig ne "NA" && $out_sig_sim eq "NA"){
		die "Error parsing clinvar annotations from : $cln_v \n";
	}
	return ($out_id , $out_dis , $out_db , $out_rev , $out_sig, $out_sig_sim); 
} 

sub extract_mssngkey{
	my ($this_info) = @_;
	my $this_ret = -1;
	#MSSNG_KEY=chr1-10145-10147-AC-A;OVCF_KEY=chr1:10146:AC:A;MULTI_
	if($this_info =~ m/$pat_mssng/ ){
		$this_ret = $1;
	}else{
		die "Critical error! Error parsing MSSNG_KEY from $this_info\n";
	}
	return $this_ret;
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




sub threshold{
	my ($s_sft, $s_poly, $s_ma, $s_mt, $s_phymam, $s_phyvert, $s_cadd) = @_;
	my $number = 0;
	$number += ($s_sft < 0.05) ? 1: 0 unless ($s_sft eq "NA");
	$number += ($s_poly >= 0.9) ? 1: 0 unless ($s_poly eq "NA");
	$number += ($s_ma >= 1.9) ? 1: 0 unless ($s_ma eq "NA");
	$number += ($s_phymam >= 2.3) ? 1: 0 unless ($s_phymam eq "NA");
	$number += ($s_phyvert >= 4) ? 1: 0 unless ($s_phyvert eq "NA");
	$number += ($s_cadd >= 15) ? 1: 0 unless ($s_cadd eq "NA");
	$number += ($s_mt >= 0.5) ? 1: 0 unless ($s_mt eq "NA");
	#$self->{'thr'} = $number;
	return $number;
}


sub high{
	my $s_self = shift;
	my $e = shift;
	my $ret = 0;
	if($e eq 'frameshift' || $e eq 'stop_gain' || $e eq 'splice_site' || $e eq 'lof'){
		$ret = 1;
	}elsif($e eq 'missense'){
		$ret = ($s_self->{'thr'}  >= 4) ? 1: 0;
	}elsif($e eq 'other'){
		$ret = 0;
	}elsif($e eq 'reg_dec_exon'){
		$ret = ($s_self->{'spx'}  < -10 || $s_self->{'dbscAda'}  > 0.6 || $s_self->{'dbscRf'} > 0.6) ? 1: 0; 
	}elsif($e eq 'reg_inc_exon'){
		$ret = ($s_self->{'spx'} > 10) ? 1: 0;
	}elsif($e eq 'non_coding'){
		$ret = ( ($s_self->{'cadd'} >= 17.5 || $s_self->{'phymam'} >= 2.5 || $s_self->{'phyvert'} >= 4.5) && $s_self->{'phastCons'} ne "NA" ) ? 1: 0;
	}else{
		$ret = 0;
	}
	return $ret;
}

sub medium{
	my $self = shift;
	my $e  = shift;
	my $ret = 0;
	if($e eq 'frameshift' || $e eq 'stop_gain' || $e eq 'splice_site' || $e eq 'lof'){
		$ret = 0;
	}elsif($e eq 'missense'){
		$ret = (!high($e) &&  $self->{'thr'} >= 2) ? 1: 0;
	}elsif($e eq 'other'){
		$ret =  (( ($self->{'phymam'} >= 2.3 || $self->{'phyvert'} >= 4 || $self->{'cadd'} >= 20) && ($self->{'dbsnpc'} eq "NA")) ||  (($self->{'phymam'} >= 1.5 || $self->{'phyvert'} >= 2.5 || $self->{'cadd'} >= 15 ) && ( $self->{'dbsnp'} eq "NA" &&  $self->{'dbsnpr'} eq "NA")) ) ? 1: 0;
	}elsif($e eq 'reg_dec_exon'){
		$ret = (!high($e) && $self->{'spx'} < -2.5 ) ? 1: 0; 
	}elsif($e eq 'reg_inc_exon'){
		$ret = (!high($e) && $self->{'spx'} > 2.5) ? 1: 0;
	}elsif($e eq 'non_coding'){
		$ret = (!high($e) &&  ($self->{'cadd'} >= 17.5 || $self->{'phymam'} >= 2.5 || $self->{'phyvert'} >= 4.5) && $self->{'phastCons'} ne "NA" ) ? 1: 0;
	}else{
		$ret = 0;
	}
	return $ret;
}

sub low{
	my $self = shift;
	my $e = shift;
	my $ret = 0;
	my $thr = $self->{'thr'};
	if($e eq 'utr'){
		$ret = 1;
	}elsif($e eq 'frameshift' || $e eq 'stop_gain' || $e eq 'splice_site' || $e eq 'lof'){
		$ret = 0;
	}elsif($e eq 'missense'){
		$ret = (!high($e) && !medium($e)) ? 1: 0;
	}elsif($e eq 'other'){
		$ret =  (!high($e) && !medium($e)) ? 1: 0;
	}elsif($e eq 'reg_dec_exon'){
		$ret = (!high($e) && !medium($e) && $self->{'spx'} < 0 ) ? 1: 0; 
	}elsif($e eq 'reg_inc_exon'){
		$ret = (!high($e) && !medium($e) && $self->{'spx'} > 0 ) ? 1: 0;
	}elsif($e eq 'non_coding'){
		$ret = (!high($e) && !medium($e)) ? 1: 0;
	}else{
		$ret = 0;
	}
	return $ret;
}


sub get_effects{
	my $self = shift;
	my @to_ret = ();
	my $coding = (exists $EXONIC{$self->{'typeseq'}}) ? 1: 0;
	my $s_frameshift = (exists $FRAMESHIFT{$self->{'eff'}}) ? 1: 0;
	my $s_stopgain = (exists $stop_gain{$self->{'eff'}}) ? 1: 0;
	my $s_splicesite = ($self->{'typeseq'} eq "splicing") ? 1: 0;
	my $s_other = (exists $other{$self->{'eff'}}) ? 1: 0;
	my $s_missense = ($self->{'eff'} eq "nonsynonymous SNV") ? 1: 0;
	
	if(exists $FRAMESHIFT{$self->{'eff'}}){
		push(@to_ret , "frameshift");
	}
	if(exists $stop_gain{$self->{'eff'}}){
		push(@to_ret , "stop_gain");
	}
	if($s_splicesite){
		push(@to_ret , "splice_site");
	}
	if($coding && ($s_frameshift || $s_stopgain || $s_splicesite)){
		push(@to_ret , "lof");
	}
	
	if($s_missense){
		push(@to_ret , "missense");
	}
 	
 	if($s_other){
		push(@to_ret , "other");
 	}
 	
 	if($self->{'spx'} < 0 || $self->{'dbscAda'} > 0.6 || $self->{'dbscRf'} > 0.6){
 		push(@to_ret , "reg_dec_exon");
 	}
 	
 	if($self->{'spx'} > 0 ){
 		push(@to_ret , "reg_inc_exon");
 	}
 	
 	if($self->{'typeseq'} =~ m/^UTR/){
 		push(@to_ret , "utr");
 	}
 	
 	if(exists $NON_CODING{$self->{'typeseq'}}){
 		push(@to_ret , "non_coding");
 	}
 	
 	if(scalar(@to_ret) == 0){
 		push(@to_ret, "Unknown");
 	} 
 	 
 	return \@to_ret;
}


sub  bitwise_impact_of_effect{
	my $e = shift;
	my $im = shift;
	if(high($im, $e)){
		return $IMPACTS{'high'};
	}elsif(medium($im, $e)){
		return $IMPACTS{'medium'};
	}elsif(low($im, $e)){
		return $IMPACTS{'low'};
	}else{
		return '';
	}
}

sub decode_bit{
	my ($bt) = @_;
	my $ln = ($bt !=0) ? floor(log($bt) / log(2)) : -1;
	my $md = $ln % 3;
	my $str2ret = "NA";
	if(exists $EFFECTS_str{$ln - $md} && exists $IMPACTS_str{$md}){
		$str2ret = $EFFECTS_str{$ln - $md} . "-" . $IMPACTS_str{$md};
	}
	return $str2ret;
}


sub print_progress{
	my ($total, $current, $unit) = @_;
	my $no_of_dots = ($current / $unit);
	my $no_of_spaces = ($total / $unit) - $no_of_dots;
	my $per_fin = ($current / $total) * 100;
	#my $ret_string = "." x $no_of_dots . " (" . floor($per_fin) . " \%)" ;
	my $this_time = localtime;
	my $ret_string = "[" . $this_time . "]" . "." x $no_of_dots .  " " x $no_of_spaces  .    " (" . floor($per_fin) . " \%)" ;
	return $ret_string;
	
}

sub fecth_values{
	my $hs_name = shift;
	my $ky = shift;
	my $ret = shift;
	if ( exists $retrieved_hashes{$hs_name}->{$ky}  ) {
	  $ret = $retrieved_hashes{$hs_name}->{$ky};  	
	}
	return $ret;
}


