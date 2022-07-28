#!/bin/perl
use strict;

our $AUTHOR =   '$Author: Thomas Nalpathamkalam <thomas.nalpathamkalam@sickkids.ca> $';

## Module for annovar operations. 


use IO::File;
use File::Basename;

package Annovar;


sub new{
	my $class = shift;
	my $gene_dir = shift;
	my $ext_dir = shift;
	my $anno_home = shift;
	my $build = shift;
	my $chr_prefix = 1;
	
	my $self = {};
	bless ($self, $class);
	$self->{"gene_dir"} = $gene_dir;
	$self->{"ext_dir"} = $ext_dir;
	$self->{"anno_home"} = $anno_home;
	$self->{"build"} = $build;
	$self->{"chr_prefix"} = $chr_prefix;
	
	return $self; 
}


sub set_chrPrefix{
	my $self = shift;
	my $prefix_1 = shift;
	$self->{"chr_prefix"} = $prefix_1;
	return $self;
}


sub filter_annotation{
	my $self = shift;
	my $ann_inp = shift;
	my $db_prefix = shift;
	my $out_prefix = shift;
	my $std_err = $out_prefix . ".stderr";
	my $std_out = $out_prefix . ".stdout";
	
	#my $dnsf_command = "$perl_home $annovar_home\/annotate_variation.pl --buildver hg19 --otherinfo --outfile $dnsf_out_file -filter -dbtype $filter_dbsnfp $dbsnp_r $annovar_ext_db > $dnsf_o 2> $dnsf_e ";
	my $anno_cmd = $self->{"anno_home"} . "/annotate_variation.pl --buildver " . $self->{"build"} . " --otherinfo --outfile " . $out_prefix . " -filter -dbtype " . $db_prefix .  " $ann_inp " . $self->{"ext_dir"} . " > " . $std_out . " 2> " . $std_err; 
	return $anno_cmd;
}

#my $phpl_command = "$perl_home $annovar_home\/annotate_variation.pl --buildver hg19 --outfile $phpl_out_file -regionanno -dbtype $region_phst $dbsnp_r $annovar_ext_db > $phpl_o  2> $phpl_e
sub region_annotation{
	my $self = shift;
	my $ann_inp = shift;
	my $db_prefix = shift;
	my $out_prefix = shift;
	my $std_err = $out_prefix . ".stderr";
	my $std_out = $out_prefix . ".stdout";
	my $anno_cmd = $self->{"anno_home"} . "/annotate_variation.pl --buildver " . $self->{"build"} . " --outfile " . $out_prefix . " -regionanno -dbtype " . $db_prefix .  " $ann_inp " . $self->{"ext_dir"} . " > " . $std_out . " 2> " . $std_err;
	return $anno_cmd;
}

#my $phyloppr_command = "$bedops/bedmap --chrom $chkey --skip-unmapped --echo --echo-map-score
# $annovar_in_bed  $phylop_100w_dir\/chr" . $cno . ".phyloP100way.bed_sorted  > $phylop_o 2> $phylop_e";
#beta
sub bedops_annotation{
	my $self = shift;
	my $in_bed = shift;
	my $db_bed = shift;
	my $out_bed = shift;
	my $out_err = shift;
	my $which_chr = shift;
}

#my $cadd_command = "cat $annovar_in_bed \| $bedops/bedmap  
#--skip-unmapped --exact --echo --echo-map  -  $bed_cadd > $cadd_out_file";
sub bedmap_filter_annotation{
	my $self = shift;
	my $annovar_in_bed = shift;
	my $bed_db = shift;
	my $bed_out = shift;
	my $bed_err = shift;
	#my $which_chr = shift;
	
	my $bed_f_cmd = "cat $annovar_in_bed \| bedmap --skip-unmapped --exact --echo --echo-map - $bed_db > $bed_out 2> $bed_err";
	return $bed_f_cmd;
}

sub bedmap_region_annotation{
	my $self = shift;
	my $annovar_in_bed = shift;
	my $bed_db = shift;
	my $bed_out = shift;
	my $bed_err = shift;
	#my $which_chr = shift;
	
	my $bed_f_cmd = "cat $annovar_in_bed \| bedmap --skip-unmapped  --echo --echo-map - $bed_db > $bed_out 2> $bed_err";
	return $bed_f_cmd;
}


1;