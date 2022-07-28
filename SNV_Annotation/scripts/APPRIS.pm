#!/usr/bin/perl
use strict;
use File::Basename;

our $AUTHOR =   '$Author: Thomas Nalpathamkalam <thomas.nalpathamkalam@sickkids.ca> $';

## To use APPRIS gene files ( ie, instead of using refseq genes ) 


package Appris;

sub new{
	my $class = shift;
	my $princ = shift;
	my $self = {};
	bless ($self, $class);
	$self->{"principal"} = $princ;
	return $self; 
}


sub make_gene_db{
	my $self = shift;
	my $gene_db = shift;
	my $ref_gen = shift;
	my $db_pre = shift;
	my $types = shift;
	my $annovar = shift;
	
	my $genome_dir_pre = "--seqdir";
	if(-f $ref_gen ){
        warn "$ref_gen is seqfile!";
        $genome_dir_pre = "--seqfile";
        
	}
	
	
	my $app_gene_out = $db_pre . ".txt";
	my $app_mrna_out = $db_pre . "Mrna.fa";

	my %ap_cat = (
    	'P1'    => 'PRINCIPAL:1',
    	'P2'   => 'PRINCIPAL:2',
    	'P3' => 'PRINCIPAL:3', 
    	'P4'  => 'PRINCIPAL:4', 
    	'P5'    => 'PRINCIPAL:5',
    	'ALT1' => 'ALTERNATIVE:1',
    	'ALT2' => 'ALTERNATIVE:2',
	);	
	
	my $apris_file = $self->{"principal"};
	
	my %in_types = ();
	foreach my $vl (split("," , $types)){
		$vl =~ s/^\s+//; #remove leading spaces
		$vl =~ s/\s+$//; #remove trailing spaces		
		
		if(exists $ap_cat{$vl}){
			$in_types{$ap_cat{$vl}} = $vl;
		}else{
			die "Error $vl cannot map to appris categories! allowed values:  " , join (", " , (sort keys %ap_cat)) ,  "\n";
		}
	}
	
	my %app_trs = ();
	open(APR, $apris_file) or die "$! $apris_file \n";
	while(my $line = <APR>){
		chomp($line);
		my @fields = split("\t" , $line);
		my $this_tr = $fields[2];
		my $this_apr = $fields[4];
		$this_apr =~ s/^\s+//; #remove leading spaces
		$this_apr =~ s/\s+$//; #remove trailing spaces		
		
		if(exists $in_types{$this_apr}){
			my $this_tr2go = (split('\.' , $this_tr))[0];
			$app_trs{$this_tr2go} = 1;
		}
	}
	close APR;
	warn "Loaded appris transcripts_ids\n";
	
	warn "Filtering gene_db $gene_db..\n";
	
	open(GDB, $gene_db) or die "$! $gene_db\n";
	open(GDO, ">$app_gene_out") or die "$! $app_gene_out\n";
	
	while(my $line = <GDB>){
		chomp($line);
		my @fields = split("\t" , $line);
		my $this_tr = $fields[1];
		my $this_tr2go = (split('\.' , $this_tr))[0];
		if(exists $app_trs{$this_tr2go}){
			print GDO $line , "\n";
		}
	}
	close GDB;
	
	warn "filtered gene db : $app_gene_out generated!\n";
	
	my $gene_mrna_cmd = "perl $annovar/retrieve_seq_from_fasta.pl $app_gene_out $genome_dir_pre $ref_gen -format refGene -outfile $app_mrna_out";
	warn "$gene_mrna_cmd\n";
	my $ext = system($gene_mrna_cmd);
	if($ext == 0){
		warn "filtered mrna generated sucessfully: $app_mrna_out\n";
	}else{
		die "Error generating mrna $app_mrna_out!\n";
	}
	
}

sub make_db_from_custom_def{
	my $self = shift;
	my $gene_db = shift;
	my $ref_gen = shift;
	my $db_pre = shift;
	my $annovar = shift;
	
	my $genome_dir_pre = "--seqdir";
	if(-f $ref_gen ){
        warn "$ref_gen is seqfile!";
        $genome_dir_pre = "--seqfile";
        
	}
	
	
	my $app_gene_out = $db_pre . ".txt";
	my $app_mrna_out = $db_pre . "Mrna.fa";

	open(GDB, $gene_db) or die "$! $gene_db\n";
	open(GDO, ">$app_gene_out") or die "$! $app_gene_out\n";
	while(my $line = <GDB>){
		chomp($line);
		print GDO $line , "\n";
	}
	close GDB;
	close GDO;
	
	my $gene_mrna_cmd = "perl $annovar/retrieve_seq_from_fasta.pl $app_gene_out $genome_dir_pre $ref_gen -format refGene -outfile $app_mrna_out";
	warn "$gene_mrna_cmd\n";
	my $ext = system($gene_mrna_cmd);
	if($ext == 0){
		warn "filtered mrna generated sucessfully: $app_mrna_out\n";
	}else{
		die "Error generating mrna $app_mrna_out!\n";
	}
	
	

	
}

1;
