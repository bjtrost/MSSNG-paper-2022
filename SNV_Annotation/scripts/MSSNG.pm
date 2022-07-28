#!/bin/perl
use strict;
package MSSNG;
use File::Basename;
use Data::Dumper;
use Parallel::ForkManager;

our $AUTHOR =   '$Author: Thomas Nalpathamkalam <thomas.nalpathamkalam@sickkids.ca> $';

## Module for MSSNG file operations; convert to vcf format, extract information etc.

sub new{
	my $class = shift;
	my $filename = shift;
	my $self = {};
	bless ($self, $class);
	$self->{"name"} = $filename;
	return $self; 
}


sub convert2vcfOld{
	my $self = shift;
	my $outvcf = shift;
	#(chr, p_start, p_end, ref, alt, id) = line.chomp.split("\t")
	open(OUT, ">$outvcf") or die "$!";
	print OUT "##fileformat=VCFv4\.1\n";
	print OUT "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
	my $input = $self->{"name"};
	open(INP, $input) or die "$!";
	<INP>;  # assume header
	while(my $line = <INP>){
		chomp($line);
		my @fields = split("\t" , $line);
		print OUT $fields[0] , "\t" , ($fields[1] + 1) , "\t" ,  $fields[5] , "\t" , $fields[3] , "\t" , $fields[4] , "\t100\tPASS\tinfo=1\n";
	}
	close INP;
	close OUT;
}



#11      55595502        55595503        T       TTTTAA  AU4121302       0,1     50,6    -12.6,-0,-1664  PASS    56      99      11-55595502-55595503-T-TTTTAA   exonic  OR5L2:NM_001004739:exon1:c.809_810insTTTAA:p.
#5       52377527        52377528        T       G       7-0024-005      1,0     5,37    -24,0,-24       VQLOW   42      24      5-52377527-52377528-T-G splicing                        ITGA2   3673    integrin, alp
#1       24390611        24390613        AC      A       AU4121302       0,1     24,27   -78.1,-0,-74.4  PASS    51      99      1-24390611-24390613-AC-A        exonic  MYOM3:NM_152372:exon30:c.3571delG:p.V1191fs
#22      18072985        18072986        T       TA
#X       100749038       100749044       TGAGGC  T  
#created on Jan 24,2019
#original key is added

sub convert2vcf{
	my $self = shift;
	my $outvcf = shift;
	my $smp = shift;
	#(chr, p_start, p_end, ref, alt, id) = line.chomp.split("\t")
	open(OUT, ">$outvcf") or die "$!";
	print OUT "##fileformat=VCFv4\.2\n";
	print OUT qq(##INFO=<ID=MSSNG_KEY,Number=.,Type=String,Description="original 0-based key from mssng">) , "\n";
	print OUT qq(##FORMAT=<ID=GT,Number=1,Type=String,Description="GenotypeDummy">) , "\n";
	print OUT "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$smp\n";
	
	my $input = $self->{"name"};
	#if($input =~ /\.gz$/){
	if($input =~ /\.gz$/){
		warn "Compressed Input..\n";
		open(INP, "gunzip -c $input |") or die "gunzip $input: $!";
	}else{
		open(INP, $input) or die "$!";
	}	
	
	<INP>;  # assume header
	while(my $line = <INP>){
		chomp($line);
		my @fields = split("\t" , $line);
		my $mkey = $fields[0] . '-' . $fields[1] . '-' . $fields[2] . '-' . $fields[3] . '-' . $fields[4];
		my @infos = ();
		push(@infos, "MSSNG_KEY=" . $mkey);
		print OUT $fields[0] , "\t" , ($fields[1] + 1) , "\t" ,  $fields[5] , "\t" , $fields[3] , "\t" , $fields[4] , "\t100\tPASS\t" , join(";" , @infos) , "\tGT\t0\/1\n";  
	}
	close INP;
	close OUT;
}



sub convertFreq2vcf{
	my $self = shift;
	my $outvcf = shift;
	#(chr, p_start, p_end, ref, alt, id) = line.chomp.split("\t")
	open(OUT, ">$outvcf") or die "$!";
	print OUT "##fileformat=VCFv4\.1\n";
	print OUT "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
	my $input = $self->{"name"};
	open(INP, $input) or die "$!";
	my $header = <INP>;  # assume header
	chomp($header);
	$header =~ s/>//;
	my @hds = split("\t" , $header);
	while(my $line = <INP>){
		chomp($line);
		my @fields = split("\t" , $line);
		my $chrom = $fields[1];
		my $pos = $fields[2] + 1;
		my $id = $fields[0];
		my $ref = $fields[4];
		my $alt = $fields[17];
		my @infos = ();
		my $hd_idx = 0;
		foreach my $inf (@fields){
			push(@infos , $hds[$hd_idx++] . "=" . $inf );
		}
		
		print OUT $chrom , "\t" , $pos , "\t" , $id , "\t" , $ref , "\t" , $alt , "\t100\tPASS\t" , join(";" , @infos) , "\n";
	}
	close INP;
	close OUT;
}

sub annovarAnnotateMssng_freq{
	my $self <- shift;
	my $in_text = shift;
	my $out_file = shift;
	my $out_dir = shift;
	my $db_dir = shift;
	my $core = shift;
	
	my $dbsnp_r = $out_dir . "/" . basename($in_text) . "_sub.txt";
	my $dbsnp_w = $out_dir . "/" . basename($in_text) . "_win.txt";
	my $chr_info = $out_dir . "/" . basename($in_text) . "_chrInfo.txt";
	
	generate_sub_wind_inputs($in_text, $chr_info ,$out_dir);
	
	#parallel annotation
	#/hpf/tools/centos6/perl/5.20.1/bin/perl /hpf/largeprojects/tcagstor/users/thomas/tools/annovar/Feb_2016//annotate_variation.pl --buildver hg19  --outfile test.txt -filter -dbtype mssng /hpf/largeprojects/tcagstor/users/thomas/CG_annotation/data/package/archives/MSSNG/test/89931_S1_hFixed.vcf.gz_fixed.tsv.annovar.sub.txt /hpf/largeprojects/tcagstor/users/thomas/CG_annotation/data/package/archives/MSSNG/process/;
	
	my @cmds = ();
	my %jobs = ();
	my %dropped = ();
	#89931_S1_hFixed.vcf.gz_fixed.tsv.annovar.sub.txt_mssng.11.hg19_mssng_dropped
	open(CHR, $chr_info) or die "$!";
	while(my $line = <CHR>){
		chomp($line);
		my ($ch, $fl ) = split("\t" , $line);
		my $ch0 = $ch;
		$ch =~ s/chr//;
		$ch = ($ch eq "M") ? "MT" : $ch;
		
		my $anoutfile = $out_dir . "/" . basename($in_text) . "_mssng." . $ch;
		
		#hg19_mssng_7.f.txt
		my $dbname =  "mssng_" . $ch . ".f";
		#3-0380.3-0380-000_hfixed.vcf_fixed_1based.txt_mssng.3.hg19_mssng_3.f_dropped
		$dropped{"annovar:" . $ch} =  $out_dir . "/" . basename($in_text) . "_mssng." . $ch .  ".hg19_" . $dbname . "_dropped";
		my $an_cmd = "/hpf/tools/centos6/perl/5.20.1/bin/perl /hpf/largeprojects/tcagstor/users/thomas/tools/annovar/Feb_2016//annotate_variation.pl --buildver hg19  --outfile $anoutfile -filter -dbtype $dbname $fl $db_dir";
		push(@cmds, $an_cmd);
		$jobs{$an_cmd} = "annovar:" . $ch;
	}
	close CHR;

  my $pm = new Parallel::ForkManager($core); 
  
  my %status = ();
  
   $pm->run_on_finish(
    sub { my ($pid, $exit_code, $ident, $signal, $core_d) = @_;
      my $prestr = ($exit_code == 0) ? "***" : "";
      my $fin_time = localtime;
      print  " $prestr $ident *** finished with PID $pid and exit code: $exit_code , signal : $signal , core_dump : $core_d  \n";
      $status{$ident} = $exit_code;
      if($exit_code != 0 || $core_d || $signal){
      	warn "$prestr $ident *** finished with PID $pid and exit code: $exit_code  (SIGNAL = $signal & core_dump = $core_d) \n";
      }
    }
  );

  $pm->run_on_start(
    sub { my ($pid,$ident)=@_;
    	my $start_time = localtime;
      print  "\[$start_time\] ** $ident ** started, pid: $pid \n";
    }
  );

  #added Nov 05, 2014

#print Q1 "echo \"Program extra-gen finished with exit code \$\? at: `date`\" >> $log_file\n"; 
  my $prg = 1;
  foreach my $c (@cmds) {
  	my $idt = $jobs{$c};
  	
    $pm->start($idt) and next; # do the fork
    exec($c) or die "exec: $!";
    $pm->finish; # do the exit in the child process
    print "DONE : " , $c , "\n"; 
  }
  $pm->wait_all_children;
  
  #merge back
  foreach (keys %status){
  	if($status{$_} !=0){
  		die "$_ is finished with exit status: " , $status{$_} , "\n";
  	}
  }
  
  
  open(AUT, ">$out_file") or die "$!";
   
  
  foreach my $achr (sort keys %dropped){
  	 my $this_dropped = $dropped{$achr};
  	 open(DRP, $this_dropped) or die "$this_dropped : $!";
  	 while(my $line = <DRP>){
  	 	chomp($line);
  	 	print AUT $line , "\n";
  	 }
  	 warn "copied $this_dropped\n";
 # 	 unlink($this_dropped);
  }
  close DRP;
  warn "Annotation Finished , outfile = $out_file\n";
	
}


sub generate_sub_wind_inputs{
	my ($in, $out_cho, $tmp_dir) = @_;
	#print "$in, $out_sub, $out_win,$out_cho, $tmp_dir\n";
	my %fh;
	my %fhw;
	my %chr_info = ();
	open(IN, $in) or die "$!\n";
	open(CHO, ">$out_cho") or die "$!\n";

	while(my $line = <IN>){
		chomp($line);
		my @fields = split("\t", $line);
		my $ref = $fields[3];
		my $obs = $fields[4];
		my $id_key = $fields[6];
		my $vcf_chr = $fields[0];
		my $vcf_chr0 = $vcf_chr ;
		my $pos_start = $fields[1];
		my $pos_end = $fields[2];
		$vcf_chr =~ s/chr//;
#		my $wno_line = make_region($line,2500);
#		print WNO $wno_line , "\n";
		my $chr_filename = $tmp_dir . "/chr_" .  $vcf_chr . ".txt";
		$chr_info{$vcf_chr0} = $chr_filename ;
		#$chr_info{$vcf_chr0} = $chr_filename;
        if (not $fh{$chr_filename}) {
           open $fh{$chr_filename}, '>', $chr_filename or die "Unable to open '$chr_filename' for output: $!";
           warn "$chr_filename created\n";
        }
        print { $fh{$chr_filename} } $vcf_chr0, "\t" , $pos_start , "\t" , $pos_end , "\t" , $ref , "\t" , $obs , "\t" , $id_key , "\n"  ;		

	}

	foreach (keys %chr_info){
		print CHO $_ , "\t" , $chr_info{$_} , "\n" unless ($_ =~ /#/ || $_ eq "MT" );
	}

	close IN;
	close CHO;	
}




1;
