#!/usr/bin/perl
use strict;

our $AUTHOR =   '$Author: Thomas Nalpathamkalam <thomas.nalpathamkalam@sickkids.ca> $';

## Module to handle pipeline jobs. This will make various pipeline jobs and process parallely using Parallel::ForkManager 


package AnnotationEngine;

use JSON::XS qw(encode_json decode_json);
use threads;
use Data::Dumper;
use Data::Dumper qw(Dumper);
use Parallel::ForkManager;
use Time::Piece;
use Time::Seconds;
use FindBin;
use lib "$FindBin::Bin";
use PipelineExtra;
use IO::File;
use Term::ANSIColor;
use Term::ANSIColor qw( colored );
use File::Basename;


sub new{
	my $class = shift;
	my $mod_ref = shift;
	my $mod_name = shift;
	my $rep_file = shift;
	my $self = {};
	bless ($self, $class);
	$self->{"modref"} = $mod_ref;
	$self->{"modname"} = $mod_name;
	$self->{"outfile"} = $rep_file;
	
	#warn "rep = $rep_file\n";
	return $self; 
}

sub set_modules{
	my $self = shift;
	my $mod_ref = shift;
	$self->{"modref"} = $mod_ref;
	return $self; 
}


sub draw_workflow{
	my $self = shift;
	my $this_ref = $self->{"modref"};
	my $this_name = $self->{"modname"};
	my $this_file = $self->{"outfile"};
	open(OUT, ">$this_file") or die "$!";
	
	foreach my $kv (@{$this_name}){
		print OUT "\n\n";
		printf OUT ('+' . '-'x40 . '{%1s}' . '-'x40 . '+' , $kv );
		print OUT "\n";
		printf OUT ("|%5s = %40s|\n", "command" ,   $this_ref->{$kv}{command} );
		my $len_kv = length($kv);
		my $len_cmd = 2;
		printf OUT ("|%5s = %" . $len_cmd ."s|\n", "Start_time" ,   $this_ref->{$kv}{start_time} );
		printf OUT ("|%5s = %" . $len_cmd . "s|\n", "End_time" ,   $this_ref->{$kv}{end_time} );
		printf OUT ("|%5s = %" . $len_cmd . "s|\n", "exit_status" ,   $this_ref->{$kv}{exit_status} );
		printf OUT ("|%5s = %" . $len_cmd . "s|\n", "out_state" ,   $this_ref->{$kv}{out_state} );
		printf OUT ('+' . '-'x40 . "%" . ($len_kv+2) .  "s" . '-'x40 . '+'   ,  '*'x($len_kv+2) );
    	#print $kv, " : ", $this_ref->{$kv}{command}, "\n";
	}
	close OUT; 
}

sub draw_workflow_json{
	my $self = shift;
	my $this_ref = $self->{"modref"};
	my $this_name = $self->{"modname"};
	my $this_file = $self->{"outfile"};
	open(OUT, ">$this_file") or die "$!";
	print  JSON->new->pretty->encode($this_ref) , "\n";
	
	close OUT; 
}

sub draw_workflow_json2{
	my $self = shift;
	my $this_ref = $self->{"modref"};
	my $this_name = $self->{"modname"};
	my $this_file = $self->{"outfile"};
	open(OUT, ">$this_file") or die "$!";
	#print  JSON->new->pretty->encode($this_ref) , "\n";
	my @outs = ();
	foreach my $kv (@{$this_name}){
        if(exists $this_ref->{$kv}){
        	eval{
        	#push(@outs ,  "\"" .  $kv . "\" :" .   JSON->new->pretty->encode($this_ref->{$kv}) );
        	push(@outs ,  "\"" .  $kv . "\" :" .   encode_json($this_ref->{$kv}) );
        	}; warn $@ if $@;
        }
	}
	
	print OUT "{\n", join("," , @outs) , "}\n";
	close OUT; 
}


sub draw_workflow_json3{
	my $self = shift;
	my $this_ref = $self->{"modref"};
	my $this_name = $self->{"modname"};
	my $this_file = $self->{"outfile"};
	open(OUT, ">$this_file") or die "$!";
	my @out = ();
	my @outA = ();
	print OUT "{\n";
	foreach my $kv (@{$this_name}){
		my @out = ();
		push (@out , "\"" . $kv . "\":\n{");
		my $cmd = $this_ref->{$kv}{command};
		$cmd =~ s/\"/\\"/g;
		push(@out, "\"Command\" : \"" . $cmd  . "\"," );
		push(@out, "\"Start_time\" : \"" . $this_ref->{$kv}->{start_time}  . "\"," );
		push(@out, "\"End_time\" : \"" . $this_ref->{$kv}->{end_time}  . "\"," );
		push(@out, "\"Exit_status\" : \"" . $this_ref->{$kv}->{exit_status}  . "\"," );
		push(@out, "\"Out_state\" : \"" . $this_ref->{$kv}->{out_state}  . "\"}" );
		push (@outA, join ("\n" , @out )  . "\n");
 	}
	print OUT join("," , @outA ) , "\n";
	
	print OUT "}\n";
	
	close OUT; 
}


sub fork_manager{
	my $self = shift;
	my $this_ref = $self->{"modref"};
	my $this_name = $self->{"modname"};
	my $this_file = $self->{"outfile"};
	my @jobs2run = ();
	
	foreach (@{$this_name}){
		if(!exists $this_ref->{$_}{exit_status}){
			push(@jobs2run , $_);
		}
	}
	
	my $cores = shift;
	warn "Forking with $cores (v1) \n";	
	
	my $pm0 = new Parallel::ForkManager($cores);

    $pm0->run_on_finish(
    	sub { my ($pid, $exit_code, $ident, $signal, $core_d) = @_;
      		my $prestr = ($exit_code == 0) ? "***" : "";
      		my $fin_time = localtime;
      		$this_ref->{$ident}{end_time} = $fin_time;
      		$this_ref->{$ident}{exit_status} = $exit_code;  
      		print  "\[$fin_time\] $prestr $ident *** finished with PID $pid and exit code: $exit_code , signal : $signal , core_dump : $core_d  \n";
      		if($exit_code != 0 || $core_d || $signal){
      			warn "\[$fin_time\] $prestr $ident *** finished with PID $pid and exit code: $exit_code  (SIGNAL = $signal & core_dump = $core_d) \n";
      		}
    		}
  	);
  	
  $pm0->run_on_start(
    sub { my ($pid,$ident)=@_;
    	my $start_time = localtime;
      print  "\[$start_time\] ** $ident ** started, pid: $pid \n";
      $this_ref->{$ident}{start_time} = $start_time;
      $this_ref->{$ident}{end_time} = "NA";
      $this_ref->{$ident}{pid} = $pid;
      $this_ref->{$ident}{exit_status} = "NA";
    }
  );
  
  $pm0->run_on_wait(
    sub {
      job_status($this_ref , $this_name , $this_file);
    },
    60
  );

   my $prg = 1;
 
  foreach my $c (@jobs2run) {
  	my $idt = $c;
  	#$process_status{"conv2annovar-" . $cov2ann_jobid}{command} 
  	my $idt_cmd = $this_ref->{$idt}{command};
  	warn "$idt  $idt_cmd\n";
    $pm0->start($idt) and next; # do the fork
    exec($idt_cmd) ;
    $pm0->finish; # do the exit in the child process
    print "DONE : " , $c , "\n"; 
  }
  $pm0->wait_all_children;

  job_status($this_ref , $this_name , $this_file);
  
	
}

#Gene
sub fork_manager2{
	my $self = shift;
	my $this_ref = $self->{"modref"};
	my $this_name = $self->{"modname"};
	my $this_file = $self->{"outfile"};
	my @jobs2run = ();
	my @jobs2run2 = ();
	my %finished = ();
	my %towatch = ();
	my @watched = ();
	my $num_towatch = 0;
	my %priority_job = ();
	$priority_job{'dbsnp'} = 1;
	
	foreach (@{$this_name}){
		if($_ eq 'Welderly' || $_ eq 'cg1KB'){
			$priority_job{$_} = 1;
		}
		if(!exists $this_ref->{$_}{exit_status}){
			if(exists $this_ref->{$_}{level} && $this_ref->{$_}{level} == 2){
				push(@jobs2run2 , $_)
			}else{
				push(@jobs2run , $_);
			}
			if(exists $this_ref->{$_}{watch}){
				print "Adding to watchlist -> , $_\n";
				$num_towatch++;
			}
			
		}
	}
	
	my @orderd_jobs = keys %priority_job;
	foreach (@jobs2run){
		if(not exists $priority_job{$_}){
			push(@orderd_jobs, $_);
		}
	}
	
	print "Job2 = " , join(" -- " , @jobs2run2) , "\n";
	print "Number of jobs in watchlist = " , $num_towatch , "\n";
	my $cores = shift;
	warn "Forking with $cores (v2) \n";	
	
	my $pm0 = new Parallel::ForkManager($cores);

    $pm0->run_on_finish(
    	sub { my ($pid, $exit_code, $ident, $signal, $core_d) = @_;
      		my $prestr = ($exit_code == 0) ? "***" : "";
      		my $fin_time = localtime;
      		$this_ref->{$ident}{end_time} = $fin_time;
      		$this_ref->{$ident}{exit_status} = $exit_code;  
      		print  "\[$fin_time\] $prestr $ident *** finished with PID $pid and exit code: $exit_code , signal : $signal , core_dump : $core_d  \n";
      		if($exit_code != 0 || $core_d || $signal){
      			warn "\[$fin_time\] $prestr $ident *** finished with PID $pid and exit code: $exit_code  (SIGNAL = $signal & core_dump = $core_d) \n";
      		}else{
      			$finished{$ident} = 1;
      		}
    		}
  	);
  	
  $pm0->run_on_start(
    sub { my ($pid,$ident)=@_;
      my $start_time = localtime;
      my $this_jb_watch = (exists $this_ref->{$ident}{watch}) ? $this_ref->{$ident}{watch} : 0;
      #if($ident eq "Gene" || $ident eq "Gene-Rlx" || $ident eq "gene_info" || $ident eq "Gene-Appris" || $ident eq "Gene-Select" ){
      if($this_jb_watch == 1){	
      	$towatch{$ident} = $pid;
      }
      print  "\[$start_time\] ** $ident ** started, pid: $pid \n";
      $this_ref->{$ident}{start_time} = $start_time;
      $this_ref->{$ident}{end_time} = "NA";
      $this_ref->{$ident}{pid} = $pid;
      $this_ref->{$ident}{exit_status} = "NA";
    }
  );
  
  $pm0->run_on_wait(
    sub {
      job_status($this_ref , $this_name , $this_file);
    },
    60
  );

   my $prg = 1;
 
 
  my $l2_added = 0;
  
  while ( @orderd_jobs ) {
  	
  	$pm0->reap_finished_children;
  	print "l2_added = $l2_added\n";
  	if($l2_added == 0){
  		my $ct_fin = 0;
  		foreach (keys %towatch){
  			print "this_watch = " , $_  , "\n";
 			if(exists $finished{$_}){
 				print "this_watch done = " , $_  , "\n";
 				$ct_fin++;
 			}
 		}
 		print "ct_fin = $ct_fin\n";
 		if($ct_fin == $num_towatch ){
 			unshift(@orderd_jobs , @jobs2run2);
 			$l2_added = 1;
 		}
 		
  	}
  	
  	my $idt = shift @orderd_jobs;
  	$this_ref->{$idt}{exit_status} = "W";
  	warn "Last submitted = $idt\n";
  	warn "# of jobs waiting = " , scalar(@orderd_jobs) , "\n";
  	my $idt_cmd = $this_ref->{$idt}{command};
    my $this_pid = $pm0->start($idt) and next; # do the fork
  
    
  	exec($idt_cmd) ;
  	$pm0->finish; 
  	
  } #new loop 
  	
#  foreach my $c (@jobs2run) {
#  	
#  	
#  	my $idt = $c;
#  	#$process_status{"conv2annovar-" . $cov2ann_jobid}{command} 
#  	my $idt_cmd = $this_ref->{$idt}{command};
#  	warn "$idt  $idt_cmd\n";
#    my $this_pid = $pm0->start($idt) and next; # do the fork
#    if($idt eq "Gene"){
#    	$towatch{$this_pid} = $idt;
#    }
#    exec($idt_cmd) ;
#    $pm0->finish; # do the exit in the child process
#    print "DONE : " , $c , "\n"; 
#  }
$pm0->wait_all_children;
job_status($this_ref , $this_name , $this_file);
my $ang_ret = 0;
my @ret_jobs = ();
if($l2_added == 0){
	@ret_jobs = @jobs2run2;
	print "Job(s) still pending to start... \n";
	$ang_ret = 1;
}



#if($l2_added == 0){
#	print "Job(s) still pending to start... \n";
#	foreach my $jb2 (@jobs2run2){
#		my $idt2 = $jb2;
#		$this_ref->{$idt2}{exit_status} = "W";
#		my $idt2_cmd = $this_ref->{$idt2}{command};
#		$this_ref->{$idt2}{exit_status} = "R";
#		my $ex_out = exec($idt2_cmd);
#		if($ex_out == 0){
#			$this_ref->{$idt2}{exit_status} = $ex_out;
#			print "Job $idt2 finished\n";
#		}else{
#			warn "Job $idt2 failed! $ex_out \n";
#		}
#	}
#}
#job_status($this_ref , $this_name , $this_file);
return @ret_jobs;

}



sub job_status{
	my $job_monitor = shift;
	my $job_pids = shift;
	my 	$monitor_out = shift; 
   
   my $this_time = localtime;
   my $fh_jsto = IO::File->new($monitor_out , q{>}) or die "$! --- $monitor_out\n";
   printf $fh_jsto ("|%26s\t|%26s\t|%26s\t|%26s|\n", "module" , "pid" , "exit_status" , "Elapsed_Time(min)"  );
   print $fh_jsto "|";
   print  $fh_jsto '-' x 122;
   print $fh_jsto "|\n";
   
   my $total_jobs = 1;
   my $total_running = 0;
   my $total_fin_error = 0;
   my $total_fin_exit0 = 0;
   foreach my $jb (@{$job_pids}){
   $total_jobs++;
   my @spaninfo = ();
     
   my $this_job_start = (exists $job_monitor->{$jb}{start_time} ) ? $job_monitor->{$jb}{start_time} : "NA";
   my $this_job_end = (exists $job_monitor->{$jb}{end_time}) ? $job_monitor->{$jb}{end_time} : "NA";
   my $this_job_pid = (exists $job_monitor->{$jb}{pid}) ? $job_monitor->{$jb}{pid} : "NA";
   my $this_job_ext = (exists $job_monitor->{$jb}{exit_status}) ? $job_monitor->{$jb}{exit_status} : "NA";
   my $this_span = ($this_time-$this_job_start) ;
   if($this_job_pid ne "NA"  && $this_job_ext eq "NA" ){
   		$this_job_ext = "R";
   		$total_running++;
   }elsif($this_job_ext eq '0'){
   		 $total_fin_exit0++;
   }elsif($this_job_ext != 0){
   		$total_fin_error++;
   }
   if(!($this_job_ext eq "NA" || $this_job_ext eq "R")){
   		$this_span = ($this_job_end-$this_job_start) ;
   }
	my $span_min =  sprintf("%.2f", int($this_span)/60);
	  
   	#print GREEN, $fh_jsto $jb , "\t" , $this_job_pid , "\t"  ,  $this_job_ext  , "\t" , $span_min , "\n" ;
   	my $col_type = 'white';
	if($this_job_ext eq "R"){
		$col_type = 'CYAN';
	}elsif($this_job_ext eq "0"){
		$col_type = 'GREEN';
	}elsif($this_job_ext eq "NA"){
		$col_type = 'white';		
	}elsif($this_job_ext eq "W"){
		$col_type = 'yellow'  ;
	}else{
		$col_type = 'RED';
	}
	
	my $jb2go = $jb;
	if(length($jb) > 20){
		$jb2go = substr($jb, 0,20) . ".."; 
	}

   	printf $fh_jsto  ("|%26s\t|%26s\t|%26s\t|%26s\n", colored($jb2go, $col_type) , $this_job_pid  ,  $this_job_ext , $span_min  );
   }
   print $fh_jsto "|";
   print  $fh_jsto '-' x 122;
   print $fh_jsto "|\n";
    
   #print JSTO " \nX************************************X\n";
   $total_jobs--;
   print $fh_jsto "\nTotal Jobs : $total_jobs\n";
   print $fh_jsto "Running : $total_running\n";
   print $fh_jsto "Finished exit 0  : $total_fin_exit0\n";
   print $fh_jsto "Total Finished error : $total_fin_error\n";
   close $fh_jsto;
}


#"$perl_home $annovar_home\/annotate_variation.pl --buildver $build  -dbtype $annovar_gene_type --separate
# --outfile $gene_out_file $dbsnp_r  $annovar_gene_db > $gene_o 2> $gene_e  ";

#my $this_output = $geneuc_out_file . ".variant_function";
#$process_status{'Gene-UCSC'}{var_output} = $this_output;
#my $this_output = $geneuc_out_file . "\.exonic_variant_function";
#$process_status{'Gene-UCSC'}{exonic_output} = $this_output;

sub create_geneAnno_pipeline{
	my $self = shift;
	my $module_set = shift;
	my $this_ref = $self->{"modref"};	
	my $annovar_home = $this_ref->{'pipeline'}{annovar_home};
	my $build_ver = $this_ref->{'pipeline'}{build_ver};
	my $log_dir = $this_ref->{'pipeline'}{log_dir};
	my $tmp_dir = $this_ref->{'pipeline'}{tmp_dir};
	my $annovar_in = $this_ref->{'extra'}{dbsnp_r} ;
	my $sp_threshold = $this_ref->{'pipeline'}{splicing_threshold} ;
	my $neargene = $this_ref->{'pipeline'}{neargene} ;
	my $add_gene_opt = $this_ref->{'pipeline'}{add_gene_opt} ;
	
	
	my $base_annovar_in = basename($annovar_in); 
	my $annovar_gene_db = $this_ref->{'pipeline'}{gene_db} ;
	
	my %ret_pipes = ();
	foreach my $mk (keys %{$module_set}){
		my $mod_db = $module_set->{$mk}{db};
		my $annovar_gene_type = $module_set->{$mk}{gene_type};
		my $mod_dir = $tmp_dir . "/" . $mk;
		system("mkdir -p $mod_dir"); 
		my $db_dir = (exists $module_set->{$mk}{db_dir}) ? $module_set->{$mk}{db_dir} : $annovar_gene_db;
		my $out = $mod_dir . "/" . $base_annovar_in ;
		my $mod_stdout = $mod_dir . "/" . $base_annovar_in . ".stdout";
		my $mod_stderr = $mod_dir . "/" . $base_annovar_in . ".stderr";
		my $pipe_cmd = "perl $annovar_home\/annotate_variation.pl --buildver $build_ver -dbtype $annovar_gene_type --separate  --splicing_threshold $sp_threshold $add_gene_opt  --neargene $neargene --outfile $out $annovar_in $db_dir > $mod_stdout 2> $mod_stderr  " ;
		$ret_pipes{$mk}{command} = $pipe_cmd;
		my $dropped = $out . "." . $build_ver . "_" . $mod_db . "_dropped";
		$ret_pipes{$mk}{var_output} = $out . "\.variant_function";
		$ret_pipes{$mk}{exonic_output} = $out . "\.exonic_variant_function";
		$ret_pipes{$mk}{load_type} = "g"; 
		$ret_pipes{$mk}{error_file} = $mod_stderr; 
	}
	
	return \%ret_pipes;
}


#"$perl_home $annovar_home\/annotate_variation.pl 
#--buildver $build --otherinfo --outfile $dbscsnv_out_file 
# -filter -dbtype $filter_dbscsnv $dbsnp_r $annovar_ext_db > $dbscsnv_o 2> $dbscsnv_e ";


sub create_filter_pipeline{
	my $self = shift;
	my $module_set = shift;
	my $this_ref = $self->{"modref"};	
	my $annovar_home = $this_ref->{'pipeline'}{annovar_home};
	my $build_ver = $this_ref->{'pipeline'}{build_ver};
	my $log_dir = $this_ref->{'pipeline'}{log_dir};
	my $tmp_dir = $this_ref->{'pipeline'}{tmp_dir};
	my $annovar_in = $this_ref->{'extra'}{dbsnp_r} ;
	my $base_annovar_in = basename($annovar_in); 
	my $annovar_ext_db = $this_ref->{'pipeline'}{ext_db} ;
	
	my %ret_pipes = ();
	foreach my $mk (keys %{$module_set}){
		my $mod_db = $module_set->{$mk}{db};
		my $annovar_in2go  = (exists $module_set->{$mk}{'annovar-in'}) ? $module_set->{$mk}{'annovar-in'} : $annovar_in;
		my $mod_dir = $tmp_dir . "/" . $mk;
		my $out_prefix = (exists $module_set->{$mk}{db_prefix} ) ? $module_set->{$mk}{db_prefix} : $mod_db; 
		system("mkdir -p $mod_dir"); 
		my $db_dir = (exists $module_set->{$mk}{db_dir}) ? $module_set->{$mk}{db_dir} : $annovar_ext_db;
		my $this_database = $db_dir . "/" . $build_ver . "_" . $mod_db . ".txt";
		my $out = $mod_dir . "/" . $base_annovar_in ;
		my $mod_stdout = $mod_dir . "/" . $base_annovar_in . ".stdout";
		my $mod_stderr = $mod_dir . "/" . $base_annovar_in . ".stderr";
		my $pipe_cmd = "perl $annovar_home\/annotate_variation.pl --buildver $build_ver --otherinfo  --outfile $out -filter -dbtype $mod_db $annovar_in2go $db_dir > $mod_stdout 2> $mod_stderr" ;
		$ret_pipes{$mk}{command} = $pipe_cmd;
		my $dropped = $out . "." . $build_ver . "_" . $out_prefix . "_dropped";
		$ret_pipes{$mk}{output} = $dropped;
		$ret_pipes{$mk}{database} = $this_database;
		$ret_pipes{$mk}{load_type} = "f"; 
		$ret_pipes{$mk}{error_file} = $mod_stderr; 
	}
	return \%ret_pipes;
	
}




#cat $annovar_in_bed \| $bedops/bedmap  --skip-unmapped --exact --echo --echo-map  -  $bed_cadd > $cadd_out_file";
#$process_status{'bedops-sort'}{output} 
#$bedops_filter_module{'spidex'}{db} = $bed_spidex;
sub create_bedops_filter_pipeline{
	my $self = shift;
	my $module_set = shift;
	my $this_ref = $self->{"modref"};	
	my $log_dir = $this_ref->{'pipeline'}{log_dir};
	my $tmp_dir = $this_ref->{'pipeline'}{tmp_dir};
	my $bed_in = $this_ref->{'bedops-sort'}{output} ;
	my $base_bed_in = basename($bed_in); 
	#my $bed_db = ->{'pipeline'}{ext_db} ;
	#my $tool_home = $this_ref->{'pipeline'}{tool_home};
	my $script_dir = $this_ref->{'pipeline'}{tool_home} . "/scripts/";
	
	 
	my %ret_pipes = ();
	foreach my $mk (keys %{$module_set}){
		my $mod_db = $module_set->{$mk}{db};
		my $delim_param = (exists $module_set->{$mk}{delim}) ? "--delim " . $module_set->{$mk}{delim} : "";
		my $pp_command = qq( \|  perl -I $script_dir  -Mbedops -e 'Bedops::var_merge()' );
		if(exists $module_set->{$mk}{delim}){
			my $this_delim = $module_set->{$mk}{delim};
			$pp_command = qq( \|  perl -I $script_dir  -Mbedops -e 'Bedops::var_mergeMulti($this_delim)' ); #added 'var_mergeMulti' on Dec09,2021
		} 
		my $mod_dir = $tmp_dir . "/" . $mk;
		system("mkdir -p $mod_dir"); 
		#my $db_dir = (exists $module_set->{$mk}{db_dir}) ? $module_set->{$mk}{db_dir} : $annovar_ext_db;
		my $out = $mod_dir . "/" . $mk . "_out.bed";
		my $mod_stdout = $mod_dir . "/" . $mk . ".stdout";
		my $mod_stderr = $mod_dir . "/" . $mk . ".stderr";
		my $pipe_cmd = "cat $bed_in  \| bedmap $delim_param --skip-unmapped --exact --echo --echo-map - $mod_db $pp_command > $out 2> $mod_stderr" ;
		$ret_pipes{$mk}{command} = $pipe_cmd;
		my $dropped = $out;
		$ret_pipes{$mk}{output} = $dropped;
		$ret_pipes{$mk}{load_type} = "f";
		$ret_pipes{$mk}{error_file} = $mod_stderr; 		
	}
	
	return \%ret_pipes;
}

#bedmap --skip-unmapped --echo --echo-map $annovar_in_bed $dng_bed_file > $dng_region_out
sub create_bedops_region_pipeline{
	my $self = shift;
	my $module_set = shift;
	my $this_ref = $self->{"modref"};	
	my $log_dir = $this_ref->{'pipeline'}{log_dir};
	my $tmp_dir = $this_ref->{'pipeline'}{tmp_dir};
	my $bed_in = $this_ref->{'bedops-sort'}{output} ;
	my $base_bed_in = basename($bed_in); 
		#my $bed_db = ->{'pipeline'}{ext_db} ;
	
	my %ret_pipes = ();
	foreach my $mk (keys %{$module_set}){
		my $mod_db = $module_set->{$mk}{db};
		my $mod_dir = $tmp_dir . "/" . $mk;
		system("mkdir -p $mod_dir"); 
		#my $db_dir = (exists $module_set->{$mk}{db_dir}) ? $module_set->{$mk}{db_dir} : $annovar_ext_db;
		my $out = $mod_dir . "/" . $mk . "_out.bed";
		my $mod_stdout = $mod_dir . "/" . $mk . ".stdout";
		my $mod_stderr = $mod_dir . "/" . $mk . ".stderr";
		my $pipe_cmd = "cat $bed_in  \| bedmap --skip-unmapped  --echo --echo-map - $mod_db > $out 2> $mod_stderr" ;
		$ret_pipes{$mk}{command} = $pipe_cmd;
		my $dropped = $out;
		$ret_pipes{$mk}{output} = $dropped;
		$ret_pipes{$mk}{load_type} = "r";
		$ret_pipes{$mk}{error_file} = $mod_stderr; 		
	}
	
	return \%ret_pipes;
}

#$bedops/bedmap --chrom $chkey --skip-unmapped --echo --echo-map-score $annovar_in_bed 
# $phylop_100w_dir\/chr" . $cno . ".phyloP100way.bed_sorted  > $phylop_o 2> $phylop_e

sub create_bedops_region_chr_pipeline{
	my $self = shift;
	my $module_set = shift;
	my $this_ref = $self->{"modref"};	
	my $log_dir = $this_ref->{'pipeline'}{log_dir};
	my $tmp_dir = $this_ref->{'pipeline'}{tmp_dir};
	my $bed_in = $this_ref->{'bedops-sort'}{output};
	my $bed_mssng_in = $this_ref->{'mssng-bedsort'}{output};
	my $script_dir = $this_ref->{'pipeline'}{tool_home} . "/scripts/";
	
	my $base_bed_in = basename($bed_in); 
	
	my %ret_pipes = ();
	foreach my $mk (keys %{$module_set}){
		my $mod_db = $module_set->{$mk}{db};
		my $mod_chrkey = $module_set->{$mk}{chrkey} ;
		my $mod_dir = $tmp_dir . "/" . $mk;
		my $echo_map = (exists $module_set->{$mk}{echo_map} ) ? $module_set->{$mk}{echo_map} : " " ;
		my $bed_range = (exists $module_set->{$mk}{bed_range} ) ? $module_set->{$mk}{bed_range} : 0 ;
		my $this_pp =  (exists $module_set->{$mk}{post_process} ) ? $module_set->{$mk}{post_process} : "-1" ;
		#warn "*** this PP = $this_pp\n";
		my $pp_command = "";
		if($this_pp eq "phylop"){
			$pp_command = qq( \|  perl -I $script_dir  -Mbedops -e 'Bedops::phylop_merge()' );  
			 
		}elsif($this_pp eq "dbsnp"){
			$pp_command = qq( \|  perl -I $script_dir  -Mbedops -e 'Bedops::dbsnp_merge()' );
		}elsif($this_pp eq "cgint"){
			#warn "** cgInt\n";
			$bed_in = $bed_mssng_in;
			$pp_command = qq( \|  perl -I $script_dir  -Mbedops -e 'Bedops::cgfreq_merge()' );
		}
		
		
		system("mkdir -p $mod_dir"); 
		#my $db_dir = (exists $module_set->{$mk}{db_dir}) ? $module_set->{$mk}{db_dir} : $annovar_ext_db;
		my $out = $mod_dir . "/" . $mk . "_out.bed";
		my $mod_stdout = $mod_dir . "/" . $mk . ".stdout";
		my $mod_stderr = $mod_dir . "/" . $mk . ".stderr";
		my $pipe_cmd = "bedmap --chrom $mod_chrkey --skip-unmapped  --echo   $echo_map --range $bed_range  $bed_in $mod_db $pp_command > $out 2> $mod_stderr" ;
		$ret_pipes{$mk}{command} = $pipe_cmd;
		my $dropped = $out;
		$ret_pipes{$mk}{output} = $dropped;
		$ret_pipes{$mk}{load_type} = "r";
		$ret_pipes{$mk}{error_file} = $mod_stderr; 		
	}
	
	return \%ret_pipes;
}



#my $bedops_closf_cmd =  "closest-features --closest  --dist $annovar_in_bed $target_bed_sorted > $target_bedout 2> $target_bederr";
sub create_bedops_closest_pipeline{
	my $self = shift;
	my $module_set = shift;
	my $this_ref = $self->{"modref"};	
	my $log_dir = $this_ref->{'pipeline'}{log_dir};
	my $tmp_dir = $this_ref->{'pipeline'}{tmp_dir};
	my $bed_in = $this_ref->{'bedops-sort'}{output};
	my $script_dir = $this_ref->{'pipeline'}{tool_home} . "/scripts/";
	
	my $base_bed_in = basename($bed_in); 
	my %ret_pipes = ();

	foreach my $mk (keys %{$module_set}){
		my $mod_db = $module_set->{$mk}{db};
		my $mod_dir = $tmp_dir . "/" . $mk;
		
		system("mkdir -p $mod_dir"); 
		#my $db_dir = (exists $module_set->{$mk}{db_dir}) ? $module_set->{$mk}{db_dir} : $annovar_ext_db;
		my $out = $mod_dir . "/" . $mk . "_out.bed";
		my $mod_stdout = $mod_dir . "/" . $mk . ".stdout";
		my $mod_stderr = $mod_dir . "/" . $mk . ".stderr";
		my $pipe_cmd = "closest-features --closest  --dist  $bed_in $mod_db > $out 2> $mod_stderr";
		$ret_pipes{$mk}{command} = $pipe_cmd;
		my $dropped = $out;
		$ret_pipes{$mk}{output} = $dropped;
		$ret_pipes{$mk}{load_type} = "r";
		$ret_pipes{$mk}{error_file} = $mod_stderr; 		
	}
	
	return \%ret_pipes;
	
}

#$annovar_region_module{'phastcon'}{db} = $region_phst;
#my $phpl_command = "$perl_home $annovar_home\/annotate_variation.pl --buildver $build
# --outfile $phpl_out_file -regionanno -dbtype $region_phst $dbsnp_r $annovar_ext_db > $phpl_o  2> $phpl_e";

#my $pfam_command = "$perl_home $annovar_home\/annotate_variation.pl --buildver $build
# --outfile $pfam_out_file -bedfile $bed_pfam -regionanno -dbtype bed --colsWanted 4 $dbsnp_r $annovar_ext_db > $pfam_o 2> $pfam_e";

sub create_region_pipeline{
	my $self = shift;
	my $module_set = shift;
	my $this_ref = $self->{"modref"};	
	my $annovar_home = $this_ref->{'pipeline'}{annovar_home};
	my $build_ver = $this_ref->{'pipeline'}{build_ver};
	my $log_dir = $this_ref->{'pipeline'}{log_dir};
	my $tmp_dir = $this_ref->{'pipeline'}{tmp_dir};
	my $annovar_in = $this_ref->{'extra'}{dbsnp_r} ;
	my $base_annovar_in = basename($annovar_in); 
	my $annovar_ext_db = $this_ref->{'pipeline'}{ext_db} ;
	
	my %ret_pipes = ();
	foreach my $mk (keys %{$module_set}){
		#my $mod_db = $module_set->{$mk}{db};
		my $mod_dir = $tmp_dir . "/" . $mk;
		system("mkdir -p $mod_dir"); 
		my $db_dir = (exists $module_set->{$mk}{db_dir}) ? $module_set->{$mk}{db_dir} : $annovar_ext_db;
		my $db_type = (exists $module_set->{$mk}{dbtype}) ? $module_set->{$mk}{dbtype} : "";
		my $out_prefix = (exists $module_set->{$mk}{db_prefix} ) ? $module_set->{$mk}{db_prefix} : $db_type;
		my $db_bed = "";
		if($db_type eq "bed"){
			$db_bed  = (exists $module_set->{$mk}{bed}) ? $module_set->{$mk}{bed} : die "Critical! Missing bed file for $mk annotation!\n";
		}
		
		my $out = $mod_dir . "/" . $base_annovar_in ;
		my $mod_stdout = $mod_dir . "/" . $base_annovar_in . ".stdout";
		my $mod_stderr = $mod_dir . "/" . $base_annovar_in . ".stderr";
		my $dropped = $out . "." . $build_ver . "_" . $out_prefix;
		#$pfam_out_file . "."  . $build ."_bed";
		#$phpl_out_file . "." . $build . "_" .  $region_phst ;
		my $pipe_cmd = "perl $annovar_home\/annotate_variation.pl --buildver $build_ver  --outfile $out -regionanno -dbtype $db_type $annovar_in $db_dir > $mod_stdout 2> $mod_stderr" ;
		if($db_type eq "bed"){
			$pipe_cmd = "perl $annovar_home\/annotate_variation.pl --buildver $build_ver  --outfile $out -regionanno -dbtype $db_type -bedfile $db_bed  --colsWanted 4 $annovar_in $db_dir > $mod_stdout 2> $mod_stderr" ;
			$dropped = $out . "." . $build_ver . "_bed";
		}
		$ret_pipes{$mk}{command} = $pipe_cmd;
		$ret_pipes{$mk}{output} = $dropped;
		$ret_pipes{$mk}{load_type} = "r";
		$ret_pipes{$mk}{error_file} = $mod_stderr; 
	}
	
	return \%ret_pipes;
	
}

#beta
sub create_compile_pipeline{
	my $self = shift;
	my $module_set = shift;
	my $this_ref = $self->{"modref"};
	my $log_dir = $this_ref->{'pipeline'}{log_dir};
	my $tmp_dir = $this_ref->{'pipeline'}{tmp_dir};
	my $cmp_script = $this_ref->{'pipeline'}{cmp_script};
	my $meta_file = $this_ref->{'pipeline'}{meta_file};
	
	my %ret_pipes = ();
	foreach my $mk (keys %{$module_set}){
		my $run_chr = $module_set->{$mk}{chrom};
		my $out_file = $tmp_dir . "/tmp_out_" . $run_chr; 
		my $err_file = $tmp_dir . "/tmp_out_" . $run_chr . "_merge.log"; 
		my $pipe_cmd = "perl $cmp_script -m $meta_file -t $tmp_dir -c $run_chr > $out_file 2> $err_file";
		$ret_pipes{$mk}{command} = $pipe_cmd;
		$ret_pipes{$mk}{output} = $out_file;
		$ret_pipes{$mk}{error_file} = $err_file; 
		
	}
	return \%ret_pipes;	
	
}

#tabix -h VenterBlood500ng_singleton.vcf.gz.norm.vcf.gz chr1 | perl  /hpf/largeprojects/tcagstor/users/thomas/Annotation/Human/2020/dev/rev27.2/scripts/tabix_search_v3.pl 
#-o out2.txt -d /hpf/largeprojects/tcagstor/users/thomas/Annotation/Human/2020/data/package/archives/gnomAD/vcf/v3.0/hg38/gnomad.genomes.r3.0.sites.vcf.bgz 
# --chunks 5000
# --out-tags AC_female,AC_male,AF_afr,AF_ami,AF_amr,AF_asj,AF_eas,AF_female,AF_fin,AF_male,AF_oth,AF_raw,AF_sas,AN_female,AN_male,faf95_adj,faf95_afr,faf95_amr,faf95_eas,faf95_nfe,faf95_adj,nhomalt_female,nhomalt_male 
# > out_testGn3.tsv

sub create_custom_pipeline{
	my $self = shift;
	my $module_set = shift;
	my $this_ref = $self->{"modref"};	
	my $log_dir = $this_ref->{'pipeline'}{log_dir};
	my $tmp_dir = $this_ref->{'pipeline'}{tmp_dir};
	my $annovar_in = $this_ref->{'extra'}{dbsnp_r} ;
	my $base_annovar_in = basename($annovar_in);
	
	my %ret_pipes = ();
	foreach my $mk (keys %{$module_set}){
		my $mod_db = $module_set->{$mk}{db};
		my $md_cmd = $module_set->{$mk}{cmd};
		#my $md_out0 =  $module_set->{$mk}{output};
		my $this_level = $module_set->{$mk}{level};
		my $this_stat = $module_set->{$mk}{stats};
		
		my $mod_dir = $tmp_dir . "/" . $mk;
		system("mkdir -p $mod_dir"); 
		my $md_out = $mod_dir . "/gene_extra.tsv";
		my $mod_stdout = $mod_dir . "/" . $base_annovar_in . ".stdout";
		my $mod_stderr = $mod_dir . "/" . $base_annovar_in . ".stderr";
		my $pipe_cmd = "$md_cmd > $md_out 2> $mod_stderr" ;
		$ret_pipes{$mk}{command} = $pipe_cmd;
		$ret_pipes{$mk}{output} = $md_out;
		$ret_pipes{$mk}{load_type} = "c"; 
		$ret_pipes{$mk}{error_file} = $mod_stderr;
		$ret_pipes{$mk}{level} = $this_level;  
		$ret_pipes{$mk}{stats} = $this_stat;
	}
	return \%ret_pipes;
}

sub create_vcfAnno_pipeline{
	my $self = shift;
	my $module_set = shift;
	my $this_ref = $self->{"modref"};	
	my $log_dir = $this_ref->{'pipeline'}{log_dir};
	my $tmp_dir = $this_ref->{'pipeline'}{tmp_dir};
	#'extra'}{dbsnp_r
	my $vcf_in = $this_ref->{'extra'}{dbsnp_r} ;
	my $base_vcf_in = basename($vcf_in);
	print "**create_vcfAnno_pipeline = $base_vcf_in\n";
	#print Dumper($module_set) , "\n";
	
	my %ret_pipes = ();
	foreach my $mk (keys %{$module_set}){
		my $mod_db = $module_set->{$mk}{db};
		my $md_cmd = $module_set->{$mk}{cmd};
		#my $md_out0 =  $module_set->{$mk}{output};
		my $this_level = $module_set->{$mk}{level};
		my $this_stat = $module_set->{$mk}{stats};
		my $this_chr = (exists $module_set->{$mk}{chrom}) ? $module_set->{$mk}{chrom} : "";
		my $this_script = $module_set->{$mk}{script};
		my $this_inf_tags =  $module_set->{$mk}{info_tags};
		my $this_in_type = $module_set->{$mk}{type};
		
		my $mod_dir = $tmp_dir . "/" .  $mk;
		system("mkdir -p $mod_dir"); 
		my $md_out = $mod_dir . "/".  $base_vcf_in . "_out.tsv";
		my $mod_stdout = $mod_dir . "/" . $base_vcf_in . ".stdout";
		my $mod_stderr = $mod_dir . "/" . $base_vcf_in . ".stderr";
		my $pipe_cmd = qq(sort -k 1,1 -k 2,2n $this_chr | perl $this_script -o out.txt -d $mod_db --chunks 8000  --out-tags $this_inf_tags --type $this_in_type > $md_out 2> $mod_stderr) ;
		$ret_pipes{$mk}{command} = $pipe_cmd;
		$ret_pipes{$mk}{output} = $md_out;
		$ret_pipes{$mk}{load_type} = "gn"; 
		$ret_pipes{$mk}{error_file} = $mod_stderr;
		$ret_pipes{$mk}{level} = $this_level;  
		$ret_pipes{$mk}{stats} = $this_stat;
		$ret_pipes{$mk}{norm_vcf} =  $base_vcf_in;
		$ret_pipes{$mk}{script} =  $this_script;
		$ret_pipes{$mk}{info_tags} =  $this_inf_tags;
		
	}
	return \%ret_pipes;
}


1;

