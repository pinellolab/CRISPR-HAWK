package format_features_25;
use strict;
use warnings;

my ($dH_ref, $dS_ref, $dG_ref, $dHi, $dSi, $dGi) = rna_dna_param();
my %dG = %{$dG_ref};
	
sub get_seq_features{
	my $seq = shift;
	my $leftflank_numbases = shift;
	my $rightflank_numbases = shift;
	################################################################################################
	# Constants/Parameters
	################################################################################################
	# my $leftflank_numbases = 0;
	my $seed_numbases = 20;
	# my $rightflank_numbases = 3;
	my $pam_numbases = 3;
	my $scaffold_seq = "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTT";

	################################################################################################
	# All Feature Variables
	################################################################################################
	my $class_n = 0;
	my $pam_loc = $leftflank_numbases + $seed_numbases;
	my @feat_name_list;
	my @feat_vals;
	my @feat_type;
	
	my %seq_parts;
	$seq_parts{'all'} = $seq;
	$seq_parts{'lflank'} = substr($seq, 0,$leftflank_numbases);
	$seq_parts{'seed'} = substr($seq, $leftflank_numbases,$seed_numbases);
	$seq_parts{'pam'} = substr($seq, $leftflank_numbases+$seed_numbases,$pam_numbases);
	$seq_parts{'rflank'} = substr($seq, $leftflank_numbases+$seed_numbases+$pam_numbases, $rightflank_numbases);
	################################################################################################
	# Feature class 1: Position-specific base composition 
	#	All nucleotides in target, pam, flanking regions except GG in PAM
	################################################################################################
	$class_n ++;
	#For all positions in the input sequence
	for (my $indx = 0; $indx < length($seq); $indx++) {
		# skip the GG in Pam NGG #nvm, filter features in later step
		#next if ($indx == $pam_loc+1) or ($indx == $pam_loc+2);

		
		#   0 is the first nuc of the seed
		# -15 is the first nuc of the left flank
		#$info_header .= "$pos:A\t$pos:C\t$pos:G\t$pos:T\t";
		
		# record the feature (binary)
		my $base = substr($seq, $indx, 1);
		my $a_base = 0;
		my $c_base = 0;
		my $g_base = 0;
		my $ut_base = 0;

		if ($base =~ /^a$/i) {
			$a_base = 1;
		}
		elsif ($base =~ /^c$/i) {
				$c_base = 1;
		}
		elsif ($base =~ /^g$/i) {
				$g_base = 1;
		}
		elsif ($base =~ /^[ut]$/i || $base =~ /^t$/i) {
				$ut_base = 1;
		}
		
		#convert position relative to the seed
		my $pos = $indx-$leftflank_numbases;
		my $part = 'seed';
		$part = 'lflank' if $pos < 0;
		$part = 'pam' if $pos > $seed_numbases-1 and $pos < $seed_numbases+$pam_numbases;
		$part = 'rflank' if $pos > $seed_numbases+$pam_numbases-1;
		# record the feature and its name
		#Feature: <pos>:A
		push @feat_vals, $a_base; push @feat_name_list, "class$class_n" . '_' . "$pos:A_$part" ; push @feat_type, 'binary';
		#Feature: <pos>:T
		push @feat_vals, $ut_base; push @feat_name_list, "class$class_n" . '_' . "$pos:UT_$part"; push @feat_type, 'binary';
		#Feature: <pos>:C
		push @feat_vals, $c_base; push @feat_name_list, "class$class_n" . '_' . "$pos:C_$part"; push @feat_type, 'binary';
		#Feature: <pos>:G
		push @feat_vals, $g_base; push @feat_name_list, "class$class_n" . '_' . "$pos:G_$part"; push @feat_type, 'binary';
		
	}
	################################################################################################
	# Feature class 2: UUU (within last six bases); PolyX in the gRNA
	################################################################################################
	$class_n ++;
	foreach my $part ('seed'){
		my $seq_part = $seq_parts{$part};
		$seq_part =~ tr/ATCGU/atcgu/; $seq_part =~ tr/t/u/;
		
		my $polyA_reject_length = 5;
		my $polyC_reject_length = 5;
		my $polyG_reject_length = 4;
		my $polyU_reject_length = 4; # this was a prefilter implemented for the NBT dataset
		
		#Feature: polyX<_part>
		my $polyX_flag = 0;
		$polyX_flag = 1 if $seq_part =~ /a{$polyA_reject_length}/i || $seq_part =~ /c{$polyC_reject_length}/i || $seq_part =~ /g{$polyG_reject_length}/i || $seq_part =~ /[ut]{$polyU_reject_length}/i;
		push @feat_vals, $polyX_flag; push @feat_name_list, "class$class_n" . '_' .  "polyX" . "_$part"; push @feat_type, 'binary';
		
		#Feature: UUU<_part>
		my $UUU_flag = 0;
		$UUU_flag = 1 if substr($seq_part, -6, 6) =~ /uuu/i;
		push @feat_vals, $UUU_flag; push @feat_name_list, "class$class_n" . '_' . "UUU" . "_$part"; push @feat_type, 'binary';
	}
	
	################################################################################################
	# Feature class 3: GC Content
	################################################################################################
	$class_n ++;
	foreach my $part ('lflank','seed','rflank'){
		my $seq_part = $seq_parts{$part};
		next if length($seq_part) == 0;
		#Feature: gc_content<_part>
		my $gc_content = gcPercent($seq_part);
	
		# 1 if GC content of seed is high (>80%), 0 otherwise
		{ 
			my $gc_content_flag = 0;
			$gc_content_flag = 1 if $gc_content > 0.8;
			push @feat_vals, $gc_content_flag; push @feat_name_list, "class$class_n" . '_' . "GC_pc_gt_80" . "_$part"; push @feat_type, 'binary';
		}
		
		# 1 if GC content of seed is low (<30%), 0 otherwise

		{ 
			my $gc_content_flag = 0;
			$gc_content_flag = 1 if $gc_content < 0.3;
			push @feat_vals, $gc_content_flag; push @feat_name_list, "class$class_n" . '_' . "GC_pc_lt_30" . "_$part"; push @feat_type, 'binary';
		}
	}
	################################################################################################
	# Feature class 4: Duplex binding stability of the gRNA-seed
	################################################################################################
	$class_n ++;
	foreach my $part ('seed'){
		my $seq_part = $seq_parts{$part};
		$seq_part =~ tr/ATCGU/atcgu/; $seq_part =~ tr/t/u/;
		
		my $binding_flag = 0;
		for (my $start = 0; $start < length($seq_part) - 4; $start++) {
			next unless $start == 0 or $start == 7;
			
			#Feature: binding<n_part>
			my $dG_binding = dG_binding(substr($seq_part, $start, length($seq_part)-$start));
			push @feat_vals, $dG_binding; push @feat_name_list, "class$class_n" . '_' . "binding" . ($start + 1) . "_$part"; push @feat_type, 'numerical';
			
			#Feature: binding_flag<_part>
			$binding_flag = 1 if ($start == 0 and $dG_binding >= -18) or ($start == 7 and $dG_binding < -22) or ($start == 7 and $dG_binding > -9);
		}
		
		push @feat_vals, $binding_flag; push @feat_name_list, "class$class_n" . '_' . "binding_flag" . "_$part"; push @feat_type, 'binary';
	}
	
	################################################################################################
	# Feature class 5: Base accessibility (alignment)
	################################################################################################
	$class_n ++;
	foreach my $part ('seed'){
		my $seq_part = $seq_parts{$part};
		$seq_part =~ tr/ATCGU/atcgu/; $seq_part =~ tr/t/u/;
		
		my $dG_folding = foldingdG($seq_part);
		my $dG_folding_flag = 0;
		$dG_folding_flag = 1 if $dG_folding < -8;
		
		#Feature: dG_folding_flag<_part>
		push @feat_vals, $dG_folding_flag; push @feat_name_list, "class$class_n" . '_' . "dG_folding_flag" . "_$part"; push @feat_type, 'binary';

		#Feature: dG_folding<_part>
		push @feat_vals, $dG_folding; push @feat_name_list, "class$class_n" . '_' . "dG_folding" . "_$part"; push @feat_type, 'numerical';
		
		#alignment block **************************************************************************************************************************************************
		my $gRNA = $seq_part . $scaffold_seq;
		my ($dG_folding1, $alignment) = RNA_fold($gRNA);
		#print "$alignment\n"; exit;
		my $ext_stem = "(((((((((.((((....))))...)))))))";
		my $aligned_stem = substr($alignment, 18, length($ext_stem));
		
		#Feature: ali_flag<_part>
		my $align_flag = 0;
		$align_flag = 1 if $aligned_stem eq $ext_stem;
		push @feat_vals, $align_flag; push @feat_name_list, "class$class_n" . '_' . "ali_flag" . "_$part"; push @feat_type, 'binary';

		#Feature: ali<n_part>
		$alignment =~ tr/\.\(\)/011/;
		my @aligned = split "", $alignment;
		for (my $i = 0; $i <= $#aligned; $i++) {
			next unless $i == 14 or $i == 17 or $i == 18 or $i == 19 or $i == 20 or $i == 50 or $i == 51 or $i == 52;
			push @feat_vals, $aligned[$i]; push @feat_name_list, "class$class_n" . '_' . "ali" . ($i + 1) . "_$part"; push @feat_type, 'binary';
		}
		#end alignment block***********************************************************************************************************************************************
	}
	
	################################################################################################
	# Feature class 6: mono-, di- and tri-nucleotide composition
	################################################################################################
	$class_n ++;
	#monomer
	foreach my $part ('lflank','seed','rflank'){ 
		next if length($seq_parts{$part}) == 0;
		my @base_ary; #keep this temporary variables in this scope 
		#Feature: A_count<_part>
		@base_ary = $seq_parts{$part} =~ /a/gi; push @feat_vals,scalar(@base_ary); push @feat_name_list, "class$class_n" . '_' . "A_count_" . $part; push @feat_type, 'numerical';
		#Feature: UT_count<_part>
		@base_ary = $seq_parts{$part} =~ /[ut]/gi; push @feat_vals,scalar(@base_ary); push @feat_name_list, "class$class_n" . '_' . "UT_count_" . $part; push @feat_type, 'numerical';
		#Feature: C_count<_part>
		@base_ary = $seq_parts{$part} =~ /c/gi; push @feat_vals,scalar(@base_ary); push @feat_name_list, "class$class_n" . '_' . "C_count_" . $part; push @feat_type, 'numerical';
		#Feature: G_count<_part>
		@base_ary = $seq_parts{$part} =~ /g/gi; push @feat_vals,scalar(@base_ary); push @feat_name_list, "class$class_n" . '_' . "G_count_" . $part; push @feat_type, 'numerical';
	}
	
	#dimer
	foreach my $part ('lflank','seed','rflank'){ 
		next if length($seq_parts{$part}) < 2; #need at least 2 bases for dimer
		my %dimer_count;
		my $seq_part = $seq_parts{$part};
		$seq_part =~ tr/ATCGU/atcgu/; $seq_part =~ tr/u/t/;
		for (my $pos = 0; $pos < length($seq_part)-1; $pos++) {
			my $dimer = substr($seq_part, $pos, 2);
			$dimer_count{$dimer}++;
		}
		
		#Feature: <dimer><_part>
		foreach my $base1 ('a','t','c','g') {
			foreach my $base2 ('a','t','c','g') {
				my $dimer = $base1 . $base2;
				$dimer_count{$dimer} = 0 if !exists $dimer_count{$dimer};
				push @feat_vals, $dimer_count{$dimer}; push @feat_name_list, "class$class_n" . '_' . $dimer . "_$part"; push @feat_type, 'numerical';
			}
		}
	}
	
	#trimer
	foreach my $part ('lflank','seed','rflank'){ 
	next if length($seq_parts{$part}) < 3; #need at least 3 bases for trimer
		my %trimer_count;
		my $seq_part = $seq_parts{$part};
		$seq_part =~ tr/ATCGU/atcgu/; $seq_part =~ tr/u/t/;
		for (my $pos = 0; $pos < length($seq_part)-2; $pos++) {
			my $trimer = substr($seq_part, $pos, 3);
			$trimer_count{$trimer}++;
		}
		
		#Feature: <trimer><_part>
		foreach my $base1 ('a','t','c','g') {
			foreach my $base2 ('a','t','c','g') {
				foreach my $base3 ('a','t','c','g') {
					my $trimer = $base1 . $base2 . $base3;
					$trimer_count{$trimer} = 0 if !exists $trimer_count{$trimer};
					push @feat_vals, $trimer_count{$trimer}; push @feat_name_list, "class$class_n" . '_' . $trimer . "_$part"; push @feat_type, 'numerical';
				}
			}
		}
	}
	
	################################################################################################
	# Feature class 7: Purine/Pyrimidine content
	################################################################################################
	$class_n ++;
	my ($feat_vals_ref,$feat_name_list_ref,$feat_type_ref) = 
	feature_pupy_content(\%seq_parts, \@feat_vals,\@feat_name_list,\@feat_type,$class_n);
	@feat_vals = @{$feat_vals_ref};
	@feat_name_list = @{$feat_name_list_ref};
	@feat_type = @{$feat_type_ref};
	
	################################################################################################
	# Feature class 8: NGGN Pam dinucleotide
	################################################################################################
	$class_n ++;
	($feat_vals_ref,$feat_name_list_ref,$feat_type_ref) = 
	feature_nggn_pam(\%seq_parts, \@feat_vals,\@feat_name_list,\@feat_type,$class_n);
	@feat_vals = @{$feat_vals_ref};
	@feat_name_list = @{$feat_name_list_ref};
	@feat_type = @{$feat_type_ref};
	
	################################################################################################
	# Output
	################################################################################################
	return (\@feat_vals,\@feat_name_list,\@feat_type);
	
}


################################################################################################
# Feature class 7: Purine/Pyrimidine content
################################################################################################
sub feature_pupy_content{
	my ($seq_parts_ref, $feat_vals_ref,$feat_name_list_ref,$feat_type_ref,$class_n) = @_;
	my %seq_parts = %{$seq_parts_ref};
	my @feat_vals = @{$feat_vals_ref};
	my @feat_name_list = @{$feat_name_list_ref};
	my @feat_type = @{$feat_type_ref};
	   $class_n = $class_n;
	
	foreach my $part ('seed'){
		my $seq = $seq_parts{$part};
		$seq =~ tr/ATCGU/atcgu/; $seq =~ tr/u/t/;
		
		my @base_ary_a = $seq_parts{$part} =~ /a/gi; 
		my @base_ary_g = $seq_parts{$part} =~ /g/gi; 
		push @feat_vals,scalar(@base_ary_a) + scalar(@base_ary_g); push @feat_name_list, "class$class_n" . '_' . "purine_count_" . $part; push @feat_type, 'numerical';
		
		my @base_ary_ut = $seq_parts{$part} =~ /[ut]/gi; 
		my @base_ary_c = $seq_parts{$part} =~ /c/gi; 
		push @feat_vals,scalar(@base_ary_ut) + scalar(@base_ary_c); push @feat_name_list, "class$class_n" . '_' . "pyrimidine_count_" . $part; push @feat_type, 'numerical';
		
	}
	return (\@feat_vals,\@feat_name_list,\@feat_type);
}
################################################################################################
# Feature class 8: NGGN Pam dinucleotide
################################################################################################
sub feature_nggn_pam{
	my ($seq_parts_ref, $feat_vals_ref,$feat_name_list_ref,$feat_type_ref,$class_n) = @_;
	my %seq_parts = %{$seq_parts_ref};
	my @feat_vals = @{$feat_vals_ref};
	my @feat_name_list = @{$feat_name_list_ref};
	my @feat_type = @{$feat_type_ref};
	   $class_n = $class_n;
	
	return (\@feat_vals,\@feat_name_list,\@feat_type) if length($seq_parts{'rflank'}) == 0; #need at least 1 rflank base for this feature 
	
	my $nggn = substr($seq_parts{'pam'},0,1) . substr($seq_parts{'rflank'},0,1);
	$nggn =~ tr/ATCGU/atcgu/; $nggn =~ tr/u/t/;

	#Feature: <dimer>_NGGN_pam
	foreach my $base1 ('a','t','c','g') {
		foreach my $base2 ('a','t','c','g') {
			my $dimer = $base1 . $base2;
			if ($dimer eq $nggn){
				push @feat_vals, 1; push @feat_name_list, "class$class_n" . '_' . $dimer . "_NGGN_pam_rflank"; push @feat_type, 'binary';
			}else{
				push @feat_vals, 0; push @feat_name_list, "class$class_n" . '_' . $dimer . "_NGGN_pam_rflank"; push @feat_type, 'binary';
			}
		}
	}
	return (\@feat_vals,\@feat_name_list,\@feat_type);
}
################################################################################################
# Feature class 9: Windowed GC content
################################################################################################
sub feature_GC_win10_content{
	my ($seq_parts_ref, $feat_vals_ref,$feat_name_list_ref,$feat_type_ref,$class_n) = @_;
	my %seq_parts = %{$seq_parts_ref};
	my @feat_vals = @{$feat_vals_ref};
	my @feat_name_list = @{$feat_name_list_ref};
	my @feat_type = @{$feat_type_ref};
	$class_n = $class_n;
	my $window_size = 10;
	
	foreach my $part ('seed'){
		my $seq = $seq_parts{$part};
		$seq =~ tr/ATCGU/atcgu/; $seq =~ tr/u/t/;
		

		for (my $i = 0; $i+$window_size-1< length($seq);$i++){
			my $start = $i;
			my $end = $i+$window_size-1;
			my $gc_content = gcPercent(substr($seq_parts{$part},$i,$window_size));
			#push @feat_vals, $gc_content; push @feat_name_list, "class$class_n" . '_' . "GC_count_" . $start . '_' . $end; push @feat_type, 'numerical';
			
			for (my $j = 10; $j<100;$j+=10){
				if ($gc_content >= $j/100 ){
					push @feat_vals, 1; push @feat_name_list, "class$class_n" . '_' . $start . '_' . $end . "_GC_pc_gt_" . $j . "_$part"; push @feat_type, 'binary';
				}else{
					push @feat_vals, 0; push @feat_name_list, "class$class_n" . '_' . $start . '_' . $end . "_GC_pc_gt_" . $j . "_$part"; push @feat_type, 'binary';
				}
			}
		}
	}
	return (\@feat_vals,\@feat_name_list,\@feat_type);
}

################################################################################################
# Misc Subroutines
################################################################################################

################################################################################################
# Subroutine: get_svm_features
# formats the seq features in to a format ready for LIBSVM input
################################################################################################
sub get_svm_features{
	my $feat_val_ref = shift;
	my @feat_vals = @{$feat_val_ref};
	my $label = shift;
	
	my @output;
	push @output, $label if $label ne 'NaN';
	my $index = 0;
	for( my $i=0;$i <= $#feat_vals;$i++){
		push @output, ($i+1) . ':' .$feat_vals[$i];
	}
	return join ' ', @output;
}
################################################################################################
# Subroutine: dG_binding
# The overall deltaG of the target binding duplex. This is similar to GC content
################################################################################################
sub dG_binding {
	my $oligo = shift;
	$oligo =~ s/u/t/g;
	my $binding_dG = 0;
	for (my $indx = 0; $indx < length($oligo) - 1; $indx++) {
		$binding_dG += $dG{substr($oligo, $indx, 2)};
	}
	$binding_dG += $dGi;

	return $binding_dG;
}

################################################################################################
# Subroutine: foldingdG
# Calculate the secondary structure delta G for a single input sequence
# Input - a single RNA oligo sequence
# Output - the deltaG value
################################################################################################
sub foldingdG {
   my $sequence = shift;

   my $tempSeq = "tempSeq$$";
   my $tempOUT = "tempOUT$$";

   $sequence =~ s/[5|3|'|\-|\s+]//g;
   $sequence =~ tr/Tt/Uu/;

   open(OLIGO, ">$tempSeq") or die "can not open $tempSeq for writing: $!\n";
   print OLIGO $sequence;
   close OLIGO;

   my $dG;
   system("./RNAfold < $tempSeq > $tempOUT");
   open(RESULT, "$tempOUT") or die "Cannot open $tempOUT for reading $!\n";
   while(my $line = <RESULT>){
      # .((((((((((((((((((((((((((....)))))))))))))))))))))))))) (-17.80)
      #print $line, "\n";
      if($line =~ /([\-|\d][\.|\d]+)\)/){
         $dG = $1;
         last;
      }
   }
   close(RESULT);
   unlink $tempSeq;
   unlink $tempOUT;
   return $dG;
}
################################################################################################
# Subroutine: RNA_fold
# Calculate the secondary structure delta G for a single input sequence
# Input - a single RNA oligo sequence
# Output - the deltaG value
################################################################################################
sub RNA_fold {

   my $sequence = shift;

   my $tempSeq = "tempSeq$$";
   my $tempOUT = "tempOUT$$";

   $sequence =~ s/[5|3|'|\-|\s+]//g;
   $sequence =~ tr/Tt/Uu/;

   open(OLIGO, ">$tempSeq") or die "can not open $tempSeq for writing: $!\n";
   print OLIGO $sequence;
   close OLIGO;

   my ($dG, $align);
   system("./RNAfold < $tempSeq > $tempOUT");
   open(RESULT, "$tempOUT") or die "Cannot open $tempOUT for reading $!\n";
   while(my $line = <RESULT>){

         # ACCUUUUGUAUUUUAGUAACUGAAUCCCCACUGUGCAGUGUUAGGGCUGCCUGGUUGUUUGCAGUAGAUUAGAGCUUU
         # ..........(((((((.((((.((..(((..(.(((((......)))))))))..))...))))..))))))).... (-14.50)
         if($line =~ /([\-|\d][\.|\d]+)\)/){
             $dG = $1;
             if($line =~ /^([\.\(\)]+)\s+/){
                 $align = $1;
             }
             else {
                   print $line, "\tUm... RNAfold alignment is empty!\n"; exit;
             }
             last;
         }
   }
   close(RESULT);
   unlink $tempSeq;
   unlink $tempOUT;
   return ($dG, $align);
}

################################################################################################
# Subroutine: print_svm_features
# Converts result_01_features file to svm format and prints to file
# Input - the feature filename, the output filename
# Output - the header, also writes to OUT
################################################################################################
sub print_svm_features{
	my $in_file = shift;
	my $out_file = shift;
	open IN, $in_file; my $header = <IN>; $header =~ s/\s+$//;
	open OUT, ">$out_file";
	while (<IN>){
		$_ =~ s/\s+$//;
		my @line = split /\t/, $_;
		my $label = shift @line;
		my $svm_line = get_svm_features(\@line, $label);
		print OUT $svm_line . "\n";
	}
	close IN;
	close OUT;
	return $header;
}

# The parameters are from the following paper:
# Sugimoto N, Nakano S, Katoh M, Matsumura A, Nakamuta H, Ohmichi T, Yoneyama M, Sasaki M.
# Thermodynamic parameters to predict stability of RNA/DNA hybrid duplexes.
# Biochemistry. 1995 Sep 5;34(35):11211-6.
sub rna_dna_param {

    my %dH = ();
    my %dS = ();
    my %dG = ();

    ($dH{'aa'}, $dS{'aa'}, $dG{'aa'}) = (-11.5, -36.4, -0.2);
    ($dH{'tt'}, $dS{'tt'}, $dG{'tt'}) = (-7.8, -21.9, -1.0);
    ($dH{'at'}, $dS{'at'}, $dG{'at'}) = (-8.3, -23.9, -0.9);
    ($dH{'ta'}, $dS{'ta'}, $dG{'ta'}) = (-7.8, -23.2, -0.6);
    ($dH{'ca'}, $dS{'ca'}, $dG{'ca'}) = (-10.4, -28.4, -1.6);
    ($dH{'tg'}, $dS{'tg'}, $dG{'tg'}) = (-9.0, -26.1, -0.9);
    ($dH{'ct'}, $dS{'ct'}, $dG{'ct'}) = (-9.1, -23.5, -1.8);
    ($dH{'ag'}, $dS{'ag'}, $dG{'ag'}) = (-7.0, -19.7, -0.9);
    ($dH{'ga'}, $dS{'ga'}, $dG{'ga'}) = (-8.6, -22.9, -1.5);
    ($dH{'tc'}, $dS{'tc'}, $dG{'tc'}) = (-5.5, -13.5, -1.3);
    ($dH{'gt'}, $dS{'gt'}, $dG{'gt'}) = (-5.9, -12.3, -2.1);
    ($dH{'ac'}, $dS{'ac'}, $dG{'ac'}) = (-7.8, -21.6, -1.1);
    ($dH{'cg'}, $dS{'cg'}, $dG{'cg'}) = (-16.3, -47.1, -1.7);
    ($dH{'gc'}, $dS{'gc'}, $dG{'gc'}) = (-8.0, -17.1, -2.7);
    ($dH{'gg'}, $dS{'gg'}, $dG{'gg'}) = (-9.3, -23.2, -2.1);
    ($dH{'cc'}, $dS{'cc'}, $dG{'cc'}) = (-12.8, -31.9, -2.9);

    my ($dHi, $dSi, $dGi) = (1.9, -3.9, 3.1);

    return (\%dH, \%dS, \%dG, $dHi, $dSi, $dGi);
}

# calculate gc pecent range 0-1
sub gcPercent {

    my ($sequence) = @_;
    $sequence =~ tr/AUTCG/autcg/;
    my $length = length($sequence);
    my $gcCount= 0;

    for (my $j = 0; $j < $length; $j++) {
         if (substr($sequence, $j, 1) eq 'c' || substr($sequence, $j, 1) eq 'g') {
              $gcCount++;
         }
    }

    return sprintf("%.2f", $gcCount / $length);

}
