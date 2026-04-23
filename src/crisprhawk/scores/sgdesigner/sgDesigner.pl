#!/usr/bin/perl

#################################################################################################################

# Copyright (C) 2019  Xiaowei Wang (email: xiaowei.wang@wustl.edu), Washington University in St. Louis
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
#################################################################################################################

use strict;
use warnings;
use FindBin 1.51 qw( $RealBin );
use lib $RealBin;
use Cwd;
use format_features_25;

###################################### SETTINGS ############################################################
my $PYTHON_BIN = `python -c "import sys; print(sys.executable)"`;
chomp($PYTHON_BIN);

my $version =         "V2.0";
my $file_dir =        "./";
# my $result_dir =      "./result";
# system("mkdir $result_dir") if !-e $result_dir;
# my $temp_dir =        "./temp";
# system("mkdir $temp_dir") if !-e $temp_dir;
my $result_dir = $ENV{'SGDESIGNER_RESULT_DIR'} // "./result";
system("mkdir -p $result_dir") if !-e $result_dir;
my $temp_dir = $ENV{'SGDESIGNER_TEMP_DIR'} // "./temp";
system("mkdir -p $temp_dir") if !-e $temp_dir;

# print "RESULT DIR PRINT: $result_dir\n";
# print "TEMP DIR PRINT: $temp_dir\n";

my $classifier_dir =  "$MODEL_DIR";
my $classifier_dir =  "$MODEL_DIR";

my $inputFile;
# my $result_file =     "./$temp_dir/gOligo_$version"."_prediction_result.xls";
# my $feature_file =    "./$temp_dir/custom_features_v2.0.txt";
# my $predict_file =    "./$temp_dir/custom_prediction_result_v2.0.txt";
# my $outputFile = 	  "$result_dir/sgDesigner_$version"."_prediction_result.txt";
my $result_file =     "$temp_dir/gOligo_$version"."_prediction_result.xls";
my $feature_file =    "$temp_dir/custom_features_v2.0.txt";
my $predict_file =    "$temp_dir/custom_prediction_result_v2.0.txt";
my $outputFile =      "$result_dir/sgDesigner_$version"."_prediction_result.txt";

my $scaffold_seq = "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTT";

# position restrictions for selection of sgRNA
my $min_pos = 0;    # minimum position in CDS #remove cutoffs
my $max_pos = 100000000; # maximum position in CDS
my $retained_portion = 1; # retained region (70%) in the 5' portion of the CDS
my $minLength= 26;

my $leftflank_numbases = 0; 
my $rightflank_numbases = 3;

my $RNAFOLD_BIN = "RNAfold";  # default: assume in PATH
my $MODEL_DIR = "./Stacking_model";

################################### USER INPUTS ############################################################

my @inputs = @ARGV;
my $option = shift(@inputs);
&helpText if !@inputs;

my @sequences;
# print "Welcome to sgDesigner.\n\n";
if ($option eq '-e' or $option eq '--example'){
        my $sampleSelection = shift @inputs;
        
        die ("Please select a valid sample option.\n") unless $sampleSelection eq 'short' or $sampleSelection eq 'long' or $sampleSelection eq 'multiple';
        $inputFile = "./samples/test_sequence_short.fasta" if $sampleSelection eq 'short';
        $inputFile = "./samples/test_sequence_long.fasta" if $sampleSelection eq 'long';
        $inputFile = "./samples/test_sequence_multiple.fasta" if $sampleSelection eq 'multiple';
        # print "Selected file: $inputFile\n";
        @sequences = importFasta($inputFile);        
}
elsif($option eq '-f' or $option eq '--file'){
        $inputFile = shift @inputs;
        # print "Selected file: $inputFile\n";
        if (-e $inputFile && -r $inputFile && -f $inputFile && -T $inputFile){ #check to ensure file exists, readable, plain text
                open (INPUTCHECK,$inputFile);
                while (<INPUTCHECK>){
                        s/\s+$//;
                        die ("Please ensure that the file is in FASTA format.\n") if $_ !~ /^>/;
                        last;
                }
                close (INPUTCHECK);
                @sequences = importFasta($inputFile);                                                                            
        }else{
                print "Error: Please check to make sure the file \"$inputFile\" exists and is a readable plain text file.\n";
                exit;       
        }
}
elsif ($option eq '-s' or $option eq '--sequence'){
        my $submission = shift @inputs;
        ${$sequences[0]}{'seq'} = $submission;
        my $seqLength = length $submission;
        ${$sequences[0]}{'id'} = "submittedSequence|length_$seqLength";
}
elsif ($option eq '--rnafold') {
        $RNAFOLD_BIN = shift @inputs;
}
elsif ($option eq '--model') {
        $MODEL_DIR = shift @inputs;
}
else{
        &helpText;
}

################################### USER PREDICTION ########################################################

my $startTime = time();

# $submittedSeq =~ tr/ATCGU/atcgu/; $submittedSeq =~ tr/u/t/;
# print "\n******************** Genome-Wide sgOligo Version $version Standard Output **************************\n";

unlink $result_file if -e $result_file;
unlink $feature_file if -e $feature_file;
unlink $predict_file if -e $predict_file;
mkdir $result_dir if !(-e $result_dir);

open(FEAT, ">$feature_file");

my ($dH_ref, $dS_ref, $dG_ref, $dHi, $dSi, $dGi) = rna_dna_param();
my %dG = %{$dG_ref};

$scaffold_seq =~ tr/ATCGU/atcgu/;  $scaffold_seq =~ tr/t/u/;

my @feature;    # feature lines for SVM prediction
my @annotation; # annotation lines for mapping sgRNAs to gene annotations
my %oligo_pos; # to remove gene location redundancy from mulitple genomic loci
my $gene_count = 0;
my $line_count = 0;
my $total_line_count = 0;
my $submittedSeq;
my $id;
my %id2Sequence;
foreach (@sequences){
        my $sequence = ${$_}{'seq'};
        $id = ${$_}{'id'};
        $id2Sequence{$id}=$sequence; 
        $submittedSeq = $sequence;
        $submittedSeq =~ tr/ATCGU/atcgu/; $submittedSeq =~ tr/u/t/;

        print "Error: Sequence contains bases other than A, T, C, G, or U. \n\tWU-CRISPR will now proceed to the next sequence.\n\n" and next if $submittedSeq =~/[^atcg]/i;
        print "Error: Sequence is shorter than $minLength bases. \n\tWU-CRISPR will now proceed to the next sequence.\n\n" and next if length ($submittedSeq)<$minLength;
        print "Error: Sequence is longer than 100,000 bases. \n\tWU-CRISPR will now now proceed to the next sequence.\n\n" and next if length ($submittedSeq)>100000;

        my $submittedSeq_rc = dnaComplement($submittedSeq);

        generate_feature($submittedSeq,"sense");
        generate_feature($submittedSeq_rc,"antisense");
}
predict(\@feature, \@annotation);

close FEAT;
# print "\n************* sgOligo selection process is done. Program completed successfully. *******************\n";

################################### USER OUTPUTS ########################################################

open(RESULT, $result_file) or die $!;
my %resultSeqs;
my %scoreList;
my $seqId;
while (<RESULT>) {
        s/\s+$//;
        next if $_ =~/^Prediction_scores/;
        my @inline = split /\t/, $_;
        my $seq = substr($inline[6],0,20);
        my $seqSearch = $seq;
        my $orient = $inline[2];
        $seqId = $inline[1];


        $resultSeqs{$seq}{'orient'} = $orient;
        my $oligoSearch = $id2Sequence{$seqId};
        my $oligoLoc;
        $seqSearch = reverse($seq) and $seqSearch =~ tr/atcg/tagc/ if $orient eq 'antisense';
        my @pos1based;
        while ($oligoSearch =~ /$seqSearch/g) {
                $oligoLoc = pos($oligoSearch)-length($seqSearch)+1;
                push @pos1based,$oligoLoc;
        }
        
        my $location = join(", ",sort{$a<=>$b} @pos1based);
        $resultSeqs{$seq}{'location'} = $location;
        my $score= int($inline[0]+0.5);
        $scoreList{$seqId}{$seq} = $score;
}
close RESULT;

open(OUT, ">$outputFile") or die "$outputFile could not be opened for writing\n";

print OUT "seqId\tScore\tSequence\tOrientation\tPosition\n";
foreach my $sequenceId (sort keys %scoreList){
        foreach my $seq (sort {$scoreList{$sequenceId}{$b} <=> $scoreList{$sequenceId}{$a}} keys %{$scoreList{$sequenceId}}){                
                my $score = $scoreList{$sequenceId}{$seq};
                # next if $score<50;
                my $direction = $resultSeqs{$seq}{'orient'};
                my $position = $resultSeqs{$seq}{'location'};
                print OUT "$sequenceId\t$score\t$seq\t$direction\t$position\n";
        }
        
}
close OUT;

my $endTime = time();
my $finalTime = $endTime-$startTime;
# # print"\nResults have been printed in $outputFile. Program completed in $finalTime seconds.\n";

########################################################################################################################################################################
sub generate_feature {

    my ($exon_plus,$orientation) = @_;
    
    my $dummyString = 'n'x6;
    $exon_plus = join("",$dummyString,$exon_plus,$dummyString);

    my $exon = substr($exon_plus, 6, length($exon_plus) - 12);

    my $seqStrand = $exon_plus;
    $seqStrand = dnaComplement($exon_plus) if $orientation eq 'antisense'; #return to the original sequence

    my @gg_pos_list;
    for (my $i =0; $i<length($exon);$i++){
        my $dinuc = substr($exon,$i,2);
        push @gg_pos_list, $i if $dinuc eq 'gg';
   }
    foreach my $gg_pos (@gg_pos_list) {
 
            my $oligo = substr($exon_plus, $gg_pos + 6 - 21, 26);
            next if $oligo =~ /n/g;
            
            my $matched_bases = $oligo;
            $matched_bases = dnaComplement($oligo) if $orientation eq 'antisense';
            
            my $cds_pos = index($seqStrand,$matched_bases)+1;

            my ($feat_vals_ref, $feat_names_ref) = format_features_25::get_seq_features($oligo,$leftflank_numbases,$rightflank_numbases);
            my @feat_vals = @{$feat_vals_ref};
            my $output = join("\t",@feat_vals);
            print FEAT "$output\n";

            if ($output) {
                 $feature[$line_count] = $output;

                 my $exon_pos = $gg_pos;
                 $exon_pos = length($exon) - $gg_pos if $orientation eq "antisense";
                 $exon_pos += 1; # 1-based index position
                #  print "$exon_pos\n";

                 $annotation[$line_count] = "$id\t$orientation\t$exon_pos\t$cds_pos\t".length($exon_plus)."\t$oligo";

                 $line_count++;
            }
            $total_line_count++;
    }
}


sub predict {

    my ($feature_ref, $annotation_ref) = @_;
    my @annotation = @{$annotation_ref};
    my $line_count = 0;

    my $exit_code = system("cd $MODEL_DIR && $PYTHON_BIN $classifier_dir/Stacking_classification.py");

    if ($exit_code != 0) {
        die "ERROR: Python script failed with exit code $exit_code. Check if model files exist in $classifier_dir\n";
    }

    open(OUT, ">$result_file");
    print OUT "Prediction_scores\tsequenceID\tOrientation\tPosition in Exon\tPosition in CDS\tCDS Length\tOligo Sequence\n";

    open(IN, "$predict_file");<IN>;
    while(<IN>){
      $_ =~ s/\s+$//;
      my @line = split /\t/, $_;
      my $out = $line[1]."\t".$annotation[$line_count];
      print OUT "$out\n";
      $line_count ++;
    }
    close IN;
    close OUT;
}

# The overall deltaG of the target binding duplex. This is similar to GC content
sub dG_binding {

    my $oligo = shift;
    $oligo =~ s/u/t/g;
    my $binding_dG = 0;
    for (my $indx = 0; $indx < length($oligo) - 1; $indx++) {
        next if substr($oligo, $indx, 2) =~ /n/;
                print (substr($oligo, $indx, 2),"\n") if !exists($dG{substr($oligo, $indx, 2)});

         $binding_dG += $dG{substr($oligo, $indx, 2)};
    }
    $binding_dG += $dGi;

    return $binding_dG;
}

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
   system("$RNAFOLD_BIN < $tempSeq > $tempOUT");
   open(RESULT, "$tempOUT") or die "Cannot open $tempOUT for reading $!\n";
   while(my $line = <RESULT>){
      # .((((((((((((((((((((((((((....)))))))))))))))))))))))))) (-17.80)
      print $line, "\n";
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

# Calculate the secondary structure delta G for a single input sequence
# Input - a single RNA oligo sequence
# Output - the deltaG value
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
   system("$RNAFOLD_BIN < $tempSeq > $tempOUT");
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


# return the self-complementary strand of the input sequence
sub dnaComplement {

    my ($sequence) = @_;
    $sequence =~ tr/atcgATCG/tagcTAGC/;
    $sequence = reverse($sequence);
    return $sequence;
}

sub importFasta {
    my ($fastaFile) = @_;
    my $tabFile = "$fastaFile $$.tab";
    fastaToTab($fastaFile, $tabFile);
    my @seq = importTabSeq($tabFile);
    unlink $tabFile if -e $tabFile;
    return @seq;
}

sub fastaToTab {
     my ($fastaFile, $tabFile) = @_;
     my $id = "";
     my $dna= "";
     my $lastLine = "";
     open (IN, "$fastaFile") || die("Can not open $fastaFile file for reading in fastaToTab sub!\n");
     open (OUT, ">$tabFile") || die("Can not open $tabFile file for writing!\n");

     while (<IN>) {
          s/\s+$//;
          next if ($_ !~ /\S/);
          if ($_ =~ /^\>/) {
               $id = $_;
               $id =~ s/^\>//;
               if ($lastLine =~ /^\>/) {   
                    $id .= $_;
               }
               else {
                    print OUT "\n" if ($dna ne ""); 
                    print OUT "$id\t";
                    $id = "";
               }
          }
          else {
               $_ =~ s/\s//g;
               $dna = $_;
               print OUT $dna;
          }
          $lastLine = $_;
     }
     close(IN);
     close(OUT);
}

sub importTabSeq {
     my ($tabFile) = @_;
     my @sequence = ();
     my $index = 0;
     open (IN, "$tabFile") || die("Cannot open $tabFile file for reading in importTab sub!\n");
     while (<IN>) {
          s/\s+$//;
          my ($id, $sequence) = split /\t/, $_;
          $sequence[$index]{'id'} = $id;
          $sequence =~ tr/A-Z/a-z/;
          $sequence[$index]{'seq'} = $sequence;
          $index++;
     }
     close(IN);
     return @sequence;
}

sub helpText{
        # print "\n";
        
        # print "USAGE:\n\tperl sgDesigner.pl [option] [path]\n\tperl sgDesigner.pl [option] [sequence]\n\n";
        # print "SEQUENCE SUBMISSION:\n\t-s|--sequence <sequence>\n\t\t";
        # # print "Identifies sgRNA oligos from a single submitted sequence \n\t\tand provides a score for all potential active oligos. \n\t\tResults for submitted sequences will be printed to a \n\t\ttab-delimited text file.\n";
        # print "\n\t\tExample: perl sgDesigner.pl -s acctgcgtggctcccctgagtggagt\n\n";
        
        # print "FILE SUBMISSION:\n\t-f|--file <file>\n\t\t";
        # print "Imports a FASTA file of sequences from <file> and \n\t\tidentifies potential sgRNA oligos for each submitted\n\t\tsequence. Resulting oligos are available in a tab-\n\t\tdelimited text file.\n";
        # print "\n\t\tExample: perl sgDesigner.pl -f mySampleFile.fasta\n\n";
        
        # print "EXAMPLE FILE SUBMISSION:\n\t-e|--example <short|long|multiple>\n\t\t";
        # print "Uses the short, long, or multiple sample sequence in the \n\t\tsamples directory to generate sgRNA oligos.\n";
        # print "\n\t\tExample: perl sgDesigner.pl -e short\n\n";
        
        # print "HELP SCREEN:\n\t-h|--help \n\t\t";
        # print "Brings up this help menu";
        # print "\n\t\tExample: perl sgDesigner.pl -h\n\n";
        
        exit;
}