#!/usr/bin/perl

#Author: Sarah Rennie, 2019, (parts adapted from https://github.com/riverlee/pileup2base)
#pulls piped output from mpile to generate potential editing hits (to be formally tested)
#this is the version of the script for stranded data

use strict;
use warnings;
use Getopt::Long;


  my %reverse_hash = ( A  => "a",
    T  => "t",
    C  => "c",
    G => "g",
    a => "A",
    t => "T",
    c => "C",
    g => "G");

while (<>) 
{
  
  my $input = $_;
  #print $input."\n";
  $input=~s/(\s\s\s+)/ * * /g;
  
  my ($chr,$loc,$ref,@line) = split (/\s+/,$input);
  
  my $str="$chr\t$loc"."\t".$ref;
  
  #print $str."\n";
  
  my $ref_reverse = $reverse_hash{$ref};
  
  my $flag_A=0;
  my $flag_T=0;
  my $flag_C=0;
  my $flag_G=0;
  
  for (my $i = 0; $i < @line; $i += 3) {
    my $dp = $line[$i];
    my $base_line = $line[$i+1];
    my $q_scores = $line[$i+2];

          if($dp==0)
          {
             #print "is.zero\n";
             	$str="\t".$str."\t0\t0\t0\t0\t0\t0\t0\t0";
	      	next;
          }
    
        	if( $base_line eq "*" | $base_line eq ""){
             	$str="\t".$str."\t0\t0\t0\t0\t0\t0\t0\t0";
	      	next;	     
	            }
    
    $base_line=~s/\^.//g;
    $base_line=~s/\$//g;
    
    
    
  #remove deletions and insertions (next few lines adapted from https://github.com/riverlee/pileup2base)
  my %hash=();
	my %deletion=();
	while($base_line=~/-(\d+)/g){
		$hash{$1}=1;
	}
	
	
	foreach my $key (keys %hash){
		while($base_line=~/-$key([ACGTNacgtn]{$key})/g){
			$deletion{$1}++;
		}
		$base_line=~s/-$key[ACGTNacgtn]{$key}//g;
	}
	
	%hash=();
	my %insertion=();
	while($base_line=~/\+(\d+)/g){
		$hash{$1}=1;
	}
	foreach my $key (keys %hash){
		while($base_line=~/\+$key([ACGTNacgtn]{$key})/g){
			$insertion{$1}++;
		}
		$base_line=~s/\+$key[ACGTNacgtn]{$key}//g;
	}
    
  ########################################
    
   my %nts = ( A  => 0,
    T => 0,
    C  => 0,
    G => 0,
    a => 0,
    t => 0,
    c => 0,
    g => 0,);
    
     
     my @base=split (//,$base_line);
     my @scores=split(//,$line[$i+2]);
     

     #print scalar @base."\t". scalar @scores."\n";
    

     for(my $j=0; $j<@base; $j++){
       my $nt=$base[$j];
       my $score=ord($scores[$j])-33;
       if($score>0)
       {
          if($nt eq "."){
             $nts{$ref}++
           }elsif($nt eq ","){
               $nts{$ref_reverse}++
         }else{
          $nts{$nt}++
         }
       }

     }
   
   #print $nts{A}."\n";
   
   my $A = $nts{A} + $nts{a};
   my $T = $nts{T} + $nts{t};
   my $C = $nts{C} + $nts{c};
   my $G = $nts{G} + $nts{g};
   
   if($A >= 1)
   {
     $flag_A++;
   }
   if($T >= 1)
   {
     $flag_T++;
   }
   if($C >= 1)
   {
     $flag_C++;
   }
   if($G >= 1)
   {
     $flag_G++;
   }
   
   $str="\t".$str."\t".$nts{A}."\t".$nts{T}."\t".$nts{C}."\t".$nts{G}."\t".$nts{a}."\t".$nts{t}."\t".$nts{c}."\t".$nts{g};
	   
   }
  
  if($flag_A>=2)
  {
    $flag_A=1;
  }else{
    $flag_A=0;
  }
  if($flag_T>=2)
  {
    $flag_T=1;
  }else{
    $flag_T=0;
  }
  if($flag_C>=2)
  {
    $flag_C=1;
  }else{
    $flag_C=0;
  }
  if($flag_G>=2)
  {
    $flag_G=1;
  }else{
    $flag_G=0;
  }
  my $flag_sum = $flag_A + $flag_T + $flag_C + $flag_G;
  if($flag_sum>1)
  {
    #print $flag_sum."\n";
    print $str."\n";
  }
}