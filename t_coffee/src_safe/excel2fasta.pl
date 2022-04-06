#!/usr/bin/env perl

open (F, $ARGV[0]);

while ( <F>)
  {
    @l=($_=~/(\S+)/g);
    
    $name=shift @l;
    
    print STDOUT "\n>$name\n";
    foreach $e (@l){$e=($e eq "0")?"O":"I";print "$e";}
  }
close (F);

		       
    
# mar 05 abr 2022 14:52:23 CEST
