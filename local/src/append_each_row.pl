#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
$\="\n";

my $usage ="$0 [-A|-B] [-F '-'] 'word' file";

my $A=1;
my $B=0;
my $F="\t";

GetOptions (
	'before|b' => \$B,
	'after|a' => \$A,
	'field-seaprator|F=s' => \$F,
) or die($usage);

my $a=shift @ARGV;
$a=~s/\\t/\t/g;

$A=0 if $B;

die("invalid argument") if !defined($a) or $a eq '';

while(<>){
	chomp;
	print $a,$F,$_ if $B;
	print $_,$F,$a if $A;
}
