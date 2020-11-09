#!/usr/bin/perl -w
use FindBin qw($Bin);
use strict;
#explanation:this program is edited to 
#edit by hewm;   Mon Jan 26 14:37:33 CST 2015
#Version 1.0    hewm@genomics.org.cn 

die  "Version 1.0\t2015-01-26;\nUsage: $0 <DepthGC.wig.gz><MaxDepth><OUT_prex>\n" unless (@ARGV ==3);

#############Befor  Start  , open the files ####################

######################swimming in the sky and flying in the sea ###########################
my $file= $ARGV[0] ;

if  ($file =~s/\.gz$/\.gz/)
{
	open FIN,"gzip -cd  $file | "  || die "input file can't open $!" ;
}
else
{
	open FIN,"$file"  || die "input file can't open $!" ;
}


open (OAtmp , ">$ARGV[2].tmp") || die " output file can't open $!";
open (OAR , ">$ARGV[2].r") || die " output file can't open $!";


my %DD=();
while ($_=<FIN>)
{
	chomp ;
	next if  ($_=~s/#/#/g);
	next if  ($_=~s/NA/NA/g);
	my @inf=split ;
	if ($inf[1]>$ARGV[1])
	{
		$inf[1]=$ARGV[1];
	}
	$DD{$inf[-1]}{$inf[1]}++;
}
close FIN ;


print OAtmp  "GC_ratio\tDepth\tNumbers\n";
foreach my $k1 ( sort {$a<=>$b} keys %DD)
{
	my $sedH=$DD{$k1};
	foreach my $k2 ( sort {$a<=>$b} keys %$sedH)
	{
		print OAtmp  "$k1\t$k2\t$DD{$k1}{$k2}\n";
	}
}
close OAtmp ;


my $Pit = <<PIG;

library(ggplot2)
library(gridExtra)
library(ggExtra)


read.table("$ARGV[2].tmp",header=T)->tab;
GC_ratio<-tab[,1]
Depth<-tab[,2]
Numbers<-tab[,3]

pdf("$ARGV[2].pdf");

df <- data.frame(GC_ratio , Depth )
p <- ggplot(df, aes(GC_ratio, Depth,colour=Numbers)) + geom_point(size=1.5) + theme_classic()
ggExtra::ggMarginal(p, type = "histogram")


dev.off();

png("$ARGV[2].png");

df <- data.frame(GC_ratio , Depth )
p <- ggplot(df, aes(GC_ratio, Depth,colour=Numbers)) + geom_point(size=1.5) + theme_classic()
ggExtra::ggMarginal(p, type = "histogram")


dev.off();
PIG
print OAR  $Pit ;
close OAR ;
my $Rscript="/usr/bin/Rscript";
if  ( !(-e $Rscript) )
{
	$Rscript="/ifswh1/BC_PUB/biosoft/BC_NQ/01.Soft/03.Soft_ALL/R-3.4.1/bin/Rscript";
	if  ( !(-e $Rscript) )
	{
		$Rscript=`which Rscript`;chomp $Rscript;
	}
}




if  ( !(-e $Rscript) )
{
	print "Can't find the [ Rscript ] bin, You shoud install the R First,then:\n";
	print "  Rscript   $ARGV[2].r  \n";
	exit(1);
}


system ("  $Rscript   $ARGV[2].r  ");
system ("  rm -rf  $ARGV[2].tmp  $ARGV[2].r ") ;

