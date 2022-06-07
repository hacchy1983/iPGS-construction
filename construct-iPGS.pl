#!/usr/bin/perl
BEGIN{ push(@INC, $ENV{LOCAL_PERL_MODULES}); }
use strict;
use warnings;
use Getopt::Std;
use vars qw($opt_h $opt_i $opt_s $opt_p $opt_o);
&getopts('hi:s:p:o:');

my $usage = <<_EOH_;

## Options ###########################################
## Required:
# -i    input mdoel weight file
# -s    input model scale file
# -p    input variant ID file 
# -o    output model file

## Optional:
# 

## Others:
# -h print help

_EOH_
    ;

die $usage if $opt_h;

# Get command line options
my $betaFile  = $opt_i or die $usage;
my $scaleFile = $opt_s or die $usage;
my $snpFile   = $opt_p or die $usage;
my $outFile   = $opt_o or die $usage;

# Load model beta & scale
print STDERR "Loading model weights & scaling factors...\n";
my $model2beta  = loadModelParam( $betaFile );
my $model2scale = loadModelParam( $scaleFile );
my $model2w     = calcModelWeight( $model2beta, $model2scale );
print STDERR "...done.\n\n";

# Load SNP list
print STDERR "Loading SNP ID list...\n";
my $SNPs = loadSNPs( $snpFile );
print STDERR "...done.\n\n";

# Construct iPGS
print STDERR "Constructing iPGS model...\n";
my $SNP2w = construct_iGRS( $SNPs, $model2w );
print STDERR "...done.\n\n";

# Output
print STDERR "Saving model data...\n";
saveModel( $outFile, $SNPs, $SNP2w );
print STDERR "...done.\n\n";

sub loadModelParam {
    my $file = shift;
    my $model2param = {};
    open IN, $file or die "Can't open the file [ $file ] to read.";
    while(<IN>) {
	chomp;
	my @cols = split(/\s+/, $_);
	if( scalar(@cols) == 2 ) {
	    my ($model, $param) = @cols;
	    $model2param->{$model} = $param;
	}
    }
    close IN;
    print STDERR "\tParameters for " . scalar(keys %$model2param) . " models were loaded from the file [ $file ].\n";
    return $model2param;
}

sub calcModelWeight {
    my $model2beta  = shift;
    my $model2scale = shift;
    my $model2w     = {};
    foreach my $model (sort keys %$model2beta) {
	my $beta  = $model2beta->{$model};
	my $scale = $model2scale->{$model};
	next if ! defined $scale;
	$model2w->{$model} = $beta / $scale;
    }
    print STDERR "\tWeights for " . scalar(keys %$model2w) . " models were calculated:\n";
    foreach my $model (sort keys %$model2w) {
	my $w = $model2w->{$model};
	print STDERR "\t\t" . $model . "\t" . $w . "\n";
    }
    return $model2w;
}

sub loadSNPs {
    my $file = shift;
    my $SNPs = [];
    open IN, $file or die "Can't open the file [ $file ] to read.";
    while(<IN>) {
	chomp;
	my @cols = split("\t", $_);
	if( scalar(@cols) == 1 ) {
	    my $SNP = $cols[0];
	    push(@$SNPs, $SNP);
	}
    }
    close IN;
    print STDERR "\t" . scalar(@$SNPs) . " variants were loaded from the file [ $file ].\n";
    return $SNPs;
}

sub construct_iGRS {
    my $SNPs    = shift;
    my $model2w = shift;
    my $SNP2w   = initSnpWeights( $SNPs );
    foreach my $model (sort keys %$model2w) {
	my $modelWeight = $model2w->{$model};
	print STDERR "\t\tUpdate model weight: $model (weight = $modelWeight) ... ";
	my $numUpdates = updateSnpWeights( $SNP2w, $model, $modelWeight );
	print STDERR "weights for " . $numUpdates . " variants were updated.\n";
    }
    return $SNP2w;
}

sub initSnpWeights {
    my $SNPs = shift;
    my $SNP2w = {};
    foreach my $SNP (@$SNPs) {
	$SNP2w->{$SNP} = 0;
    }
    print STDERR "\t" . "Weights for " . scalar(keys %$SNP2w) . " variants were initialized.\n";
    return $SNP2w;
}

sub updateSnpWeights {
    my $SNP2w  = shift;
    my $file   = shift;
    my $modelW = shift;
    if( $file =~ /\.gz$/ ) {
	open IN, "zcat $file 2>/dev/null |" or die "Can't open the file [ $file ] to read.";
    } else {
	open IN, $file or die "Can't open the file [ $file ] to read.";
    }
    my $header     = <IN>; chomp $header;
    my $var2hidx   = parseHeader( $header );
    my $numUpdates = 0;
    while(<IN>) {
	chomp;
    	my @cols = split(/\s+/, $_);
    	if( scalar(@cols) == scalar(keys %$var2hidx) ) {
	    my $chr  = $cols[$var2hidx->{chr_name}];
	    my $bp   = $cols[$var2hidx->{chr_position}];
	    my $ref  = $cols[$var2hidx->{reference_allele}];
	    my $alt  = $cols[$var2hidx->{effect_allele}];
	    my $beta = $cols[$var2hidx->{effect_weight}];
	    my $ea   = $alt;
	    my $oa   = $ref;

	    next if $ea eq "A" and $oa eq "T";
	    next if $ea eq "T" and $oa eq "A";
	    next if $ea eq "G" and $oa eq "C";
	    next if $ea eq "C" and $oa eq "G";
	    next if $beta eq "NA";
	    next if $beta eq "NULL";

	    my $ID1 = $chr . ':' . $bp . ':' . $oa . ':' . $ea;
	    my $ID2 = $chr . ':' . $bp . ':' . cmpl($oa) . ':' . cmpl($ea);
	    my $ID3 = $chr . ':' . $bp . ':' . $ea . ':' . $oa;
	    my $ID4 = $chr . ':' . $bp . ':' . cmpl($ea) . ':' . cmpl($oa);
	    if( defined $SNP2w->{$ID1} ) {
		$SNP2w->{$ID1} += $modelW * $beta;
		++$numUpdates;
	    } elsif( defined $SNP2w->{$ID2} ) {
		$SNP2w->{$ID2} += $modelW * $beta;
		++$numUpdates;
	    } elsif( defined $SNP2w->{$ID3} ) {
		$SNP2w->{$ID3} += -1 * $modelW * $beta;
		++$numUpdates;
	    } elsif( defined $SNP2w->{$ID4} ) {
		$SNP2w->{$ID4} += -1 * $modelW * $beta;
		++$numUpdates;
	    }
	}
    }
    close IN;
    return $numUpdates;
}

sub parseHeader {
    my $header = shift;
    my @hcols = split(/\s+/, $header);
    my $var2hidx = {};
    #print STDERR "\tParsing header:\n";
    for( my $i=0 ; $i<scalar(@hcols) ; ++$i ) {
	$var2hidx->{$hcols[$i]} = $i;
	#print STDERR "\t\t$i\t$hcols[$i]\n";
    }
    return $var2hidx;
}

sub capitalizeBase {
    my $base = shift;
    return "A" if $base eq "a";
    return "G" if $base eq "g";
    return "T" if $base eq "t";
    return "C" if $base eq "c";
    return $base;
}

sub saveModel {
    my $file  = shift;
    my $SNPs  = shift;
    my $SNP2w = shift;
    open OUT, ">$file" or die "Can't open the file [ $file ] to write.";
    print OUT "chr_name\tchr_position\treference_allele\teffect_allele\teffect_weight\n";
    my $cnt = 0;
    foreach my $SNP (@$SNPs) {
	my $w = $SNP2w->{$SNP};
	next if ! defined $w;
	next if $w == 0;
	my ($chr, $bp, $ref, $alt) = split(':', $SNP);
	print OUT $chr . "\t"
	    . $bp . "\t"
	    . $ref . "\t"
	    . $alt . "\t"
	    . $w . "\n";
	++$cnt;
    }
    close OUT;
    print STDERR "\t" . "Model data for " . $cnt . " variants were saved to the file [ $file ].\n";
}

sub cmpl {
    my $b = shift;
    return "A" if $b eq "T";
    return "T" if $b eq "A";
    return "G" if $b eq "C";
    return "C" if $b eq "G";
    return "NA";
}
