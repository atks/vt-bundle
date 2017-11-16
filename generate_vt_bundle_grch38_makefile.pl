#!/usr/bin/perl -w

use warnings;
use strict;
use POSIX;
use Getopt::Long;
use File::Path;
use File::Basename;
use Pod::Usage;

=head1 NAME

generate_vt_bundle_grch38_makefile

=head1 SYNOPSIS

 generate_vt_bundles_grch38_makefile [options]

 --cluster cluster (default main)
  -m       output make file

 example: codes/bin/generate_vt_bundles_makefile -m make_ref.mk --chr 20 --cluster 1000g

=head1 DESCRIPTION

=cut

#option variables
my $help;
my $verbose;
my $debug;
my $makeFile = "create_bundle_grch38.mk";
my $ext = "bcf";

#initialize options
Getopt::Long::Configure ('bundling');

if(!GetOptions ('h'=>\$help, 'v'=>\$verbose, 'd'=>\$debug,
                'e:s' => \$ext)
  || scalar(@ARGV)!=0)
{
    if ($help)
    {
        pod2usage(-verbose => 2);
    }
    else
    {
        pod2usage(1);
    }
}

my @tgts = ();
my @deps = ();
my @cmds = ();
my $cmd;

my $tgt;
my $dep;
my @cmd;

#directories
my $assembly = "grch38";
my $homeDir = "/net/fantasia/home/atks";
my $binDir = "$homeDir/dev/vt/bundle/bin";
my $outputDir = "$homeDir/dev/vt/bundle/public/$assembly";
my $grch37OutputDir = "$homeDir/dev/vt/bundle/public/grch37";
my $grch38OutputDir = "$homeDir/dev/vt/bundle/public/grch38";
my $logDir = "$homeDir/dev/vt/bundle/log/$assembly";
my $grch37LogDir = "$homeDir/dev/vt/bundle/log/grch37";
my $grch38LogDir = "$homeDir/dev/vt/bundle/log/grch38";
my $auxDir = "$homeDir/dev/vt/bundle/aux/$assembly";
my $grch37AuxDir = "$homeDir/dev/vt/bundle/aux/grch37";
my $grch38AuxDir = "$homeDir/dev/vt/bundle/aux/grch38";

mkpath($outputDir);
mkpath($logDir);
mkpath($auxDir);

#programs
my $vt = "/net/fantasia/home/atks/vt/vt";
my $bedtools = "/net/fantasia/home/atks/programs/bedtools2/bin/bedtools";
my $picard = "java -jar /net/fantasia/home/atks/programs/picard-tools-2.14.0/picard.jar";
my $mdust = "$homeDir/programs/mdust/mdust";

#reference sequence
my $refFASTAFile = "/net/fantasia/home/atks/ref/genome/hs38DH.fa";

#liftover from hg19
my $hg19ToHg38ChainFile = "/net/fantasia/home/atks/ref/ucsc/hg19/hg19ToHg38.over.chain.gz";
my $hs37ToHs38ChainFile = "/net/fantasia/home/atks/dev/vt/bundle/bin/chain/hs37ToHs38.over.chain.gz";

my $indexExt = $ext eq "bcf" ? "csi" : "tbi";

my $srcVCFFile;
my $data;
my $destVCFFile;
my $destSitesVCFFile;
my $rejectedVCFFile;
my $rejectedSitesVCFFile;
my $srcBEDFile;
my $srcUCSCSQLTableFile;
my $type;

###################
#GENCODE Annotation
###################
my $srcGTFFile = "/net/fantasia/home/atks/ref/encode/grch38/gencode.v27.annotation.gtf.gz";
my $destGTFFile = "gencode.v27.annotation.gtf.gz";

#sort and bgzip
$tgt = "$logDir/$destGTFFile.OK";
$dep = "$srcGTFFile";
@cmd = ("$binDir/clean_gencode_annotation $srcGTFFile | grep -v \"^#\" | sort -V -k1,1 -k4,5 | bgzip -c > $outputDir/$destGTFFile");
makeStep($tgt, $dep, @cmd);

#index
$tgt = "$logDir/$destGTFFile.tbi.OK";
$dep = "$logDir/$destGTFFile.OK";
@cmd = ("tabix -pgff $outputDir/$destGTFFile");
makeStep($tgt, $dep, @cmd);

####################
#GENCODE CDS regions
####################
$srcGTFFile = "/net/fantasia/home/atks/ref/encode/grch38/gencode.v27.annotation.gtf.gz";
my $destBEDFile = "gencode.v27.cds.bed.gz";

#sort and bgzip
$tgt = "$logDir/$destBEDFile.OK";
$dep = "$srcGTFFile";
@cmd = ("zcat $srcGTFFile | grep -v \"^#\" | sort -V -k1,1 -k4,5 | cut -f1,3,4,5 | grep CDS | cut -f1,3,4 | perl -lane '{--\$\$F[1]; print join(\"\\t\", \@F);}' | tr '\\t' ':' | uniq | tr ':' '\\t' | bgzip -c > $outputDir/$destBEDFile");
makeStep($tgt, $dep, @cmd);

#index
$tgt = "$logDir/$destBEDFile.tbi.OK";
$dep = "$logDir/$destBEDFile.OK";
@cmd = ("tabix -pbed $outputDir/$destBEDFile");
makeStep($tgt, $dep, @cmd);

###################
#centromere regions
###################
$srcUCSCSQLTableFile = "/net/fantasia/home/atks/ref/ucsc/hg38/centromeres.txt.gz";
$destBEDFile = "centromeres.bed.gz";

#sort and bgzip
$tgt = "$logDir/$destBEDFile.OK";
$dep = "$srcUCSCSQLTableFile";
@cmd = ("zcat $srcUCSCSQLTableFile | cut -f2-4 | " .
        " perl -lane '{if (!/_/) {print \"\$\$F[0]\\t\$\$F[1]\\t\$\$F[2]\\n\";}}' | " .
        " $bedtools sort -i - | " .
        " sort -V | " .
        " bgzip -c > $outputDir/$destBEDFile");
makeStep($tgt, $dep, @cmd);

#index
$tgt = "$logDir/$destBEDFile.tbi.OK";
$dep = "$logDir/$destBEDFile.OK";
@cmd = ("tabix -pbed $outputDir/$destBEDFile");
makeStep($tgt, $dep, @cmd);

########################
#pseudoautosomal regions
########################
$destBEDFile = "pseudoautosomal.bed.gz";

#sort and bgzip
$tgt = "$logDir/$destBEDFile.OK";
$dep = "";
@cmd = ("$binDir/create_grch38_pseudoautosomal_regions_bed -o $outputDir/$destBEDFile");
makeStep($tgt, $dep, @cmd);

#index
$tgt = "$logDir/$destBEDFile.tbi.OK";
$dep = "$logDir/$destBEDFile.OK";
@cmd = ("tabix -pbed $outputDir/$destBEDFile");
makeStep($tgt, $dep, @cmd);

##############
#mdust regions
##############
$refFASTAFile = "/net/fantasia/home/atks/ref/genome/hs38DH.fa";
$destBEDFile = "mdust.bed.gz";
my $dustBEDFile = "$outputDir/$destBEDFile";

#sort and bgzip
$tgt = "$logDir/$destBEDFile.OK";
$dep = "$refFASTAFile";
@cmd = ("$binDir/mdust/mdust $refFASTAFile -c | cut -f1,3,4 | perl -lane '{--\$F[1]; print join(\"\\t\", \@F);}' | bgzip -c > $outputDir/$destBEDFile");
makeStep($tgt, $dep, @cmd);

#index
$tgt = "$logDir/$destBEDFile.tbi.OK";
$dep = "$logDir/$destBEDFile.OK";
@cmd = ("tabix -pbed $outputDir/$destBEDFile");
makeStep($tgt, $dep, @cmd);

#####################
#repeatmasker regions
#####################
$srcUCSCSQLTableFile = "/net/fantasia/home/atks/ref/ucsc/hg38/rmsk.txt.gz";
$destBEDFile = "rmsk.bed.gz";

#sort and bgzip
$tgt = "$logDir/$destBEDFile.OK";
$dep = "$srcUCSCSQLTableFile";
@cmd = ("zcat $srcUCSCSQLTableFile | cut -f6-8 | " .
        " perl -lane '{if (!/_/) {print \"\$\$F[0]\\t\$\$F[1]\\t\$\$F[2]\\n\";}}' | " .
        " $bedtools sort -i - | " .
        " sort -V | " .
        " bgzip -c > $outputDir/$destBEDFile");
makeStep($tgt, $dep, @cmd);

#index
$tgt = "$logDir/$destBEDFile.tbi.OK";
$dep = "$logDir/$destBEDFile.OK";
@cmd = ("tabix -pbed $outputDir/$destBEDFile");
makeStep($tgt, $dep, @cmd);

#############
##trf regions
#############
#$data = "trf";
#$srcVCFFile = "$data.sites.$ext";
#$destBEDFile = "$data.bed.gz";
#
##convert from VCF to bed
#$tgt = "$grch38LogDir/$destBEDFile.OK";
#$dep = "$grch38LogDir/$srcVCFFile.OK";
#@cmd = ("$binDir/convert_vcf_2_bed $grch38OutputDir/$srcVCFFile -o $grch38OutputDir/$destBEDFile");
#makeStep($tgt, $dep, @cmd);
#
##index
#$tgt = "$grch38LogDir/$destBEDFile.tbi.OK";
#$dep = "$grch38LogDir/$destBEDFile.OK";
#@cmd = ("tabix -pbed $grch38OutputDir/$destBEDFile");
#makeStep($tgt, $dep, @cmd);
#
#####################
##trf reference panel
#####################
#$data = "trf";
#
##convert to vcf.gz as picard cannot read bcf files
#$tgt = "$grch38LogDir/$data.grch37.sites.vcf.gz.OK";
#$dep = "$grch37LogDir/$data.sites.$ext.OK";
#@cmd = ("vt view -h $grch37OutputDir/$data.sites.$ext -o $grch38AuxDir/$data.grch37.sites.vcf.gz");
#makeStep($tgt, $dep, @cmd);
#
##liftover
#$tgt = "$grch38AuxDir/$data.sites.vcf.gz.OK";
#$dep = "$refFASTAFile $grch38LogDir/$data.grch37.sites.vcf.gz.OK";
#@cmd = ("$picard LiftoverVcf CHAIN=$hs37ToHs38ChainFile " .
#                           " R=$refFASTAFile " .
#                           " I=$grch38AuxDir/$data.grch37.sites.vcf.gz " .
#                           " O=$grch38AuxDir/$data.sites.vcf.gz " .
#                           " REJECT=$grch38AuxDir/$data.rejected.sites.vcf.gz " .
#                           " 2> $grch38LogDir/$data.liftover.log");
#makeStep($tgt, $dep, @cmd);
#
##convert to ext
#$tgt = "$grch38LogDir/$data.sites.$ext.OK";
#$dep = "$grch38AuxDir/$data.sites.vcf.gz.OK";
#@cmd = ("$vt view $grch38AuxDir/$data.sites.vcf.gz -o $grch38OutputDir/$data.sites.$ext");
#makeStep($tgt, $dep, @cmd);
#
##index sites file
#$tgt = "$grch38LogDir/$data.sites.$ext.$indexExt.OK";
#$dep = "$grch38LogDir/$data.sites.$ext.OK";
#@cmd = ("$vt index $grch38OutputDir/$data.sites.$ext 2> $grch38LogDir/$data.sites.$ext.index.log");
#makeStep($tgt, $dep, @cmd);

###################
#trf/lobstr regions
###################
$data = "trf.lobstr";
$srcVCFFile = "$data.sites.$ext";
$destBEDFile = "$data.bed.gz";

#convert from VCF to bed
$tgt = "$grch38LogDir/$destBEDFile.OK";
$dep = "$grch38LogDir/$srcVCFFile.OK";
@cmd = ("$binDir/convert_vcf_2_bed $grch38OutputDir/$srcVCFFile -o $grch38OutputDir/$destBEDFile");
makeStep($tgt, $dep, @cmd);

#index
$tgt = "$grch38LogDir/$destBEDFile.tbi.OK";
$dep = "$grch38LogDir/$destBEDFile.OK";
@cmd = ("tabix -pbed $grch38OutputDir/$destBEDFile");
makeStep($tgt, $dep, @cmd);

###########################
#trf/lobstr reference panel
###########################
$data = "trf.lobstr";

#convert to vcf.gz as picard cannot read bcf files
$tgt = "$grch38LogDir/$data.grch37.sites.vcf.gz.OK";
$dep = "$grch37LogDir/$data.sites.$ext.OK";
@cmd = ("vt view -h $grch37OutputDir/$data.sites.$ext -o $grch38AuxDir/$data.grch37.sites.vcf.gz");
makeStep($tgt, $dep, @cmd);

#liftover
$tgt = "$grch38AuxDir/$data.sites.vcf.gz.OK";
$dep = "$refFASTAFile $grch38LogDir/$data.grch37.sites.vcf.gz.OK";
@cmd = ("$picard LiftoverVcf CHAIN=$hs37ToHs38ChainFile " .
                           " R=$refFASTAFile " .
                           " I=$grch38AuxDir/$data.grch37.sites.vcf.gz " .
                           " O=$grch38AuxDir/$data.sites.vcf.gz " .
                           " REJECT=$grch38AuxDir/$data.rejected.sites.vcf.gz " .
                           " 2> $grch38LogDir/$data.liftover.log");
makeStep($tgt, $dep, @cmd);

#convert to ext
$tgt = "$grch38LogDir/$data.sites.$ext.OK";
$dep = "$grch38AuxDir/$data.sites.vcf.gz.OK";
@cmd = ("$vt view $grch38AuxDir/$data.sites.vcf.gz -o $grch38OutputDir/$data.sites.$ext");
makeStep($tgt, $dep, @cmd);

#index sites file
$tgt = "$grch38LogDir/$data.sites.$ext.$indexExt.OK";
$dep = "$grch38LogDir/$data.sites.$ext.OK";
@cmd = ("$vt index $grch38OutputDir/$data.sites.$ext 2> $grch38LogDir/$data.sites.$ext.index.log");
makeStep($tgt, $dep, @cmd);

#################################
#trf/lobstr 1000G reference panel
#################################
#$data = "1000g.trf.lobstr";
#$type = "sites";
#
##convert to vcf.gz as picard cannot read bcf files even if piped via STDIN
#$tgt = "$grch38LogDir/$data.grch37.$type.vcf.gz.OK";
#$dep = "$grch37LogDir/$data.$type.$ext.OK";
#@cmd = ("vt view -h $grch37OutputDir/$data.$type.$ext -o $grch38AuxDir/$data.grch37.$type.vcf.gz");
#makeStep($tgt, $dep, @cmd);
#
##liftover
#$tgt = "$grch38AuxDir/$data.$type.vcf.gz.OK";
#$dep = "$refFASTAFile $grch38LogDir/$data.grch37.$type.vcf.gz.OK";
#@cmd = ("$picard LiftoverVcf CHAIN=$hs37ToHs38ChainFile " .
#                           " R=$refFASTAFile " .
#                           " I=$grch38AuxDir/$data.grch37.$type.vcf.gz " .
#                           " O=$grch38AuxDir/$data.$type.vcf.gz " .
#                           " REJECT=$grch38AuxDir/$data.rejected.$type.vcf.gz " .
#                           " 2> $grch38LogDir/$data.$type.liftover.log");
#makeStep($tgt, $dep, @cmd);
#
##convert to ext
#$tgt = "$grch38LogDir/$data.$type.$ext.OK";
#$dep = "$grch38AuxDir/$data.$type.vcf.gz.OK";
#@cmd = ("$vt view $grch38AuxDir/$data.$type.vcf.gz -o $grch38OutputDir/$data.$type.$ext");
#makeStep($tgt, $dep, @cmd);
#
##index sites file
#$tgt = "$grch38LogDir/$data.$type.$ext.$indexExt.OK";
#$dep = "$grch38LogDir/$data.$type.$ext.OK";
#@cmd = ("$vt index $grch38OutputDir/$data.$type.$ext 2> $grch38LogDir/$data.$type.$ext.index.log");
#makeStep($tgt, $dep, @cmd);

#$refBEDFile = "/net/fantasia/home/atks/1000g/20130822_phase3_comparison/20130723_phase3_wg/ALL.wgs.v3_lobSTR.20130502.microsat.integrated_light_weight.genotypes.vcf.gz";
#$data = "1000g.trf.lobstr";
#$destVCFFile = "$data.sites.$ext";

#############################
#trf/VNTRseek reference panel
#############################
$data = "trf.vntrseek";

#convert to vcf.gz as picard cannot read bcf files
$tgt = "$grch38LogDir/$data.grch37.sites.vcf.gz.OK";
$dep = "$grch37LogDir/$data.sites.$ext.OK";
@cmd = ("vt view -h $grch37OutputDir/$data.sites.$ext -o $grch38AuxDir/$data.grch37.sites.vcf.gz");
makeStep($tgt, $dep, @cmd);

#liftover
$tgt = "$grch38AuxDir/$data.sites.vcf.gz.OK";
$dep = "$refFASTAFile $grch38LogDir/$data.grch37.sites.vcf.gz.OK";
@cmd = ("$picard LiftoverVcf CHAIN=$hs37ToHs38ChainFile " .
                           " R=$refFASTAFile " .
                           " I=$grch38AuxDir/$data.grch37.sites.vcf.gz " .
                           " O=$grch38AuxDir/$data.sites.vcf.gz " .
                           " REJECT=$grch38AuxDir/$data.rejected.sites.vcf.gz " .
                           " 2> $grch38LogDir/$data.liftover.log");
makeStep($tgt, $dep, @cmd);

#convert to ext
$tgt = "$grch38LogDir/$data.sites.$ext.OK";
$dep = "$grch38AuxDir/$data.sites.vcf.gz.OK";
@cmd = ("$vt view $grch38AuxDir/$data.sites.vcf.gz -o $grch38OutputDir/$data.sites.$ext");
makeStep($tgt, $dep, @cmd);

#index sites file
$tgt = "$grch38LogDir/$data.sites.$ext.$indexExt.OK";
$dep = "$grch38LogDir/$data.sites.$ext.OK";
@cmd = ("$vt index $grch38OutputDir/$data.sites.$ext 2> $grch38LogDir/$data.sites.$ext.index.log");
makeStep($tgt, $dep, @cmd);

#####################
#trf/vntrseek regions
#####################
$data = "trf.vntrseek";
$srcVCFFile = "$data.sites.$ext";
$destBEDFile = "$data.bed.gz";

#convert from VCF to bed
$tgt = "$grch38LogDir/$destBEDFile.OK";
$dep = "$grch38LogDir/$srcVCFFile.OK";
@cmd = ("$binDir/convert_vcf_2_bed $grch38OutputDir/$srcVCFFile -o $grch38OutputDir/$destBEDFile");
makeStep($tgt, $dep, @cmd);

#index
$tgt = "$grch38LogDir/$destBEDFile.tbi.OK";
$dep = "$grch38LogDir/$destBEDFile.OK";
@cmd = ("tabix -pbed $grch38OutputDir/$destBEDFile");
makeStep($tgt, $dep, @cmd);

################
#Broad OMNI chip
################
$data = "1000G.omni.chip";
$type = "snps.indels.genotypes";

#convert to vcf.gz as picard cannot read bcf files even if piped via STDIN
$tgt = "$grch38LogDir/$data.grch37.$type.vcf.gz.OK";
$dep = "$grch37LogDir/$data.$type.$ext.OK";
@cmd = ("vt view -h $grch37OutputDir/$data.$type.$ext -o $grch38AuxDir/$data.grch37.$type.vcf.gz");
makeStep($tgt, $dep, @cmd);

#liftover
$tgt = "$grch38AuxDir/$data.$type.vcf.gz.OK";
$dep = "$refFASTAFile $grch38LogDir/$data.grch37.$type.vcf.gz.OK";
@cmd = ("$picard LiftoverVcf CHAIN=$hs37ToHs38ChainFile " .
                           " R=$refFASTAFile " .
                           " I=$grch38AuxDir/$data.grch37.$type.vcf.gz " .
                           " O=$grch38AuxDir/$data.$type.vcf.gz " .
                           " REJECT=$grch38AuxDir/$data.rejected.$type.vcf.gz " .
                           " 2> $grch38LogDir/$data.$type.liftover.log");
makeStep($tgt, $dep, @cmd);

#convert to ext
$tgt = "$grch38LogDir/$data.$type.$ext.OK";
$dep = "$grch38AuxDir/$data.$type.vcf.gz.OK";
@cmd = ("$vt view $grch38AuxDir/$data.$type.vcf.gz -o $grch38OutputDir/$data.$type.$ext");
makeStep($tgt, $dep, @cmd);

#index sites file
$tgt = "$grch38LogDir/$data.$type.$ext.$indexExt.OK";
$dep = "$grch38LogDir/$data.$type.$ext.OK";
@cmd = ("$vt index $grch38OutputDir/$data.$type.$ext 2> $grch38LogDir/$data.$type.$ext.index.log");
makeStep($tgt, $dep, @cmd);

####################
#Broad Knowledgebase
####################
$data = "NA12878.broad.kb";
$type = "snps.indels.complex.genotypes";

#convert to vcf.gz as picard cannot read bcf files even if piped via STDIN
$tgt = "$grch38LogDir/$data.grch37.$type.vcf.gz.OK";
$dep = "$grch37LogDir/$data.$type.$ext.OK";
@cmd = ("vt view -h $grch37OutputDir/$data.$type.$ext -o $grch38AuxDir/$data.grch37.$type.vcf.gz");
makeStep($tgt, $dep, @cmd);

#liftover
$tgt = "$grch38AuxDir/$data.$type.vcf.gz.OK";
$dep = "$refFASTAFile $grch38LogDir/$data.grch37.$type.vcf.gz.OK";
@cmd = ("$picard LiftoverVcf CHAIN=$hs37ToHs38ChainFile " .
                           " R=$refFASTAFile " .
                           " I=$grch38AuxDir/$data.grch37.$type.vcf.gz " .
                           " O=$grch38AuxDir/$data.$type.vcf.gz " .
                           " REJECT=$grch38AuxDir/$data.rejected.$type.vcf.gz " .
                           " 2> $grch38LogDir/$data.$type.liftover.log");
makeStep($tgt, $dep, @cmd);

#convert to ext
$tgt = "$grch38LogDir/$data.$type.$ext.OK";
$dep = "$grch38AuxDir/$data.$type.vcf.gz.OK";
@cmd = ("$vt view $grch38AuxDir/$data.$type.vcf.gz -o $grch38OutputDir/$data.$type.$ext");
makeStep($tgt, $dep, @cmd);

#index file
$tgt = "$grch38LogDir/$data.$type.$ext.$indexExt.OK";
$dep = "$grch38LogDir/$data.$type.$ext.OK";
@cmd = ("$vt index $grch38OutputDir/$data.$type.$ext 2> $grch38LogDir/$data.$type.$ext.index.log");
makeStep($tgt, $dep, @cmd);

$type = "snps.indels.complex";

#extract sites
$tgt = "$grch38LogDir/$data.$type.sites.$ext.OK";
$dep = "$grch38LogDir/$data.$type.genotypes.$ext.OK";
@cmd = ("$vt view $grch38OutputDir/$data.$type.genotypes.$ext -o $grch38OutputDir/$data.$type.sites.$ext");
makeStep($tgt, $dep, @cmd);

#index file
$tgt = "$grch38LogDir/$data.$type.sites.$ext.$indexExt.OK";
$dep = "$grch38LogDir/$data.$type.sites.$ext.OK";
@cmd = ("$vt index $grch38OutputDir/$data.$type.sites.$ext 2> $grch38LogDir/$data.$type.sites.$ext.index.log");
makeStep($tgt, $dep, @cmd);

#$srcVCFFile = "/net/fantasia/home/atks/ref/gatk/2.8_b37/NA12878.knowledgebase.snapshot.20131119.b37.vcf.gz";
#$data = "NA12878.broad.kb";
#$destVCFFile = "$data.snps.indels.complex.genotypes.$ext";
#$destSitesVCFFile = "$data.snps.indels.complex.sites.$ext";

#####################
#Illumina Platinum v7
#####################
$data = "NA12878.illumina.platinum.v7";
$type = "snps.indels.complex.genotypes";

#convert to vcf.gz as picard cannot read bcf files even if piped via STDIN
$tgt = "$grch38LogDir/$data.grch37.$type.vcf.gz.OK";
$dep = "$grch37LogDir/$data.$type.$ext.OK";
@cmd = ("vt view -h $grch37OutputDir/$data.$type.$ext -o $grch38AuxDir/$data.grch37.$type.vcf.gz");
makeStep($tgt, $dep, @cmd);

#liftover
$tgt = "$grch38AuxDir/$data.$type.vcf.gz.OK";
$dep = "$refFASTAFile $grch38LogDir/$data.grch37.$type.vcf.gz.OK";
@cmd = ("$picard LiftoverVcf CHAIN=$hs37ToHs38ChainFile " .
                           " R=$refFASTAFile " .
                           " I=$grch38AuxDir/$data.grch37.$type.vcf.gz " .
                           " O=$grch38AuxDir/$data.$type.vcf.gz " .
                           " REJECT=$grch38AuxDir/$data.rejected.$type.vcf.gz " .
                           " 2> $grch38LogDir/$data.$type.liftover.log");
makeStep($tgt, $dep, @cmd);

#convert to ext
$tgt = "$grch38LogDir/$data.$type.$ext.OK";
$dep = "$grch38AuxDir/$data.$type.vcf.gz.OK";
@cmd = ("$vt view $grch38AuxDir/$data.$type.vcf.gz -o $grch38OutputDir/$data.$type.$ext");
makeStep($tgt, $dep, @cmd);

#index sites file
$tgt = "$grch38LogDir/$data.$type.$ext.$indexExt.OK";
$dep = "$grch38LogDir/$data.$type.$ext.OK";
@cmd = ("$vt index $grch38OutputDir/$data.$type.$ext 2> $grch38LogDir/$data.$type.$ext.index.log");
makeStep($tgt, $dep, @cmd);

$type = "snps.indels.complex";

#extract sites
$tgt = "$grch38LogDir/$data.$type.sites.$ext.OK";
$dep = "$grch38LogDir/$data.$type.genotypes.$ext.OK";
@cmd = ("$vt view $grch38OutputDir/$data.$type.genotypes.$ext -o $grch38OutputDir/$data.$type.sites.$ext");
makeStep($tgt, $dep, @cmd);

#index file
$tgt = "$grch38LogDir/$data.$type.sites.$ext.$indexExt.OK";
$dep = "$grch38LogDir/$data.$type.sites.$ext.OK";
@cmd = ("$vt index $grch38OutputDir/$data.$type.sites.$ext 2> $grch38LogDir/$data.$type.sites.$ext.index.log");
makeStep($tgt, $dep, @cmd);

#ConfidentRegions.bed.gz
#$srcVCFFile = "/net/fantasia/home/atks/ref/platinum/v7/NA12878.vcf.gz";
#
#$data = "NA12878.illumina.platinum.v7";
#$destVCFFile = "$data.snps.indels.complex.genotypes.$ext";
#$destSitesVCFFile = "$data.snps.indels.complex.sites.$ext";

#####################
#Illumina Platinum v8
#####################
$data = "NA12878.illumina.platinum.v8";
$type = "snps.indels.complex.genotypes";

#convert to vcf.gz as picard cannot read bcf files even if piped via STDIN
$tgt = "$grch38LogDir/$data.grch37.$type.vcf.gz.OK";
$dep = "$grch37LogDir/$data.$type.$ext.OK";
@cmd = ("vt view -h $grch37OutputDir/$data.$type.$ext -o $grch38AuxDir/$data.grch37.$type.vcf.gz");
makeStep($tgt, $dep, @cmd);

#liftover
$tgt = "$grch38AuxDir/$data.$type.vcf.gz.OK";
$dep = "$refFASTAFile $grch38LogDir/$data.grch37.$type.vcf.gz.OK";
@cmd = ("$picard LiftoverVcf CHAIN=$hs37ToHs38ChainFile " .
                           " R=$refFASTAFile " .
                           " I=$grch38AuxDir/$data.grch37.$type.vcf.gz " .
                           " O=$grch38AuxDir/$data.$type.vcf.gz " .
                           " REJECT=$grch38AuxDir/$data.rejected.$type.vcf.gz " .
                           " 2> $grch38LogDir/$data.$type.liftover.log");
makeStep($tgt, $dep, @cmd);

#convert to ext
$tgt = "$grch38LogDir/$data.$type.$ext.OK";
$dep = "$grch38AuxDir/$data.$type.vcf.gz.OK";
@cmd = ("$vt view $grch38AuxDir/$data.$type.vcf.gz -o $grch38OutputDir/$data.$type.$ext");
makeStep($tgt, $dep, @cmd);

#index sites file
$tgt = "$grch38LogDir/$data.$type.$ext.$indexExt.OK";
$dep = "$grch38LogDir/$data.$type.$ext.OK";
@cmd = ("$vt index $grch38OutputDir/$data.$type.$ext 2> $grch38LogDir/$data.$type.$ext.index.log");
makeStep($tgt, $dep, @cmd);

$type = "snps.indels.complex";

#extract sites
$tgt = "$grch38LogDir/$data.$type.sites.$ext.OK";
$dep = "$grch38LogDir/$data.$type.genotypes.$ext.OK";
@cmd = ("$vt view $grch38OutputDir/$data.$type.genotypes.$ext -o $grch38OutputDir/$data.$type.sites.$ext");
makeStep($tgt, $dep, @cmd);

#index file
$tgt = "$grch38LogDir/$data.$type.sites.$ext.$indexExt.OK";
$dep = "$grch38LogDir/$data.$type.sites.$ext.OK";
@cmd = ("$vt index $grch38OutputDir/$data.$type.sites.$ext 2> $grch38LogDir/$data.$type.sites.$ext.index.log");
makeStep($tgt, $dep, @cmd);

#$srcVCFFile = "/net/fantasia/home/atks/ref/platinum/v8/NA12878.vcf.gz";
#
#$data = "NA12878.illumina.platinum.v8";
#$destVCFFile = "$data.snps.indels.complex.genotypes.$ext";
#$destSitesVCFFile = "$data.snps.indels.complex.sites.$ext";

##############################
#NIST Genome in a Bottle v2.19
##############################
$data = "NA12878.nist.giab.v2.19";
$type = "snps.indels.complex.genotypes";

#convert to vcf.gz as picard cannot read bcf files even if piped via STDIN
$tgt = "$grch38LogDir/$data.grch37.$type.vcf.gz.OK";
$dep = "$grch37LogDir/$data.$type.$ext.OK";
@cmd = ("vt view -h $grch37OutputDir/$data.$type.$ext -o $grch38AuxDir/$data.grch37.$type.vcf.gz");
makeStep($tgt, $dep, @cmd);

#liftover
$tgt = "$grch38AuxDir/$data.$type.vcf.gz.OK";
$dep = "$refFASTAFile $grch38LogDir/$data.grch37.$type.vcf.gz.OK";
@cmd = ("$picard LiftoverVcf CHAIN=$hs37ToHs38ChainFile " .
                           " R=$refFASTAFile " .
                           " I=$grch38AuxDir/$data.grch37.$type.vcf.gz " .
                           " O=$grch38AuxDir/$data.$type.vcf.gz " .
                           " REJECT=$grch38AuxDir/$data.rejected.$type.vcf.gz " .
                           " 2> $grch38LogDir/$data.$type.liftover.log");
makeStep($tgt, $dep, @cmd);

#convert to ext
$tgt = "$grch38LogDir/$data.$type.$ext.OK";
$dep = "$grch38AuxDir/$data.$type.vcf.gz.OK";
@cmd = ("$vt view $grch38AuxDir/$data.$type.vcf.gz -o $grch38OutputDir/$data.$type.$ext");
makeStep($tgt, $dep, @cmd);

#index sites file
$tgt = "$grch38LogDir/$data.$type.$ext.$indexExt.OK";
$dep = "$grch38LogDir/$data.$type.$ext.OK";
@cmd = ("$vt index $grch38OutputDir/$data.$type.$ext 2> $grch38LogDir/$data.$type.$ext.index.log");
makeStep($tgt, $dep, @cmd);

$type = "snps.indels.complex";

#extract sites
$tgt = "$grch38LogDir/$data.$type.sites.$ext.OK";
$dep = "$grch38LogDir/$data.$type.genotypes.$ext.OK";
@cmd = ("$vt view $grch38OutputDir/$data.$type.genotypes.$ext -o $grch38OutputDir/$data.$type.sites.$ext");
makeStep($tgt, $dep, @cmd);

#index file
$tgt = "$grch38LogDir/$data.$type.sites.$ext.$indexExt.OK";
$dep = "$grch38LogDir/$data.$type.sites.$ext.OK";
@cmd = ("$vt index $grch38OutputDir/$data.$type.sites.$ext 2> $grch38LogDir/$data.$type.sites.$ext.index.log");
makeStep($tgt, $dep, @cmd);

#$srcVCFFile = "/net/fantasia/home/atks/ref/giab/v2.19/NISTIntegratedCalls_14datasets_131103_allcall_UGHapMerge_HetHomVarPASS_VQSRv2.19_2mindatasets_5minYesNoRatio_all_nouncert_excludesimplerep_excludesegdups_excludedecoy_excludeRepSeqSTRs_noCNVs.vcf.gz";
#
#$data = "NA12878.nist.giab.v2.19";
#$destVCFFile = "$data.snps.indels.complex.genotypes.$ext";
#$destSitesVCFFile = "$data.snps.indels.complex.sites.$ext";

##################
#DBSNP (from GATK)
##################
$data = "dbSNP138";
$type = "snps.indels.complex.sites";

#convert to vcf.gz as picard cannot read bcf files even if piped via STDIN
$tgt = "$grch38LogDir/$data.grch37.$type.vcf.gz.OK";
$dep = "$grch37LogDir/$data.$type.$ext.OK";
@cmd = ("vt view -h $grch37OutputDir/$data.$type.$ext -o $grch38AuxDir/$data.grch37.$type.vcf.gz");
makeStep($tgt, $dep, @cmd);

#liftover
$tgt = "$grch38AuxDir/$data.$type.vcf.gz.OK";
$dep = "$refFASTAFile $grch38LogDir/$data.grch37.$type.vcf.gz.OK";
@cmd = ("$picard LiftoverVcf CHAIN=$hs37ToHs38ChainFile " .
                           " R=$refFASTAFile " .
                           " I=$grch38AuxDir/$data.grch37.$type.vcf.gz " .
                           " O=$grch38AuxDir/$data.$type.vcf.gz " .
                           " REJECT=$grch38AuxDir/$data.rejected.$type.vcf.gz " .
                           " 2> $grch38LogDir/$data.$type.liftover.log");
makeStep($tgt, $dep, @cmd);

#convert to ext
$tgt = "$grch38LogDir/$data.$type.$ext.OK";
$dep = "$grch38AuxDir/$data.$type.vcf.gz.OK";
@cmd = ("$vt view $grch38AuxDir/$data.$type.vcf.gz -o $grch38OutputDir/$data.$type.$ext");
makeStep($tgt, $dep, @cmd);

#index sites file
$tgt = "$grch38LogDir/$data.$type.$ext.$indexExt.OK";
$dep = "$grch38LogDir/$data.$type.$ext.OK";
@cmd = ("$vt index $grch38OutputDir/$data.$type.$ext 2> $grch38LogDir/$data.$type.$ext.index.log");
makeStep($tgt, $dep, @cmd);

#$srcVCFFile = "/net/fantasia/home/atks/ref/gatk/2.3_b37/dbsnp_138.b37.excluding_sites_after_129.vcf.gz";
#$data = "dbSNP138";
#$destVCFFile = "$data.snps.indels.complex.sites.$ext";

###########
#Mills Chip
###########
$data = "mills.chip";
$type = "indels.genotypes";

#convert to vcf.gz as picard cannot read bcf files even if piped via STDIN
$tgt = "$grch38LogDir/$data.grch37.$type.vcf.gz.OK";
$dep = "$grch37LogDir/$data.$type.$ext.OK";
@cmd = ("vt view -h $grch37OutputDir/$data.$type.$ext -o $grch38AuxDir/$data.grch37.$type.vcf.gz");
makeStep($tgt, $dep, @cmd);

#liftover
$tgt = "$grch38AuxDir/$data.$type.vcf.gz.OK";
$dep = "$refFASTAFile $grch38LogDir/$data.grch37.$type.vcf.gz.OK";
@cmd = ("$picard LiftoverVcf CHAIN=$hs37ToHs38ChainFile " .
                           " R=$refFASTAFile " .
                           " I=$grch38AuxDir/$data.grch37.$type.vcf.gz " .
                           " O=$grch38AuxDir/$data.$type.vcf.gz " .
                           " REJECT=$grch38AuxDir/$data.rejected.$type.vcf.gz " .
                           " 2> $grch38LogDir/$data.$type.liftover.log");
makeStep($tgt, $dep, @cmd);

#convert to ext
$tgt = "$grch38LogDir/$data.$type.$ext.OK";
$dep = "$grch38AuxDir/$data.$type.vcf.gz.OK";
@cmd = ("$vt view $grch38AuxDir/$data.$type.vcf.gz -o $grch38OutputDir/$data.$type.$ext");
makeStep($tgt, $dep, @cmd);

#index sites file
$tgt = "$grch38LogDir/$data.$type.$ext.$indexExt.OK";
$dep = "$grch38LogDir/$data.$type.$ext.OK";
@cmd = ("$vt index $grch38OutputDir/$data.$type.$ext 2> $grch38LogDir/$data.$type.$ext.index.log");
makeStep($tgt, $dep, @cmd);

$type = "indels";

#extract sites
$tgt = "$grch38LogDir/$data.$type.sites.$ext.OK";
$dep = "$grch38LogDir/$data.$type.genotypes.$ext.OK";
@cmd = ("$vt view $grch38OutputDir/$data.$type.genotypes.$ext -o $grch38OutputDir/$data.$type.sites.$ext");
makeStep($tgt, $dep, @cmd);

#index file
$tgt = "$grch38LogDir/$data.$type.sites.$ext.$indexExt.OK";
$dep = "$grch38LogDir/$data.$type.sites.$ext.OK";
@cmd = ("$vt index $grch38OutputDir/$data.$type.sites.$ext 2> $grch38LogDir/$data.$type.sites.$ext.index.log");
makeStep($tgt, $dep, @cmd);

#$srcVCFFile = "/net/fantasia/home/atks/ref/mills/src/mills_indel_genotypes_hg19.vcf";
#$data = "mills.chip";
#$destVCFFile = "$data.indels.genotypes.$ext";
#$destSitesVCFFile = "$data.indels.sites.$ext";

##################
#Mills (from GATK)
##################
$data = "mills";
$type = "indels.sites";

#convert to vcf.gz as picard cannot read bcf files even if piped via STDIN
$tgt = "$grch38LogDir/$data.grch37.$type.vcf.gz.OK";
$dep = "$grch37LogDir/$data.$type.$ext.OK";
@cmd = ("vt view -h $grch37OutputDir/$data.$type.$ext -o $grch38AuxDir/$data.grch37.$type.vcf.gz");
makeStep($tgt, $dep, @cmd);

#liftover
$tgt = "$grch38AuxDir/$data.$type.vcf.gz.OK";
$dep = "$refFASTAFile $grch38LogDir/$data.grch37.$type.vcf.gz.OK";
@cmd = ("$picard LiftoverVcf CHAIN=$hs37ToHs38ChainFile " .
                           " R=$refFASTAFile " .
                           " I=$grch38AuxDir/$data.grch37.$type.vcf.gz " .
                           " O=$grch38AuxDir/$data.$type.vcf.gz " .
                           " REJECT=$grch38AuxDir/$data.rejected.$type.vcf.gz " .
                           " 2> $grch38LogDir/$data.$type.liftover.log");
makeStep($tgt, $dep, @cmd);

#convert to ext
$tgt = "$grch38LogDir/$data.$type.$ext.OK";
$dep = "$grch38AuxDir/$data.$type.vcf.gz.OK";
@cmd = ("$vt view $grch38AuxDir/$data.$type.vcf.gz -o $grch38OutputDir/$data.$type.$ext");
makeStep($tgt, $dep, @cmd);

#index sites file
$tgt = "$grch38LogDir/$data.$type.$ext.$indexExt.OK";
$dep = "$grch38LogDir/$data.$type.$ext.OK";
@cmd = ("$vt index $grch38OutputDir/$data.$type.$ext 2> $grch38LogDir/$data.$type.$ext.index.log");
makeStep($tgt, $dep, @cmd);

#$srcVCFFile = "/net/fantasia/home/atks/ref/gatk/2.3_b37/Mills_and_1000G_gold_standard.indels.b37.vcf.gz";
#$data = "mills";
#$destVCFFile = "$data.indels.sites.$ext";

################
#1000 Genomes v5
################
$data = "1000G.v5";
$type = "snps.indels.complex.sites";

#convert to vcf.gz as picard cannot read bcf files even if piped via STDIN
$tgt = "$grch38LogDir/$data.grch37.$type.vcf.gz.OK";
$dep = "$grch37LogDir/$data.snps.indels.complex.svs.sites.$ext.OK";
@cmd = ("vt view -h $grch37OutputDir/$data.snps.indels.complex.svs.sites.$ext -o $grch38AuxDir/$data.grch37.$type.vcf.gz -f \"VTYPE!=SV\"");
makeStep($tgt, $dep, @cmd);

#lift over
$tgt = "$grch38AuxDir/$data.$type.vcf.gz.OK";
$dep = "$refFASTAFile $grch38LogDir/$data.grch37.$type.vcf.gz.OK";
@cmd = ("$picard LiftoverVcf CHAIN=$hs37ToHs38ChainFile " .
                           " R=$refFASTAFile " .
                           " I=$grch38AuxDir/$data.grch37.$type.vcf.gz " .
                           " O=$grch38AuxDir/$data.$type.vcf.gz " .
                           " REJECT=$grch38AuxDir/$data.rejected.$type.vcf.gz " .
                           " 2> $grch38LogDir/$data.$type.liftover.log");
makeStep($tgt, $dep, @cmd);

#convert to ext
$tgt = "$grch38LogDir/$data.$type.$ext.OK";
$dep = "$grch38AuxDir/$data.$type.vcf.gz.OK";
@cmd = ("$vt view $grch38AuxDir/$data.$type.vcf.gz -o $grch38OutputDir/$data.$type.$ext");
makeStep($tgt, $dep, @cmd);

#index sites file
$tgt = "$grch38LogDir/$data.$type.$ext.$indexExt.OK";
$dep = "$grch38LogDir/$data.$type.$ext.OK";
@cmd = ("$vt index $grch38OutputDir/$data.$type.$ext 2> $grch38LogDir/$data.$type.$ext.index.log");
makeStep($tgt, $dep, @cmd);

#$srcVCFFile = "/net/1000g/1000g/release/20130502/supporting/bcf_files/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.bcf";
#$data = "1000G.v5";
#$destVCFFile = "$data.snps.indels.complex.svs.sites.$ext";

###############
#UK10K 20140722
###############
$data = "UK10K.20140722";
$type = "snps.indels.sites";

#convert to vcf.gz as picard cannot read bcf files even if piped via STDIN
$tgt = "$grch38LogDir/$data.grch37.$type.vcf.gz.OK";
$dep = "$grch37LogDir/$data.$type.$ext.OK";
@cmd = ("vt view -h $grch37OutputDir/$data.$type.$ext -o $grch38AuxDir/$data.grch37.$type.vcf.gz");
makeStep($tgt, $dep, @cmd);

#lift over
$tgt = "$grch38AuxDir/$data.$type.vcf.gz.OK";
$dep = "$refFASTAFile $grch38LogDir/$data.grch37.$type.vcf.gz.OK";
@cmd = ("$picard LiftoverVcf CHAIN=$hs37ToHs38ChainFile " .
                           " R=$refFASTAFile " .
                           " I=$grch38AuxDir/$data.grch37.$type.vcf.gz " .
                           " O=$grch38AuxDir/$data.$type.vcf.gz " .
                           " REJECT=$grch38AuxDir/$data.rejected.$type.vcf.gz " .
                           " 2> $grch38LogDir/$data.$type.liftover.log");
makeStep($tgt, $dep, @cmd);

#convert to ext
$tgt = "$grch38LogDir/$data.$type.$ext.OK";
$dep = "$grch38AuxDir/$data.$type.vcf.gz.OK";
@cmd = ("$vt view $grch38AuxDir/$data.$type.vcf.gz -o $grch38OutputDir/$data.$type.$ext");
makeStep($tgt, $dep, @cmd);

#index sites file
$tgt = "$grch38LogDir/$data.$type.$ext.$indexExt.OK";
$dep = "$grch38LogDir/$data.$type.$ext.OK";
@cmd = ("$vt index $grch38OutputDir/$data.$type.$ext 2> $grch38LogDir/$data.$type.$ext.index.log");
makeStep($tgt, $dep, @cmd);

#$srcVCFFile = "/net/fantasia/home/atks/ref/uk10k/UK10K_COHORT.20140722.sites.vcf.gz";
#$data = "UK10K.20140722";
#$destVCFFile = "$data.snps.indels.sites.$ext";

######################
#Affymetrix exome chip
######################
$data = "affy.exome.chip";
$type = "snps.indels.complex.genotypes";

#convert to vcf.gz as picard cannot read bcf files even if piped via STDIN
$tgt = "$grch38LogDir/$data.grch37.$type.vcf.gz.OK";
$dep = "$grch37LogDir/$data.$type.$ext.OK";
@cmd = ("vt view -h $grch37OutputDir/$data.$type.$ext -o $grch38AuxDir/$data.grch37.$type.vcf.gz");
makeStep($tgt, $dep, @cmd);

#liftover
$tgt = "$grch38AuxDir/$data.$type.vcf.gz.OK";
$dep = "$refFASTAFile $grch38LogDir/$data.grch37.$type.vcf.gz.OK";
@cmd = ("$picard LiftoverVcf CHAIN=$hs37ToHs38ChainFile " .
                           " R=$refFASTAFile " .
                           " I=$grch38AuxDir/$data.grch37.$type.vcf.gz " .
                           " O=$grch38AuxDir/$data.$type.vcf.gz " .
                           " REJECT=$grch38AuxDir/$data.rejected.$type.vcf.gz " .
                           " 2> $grch38LogDir/$data.$type.liftover.log");
makeStep($tgt, $dep, @cmd);

#convert to ext
$tgt = "$grch38LogDir/$data.$type.$ext.OK";
$dep = "$grch38AuxDir/$data.$type.vcf.gz.OK";
@cmd = ("$vt view $grch38AuxDir/$data.$type.vcf.gz -o $grch38OutputDir/$data.$type.$ext");
makeStep($tgt, $dep, @cmd);

#index sites file
$tgt = "$grch38LogDir/$data.$type.$ext.$indexExt.OK";
$dep = "$grch38LogDir/$data.$type.$ext.OK";
@cmd = ("$vt index $grch38OutputDir/$data.$type.$ext 2> $grch38LogDir/$data.$type.$ext.index.log");
makeStep($tgt, $dep, @cmd);

$type = "snps.indels.complex";

#extract sites
$tgt = "$grch38LogDir/$data.$type.sites.$ext.OK";
$dep = "$grch38LogDir/$data.$type.genotypes.$ext.OK";
@cmd = ("$vt view $grch38OutputDir/$data.$type.genotypes.$ext -o $grch38OutputDir/$data.$type.sites.$ext");
makeStep($tgt, $dep, @cmd);

#index file
$tgt = "$grch38LogDir/$data.$type.sites.$ext.$indexExt.OK";
$dep = "$grch38LogDir/$data.$type.sites.$ext.OK";
@cmd = ("$vt index $grch38OutputDir/$data.$type.sites.$ext 2> $grch38LogDir/$data.$type.sites.$ext.index.log");
makeStep($tgt, $dep, @cmd);

#$srcVCFFile = "/net/fantasia/home/atks/ref/1000g/working/20120208_axiom_genotypes/src/ALL.wex.axiom.20120206.snps_and_indels.genotypes.vcf.gz";
#$data = "affy.exome.chip";
#$destVCFFile = "$data.snps.indels.complex.genotypes.$ext";
#$destSitesVCFFile = "$data.snps.indels.complex.sites.$ext";

###################
#CODIS STR data set
###################
#$data = "codis";
#$type = "strs.sites";
#
##convert to vcf.gz as picard cannot read bcf files even if piped via STDIN
#$tgt = "$grch38LogDir/$data.grch37.$type.vcf.gz.OK";
#$dep = "$grch37LogDir/$data.$type.$ext.OK";
#@cmd = ("vt view -h $grch37OutputDir/$data.$type.$ext -o $grch38AuxDir/$data.grch37.$type.vcf.gz");
#makeStep($tgt, $dep, @cmd);
#
##liftover
#$tgt = "$grch38AuxDir/$data.$type.vcf.gz.OK";
#$dep = "$refFASTAFile $grch38LogDir/$data.grch37.$type.vcf.gz.OK";
#@cmd = ("$picard LiftoverVcf CHAIN=$hs37ToHs38ChainFile " .
#                           " R=$refFASTAFile " .
#                           " I=$grch38AuxDir/$data.grch37.$type.vcf.gz " .
#                           " O=$grch38AuxDir/$data.$type.vcf.gz " .
#                           " REJECT=$grch38AuxDir/$data.rejected.$type.vcf.gz " .
#                           " 2> $grch38LogDir/$data.$type.liftover.log");
#makeStep($tgt, $dep, @cmd);
#
##convert to ext
#$tgt = "$grch38LogDir/$data.$type.$ext.OK";
#$dep = "$grch38AuxDir/$data.$type.vcf.gz.OK";
#@cmd = ("$vt view $grch38AuxDir/$data.$type.vcf.gz -o $grch38OutputDir/$data.$type.$ext");
#makeStep($tgt, $dep, @cmd);
#
##index sites file
#$tgt = "$grch38LogDir/$data.$type.$ext.$indexExt.OK";
#$dep = "$grch38LogDir/$data.$type.$ext.OK";
#@cmd = ("$vt index $grch38OutputDir/$data.$type.$ext 2> $grch38LogDir/$data.$type.$ext.index.log");
#makeStep($tgt, $dep, @cmd);

#$srcVCFFile = "/net/fantasia/home/atks/ref/codis/codis.sites.vcf";
#$data = "codis";
#$destVCFFile = "$data.strs.sites.$ext";

#*******************
#Write out make file
#*******************
open(MAK,">$makeFile") || die "Cannot open $makeFile\n";
print MAK ".DELETE_ON_ERROR:\n\n";
print MAK "all: @tgts\n\n";

#clean
push(@tgts,"clean");
$dep = push(@deps, "");
$cmd = ("\t-rm -rf $outputDir/*.bcf $outputDir/*.gz $outputDir/*.bcf.csi $outputDir/*.gz.tbi $logDir");
push(@cmds, $cmd);

for(my $i=0; $i < @tgts; ++$i)
{
    print MAK "$tgts[$i]: $deps[$i]\n";
    print MAK "$cmds[$i]\n";
}

##########
#functions
##########
sub makeMos
{
    my $cmd = shift;
    return ("mosbatch -E/tmp -i -r`/net/fantasia/home/atks/programs/cluster/pick_mini_node` /bin/bash -c 'set pipefail; $cmd'");
}

sub makeStep
{
    my ($tgt, $dep, @cmd) = @_;

    push(@tgts, $tgt);
    push(@deps, $dep);
    my $cmd = "";
    for my $c (@cmd)
    {
        $cmd .= "\t" . $c . "\n";
    }
    $cmd .= "\ttouch $tgt\n";
    push(@cmds, $cmd);
}