#!/usr/bin/perl -w

use warnings;
use strict;
use POSIX;
use Getopt::Long;
use File::Path;
use File::Basename;
use Pod::Usage;

=head1 NAME

generate_vt_bundle_grch37_makefile

=head1 SYNOPSIS

 generate_vt_bundles_grch37_makefile [options]  
 
 --cluster cluster (default main)
  -m       output make file
   
 example: codes/bin/generate_vt_bundles_makefile -m make_ref.mk --chr 20 --cluster 1000g
 
=head1 DESCRIPTION
 
=cut

#option variables
my $help;
my $verbose;
my $debug;
my $makeFile = "create_bundle_grch37.mk";
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
my $assembly = "grch37";
my $homeDir = "/net/fantasia/home/atks";
my $outputDir = "$homeDir/dev/vt/bundle/public/$assembly";
my $binDir = "$homeDir/dev/vt/bundle/bin";
my $logDir = "$homeDir/dev/vt/bundle/log/$assembly";

mkpath($outputDir);
mkpath($logDir);

#programs
my $vt = "/net/fantasia/home/atks/vt/vt";
my $bedtools = "/net/fantasia/home/atks/programs/bedtools2/bin/bedtools";

#reference sequence
my $refFASTAFile = "/net/fantasia/home/atks/ref/genome/hs37d5.fa";

my $indexExt = $ext eq "bcf" ? "csi" : "tbi";

my $srcVCFFile;
my $data;
my $destVCFFile;
my $destSitesVCFFile;
my $srcBEDFile;

###################
#GENCODE Annotation
###################
my $srcGTFFile = "/net/fantasia/home/atks/ref/encode/grch37/gencode.v19.annotation.gtf.gz";
my $destGTFFile = "gencode.v19.annotation.gtf.gz";

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
$srcGTFFile = "/net/fantasia/home/atks/ref/encode/grch37/gencode.v19.annotation.gtf.gz";
my $destBEDFile = "gencode.v19.cds.bed.gz";

#sort and bgzip
$tgt = "$logDir/$destBEDFile.OK";
$dep = "$srcGTFFile";
@cmd = ("zcat $srcGTFFile | grep -v \"^#\" | sort -V -k1,1 -k4,5 | cut -f1,3,4,5 | grep CDS | cut -f1,3,4 | tr '\\t' ':' | uniq | tr ':' '\\t' | bgzip -c > $outputDir/$destBEDFile");
makeStep($tgt, $dep, @cmd);

#index
$tgt = "$logDir/$destBEDFile.tbi.OK";
$dep = "$logDir/$destBEDFile.OK";
@cmd = ("tabix -pbed $outputDir/$destBEDFile");
makeStep($tgt, $dep, @cmd);

##############
#mdust regions
##############
$refFASTAFile = "/net/fantasia/home/atks/ref/genome/hs37d5.fa";
$destBEDFile = "mdust.bed.gz";
my $dustBEDFile = "$outputDir/$destBEDFile";

#sort and bgzip
$tgt = "$logDir/$destBEDFile.OK";
$dep = "$refFASTAFile";
@cmd = ("$binDir/mdust/mdust $refFASTAFile -c | cut -f1,3,4 | bgzip -c > $outputDir/$destBEDFile");
makeStep($tgt, $dep, @cmd);

#index
$tgt = "$logDir/$destBEDFile.tbi.OK";
$dep = "$logDir/$destBEDFile.OK";
@cmd = ("tabix -pbed $outputDir/$destBEDFile");
makeStep($tgt, $dep, @cmd);

#####################
#repeatmasker regions
#####################
my $srcUCSCSQLTableFile = "/net/fantasia/home/atks/ref/ucsc/hg19/rmsk.txt.gz";
$destBEDFile = "rmsk.bed.gz";

#sort and bgzip
$tgt = "$logDir/$destBEDFile.OK";
$dep = "$srcUCSCSQLTableFile";
@cmd = ("zcat $srcUCSCSQLTableFile | cut -f6-8 | " . 
        " perl -lane '{\$\$F[0]=~s/chr//; if (!/_/) {print \"\$\$F[0]\\t\$\$F[1]\\t\$\$F[2]\\n\";}}' | " . 
        " $bedtools sort -i - | " . 
        " sort -V | " .
        " bgzip -c > $outputDir/$destBEDFile");
makeStep($tgt, $dep, @cmd);

#index
$tgt = "$logDir/$destBEDFile.tbi.OK";
$dep = "$logDir/$destBEDFile.OK";
@cmd = ("tabix -pbed $outputDir/$destBEDFile");
makeStep($tgt, $dep, @cmd);

###################
#trf/lobstr regions
###################
my $refBEDFile = "/net/fantasia/home/atks/ref/lobstr/lobstr_v3.0.2_hg19_ref.bed";
$destBEDFile = "trf.lobstr.bed.gz";
my $TRFLobstrBEDFile = "$outputDir/$destBEDFile";

#sort and bgzip
$tgt = "$logDir/$destBEDFile.OK";
$dep = "$refBEDFile";
@cmd = ("$binDir/clean_lobstr_trf_bed $refBEDFile | $bedtools sort -i - | bgzip -c > $outputDir/$destBEDFile");
makeStep($tgt, $dep, @cmd);

#index
$tgt = "$logDir/$destBEDFile.tbi.OK";
$dep = "$logDir/$destBEDFile.OK";
@cmd = ("tabix -pbed $outputDir/$destBEDFile");
makeStep($tgt, $dep, @cmd);

###########################
#trf/lobstr reference panel
###########################
$refBEDFile = "/net/fantasia/home/atks/ref/lobstr/lobstr_v3.0.2_hg19_ref.bed";
$data = "trf.lobstr";
$destVCFFile = "$data.sites.$ext";

#sort and bgzip
$tgt = "$logDir/$destVCFFile.OK";
$dep = "$refFASTAFile $logDir/$destBEDFile.OK";
@cmd = ("$binDir/convert_lobstr_trf_bed_to_vcf $refBEDFile -r $refFASTAFile | $vt sort - | $vt annotate_regions - -b $dustBEDFile -t DUST -d \"Low complexity sequence annotated by mdust\" -o $outputDir/$destVCFFile 2> $logDir/$data.annotate_regions.log");
makeStep($tgt, $dep, @cmd);

#index sites file
$tgt = "$logDir/$destVCFFile.$indexExt.OK";
$dep = "$logDir/$destVCFFile.OK";
@cmd = ("$vt index $outputDir/$destVCFFile 2> $logDir/$data.sites.index.log");
makeStep($tgt, $dep, @cmd);

#################################
#trf/lobstr 1000G reference panel
#################################
#$refBEDFile = "/net/fantasia/home/atks/1000g/20130822_phase3_comparison/20130723_phase3_wg/ALL.wgs.v3_lobSTR.20130502.microsat.integrated_light_weight.genotypes.vcf.gz";
#$data = "1000g.trf.lobstr";
#$destVCFFile = "$data.sites.$ext";
#
##sort and bgzip
#$tgt = "$logDir/$destVCFFile.OK";
#$dep = "$refFASTAFile $logDir/$destBEDFile.OK";
#@cmd = ("$binDir/convert_lobstr_trf_bed_to_vcf $refBEDFile -r $refFASTAFile | $vt sort - | $vt annotate_regions - -b $dustBEDFile -t DUST -d \"Low complexity sequence annotated by mdust\" -o $outputDir/$destVCFFile 2> $logDir/$data.annotate_regions.log");
#makeStep($tgt, $dep, @cmd);

##############################
#trf/repeatseq reference panel
##############################
my $refFile = "/net/fantasia/home/atks/programs/repeatseq/regions/hg19.2014.regions.gz";
$data = "trf.repeatseq";
$destVCFFile = "$data.sites.$ext";

#sort and bgzip
$tgt = "$logDir/$destVCFFile.OK";
$dep = "$refFASTAFile $logDir/$destBEDFile.OK";
@cmd = ("$binDir/convert_repeatseq_resource_2_vcf $refFile | $vt sort - -o +| $vt annotate_regions + -b $dustBEDFile -t DUST -d \"Low complexity sequence annotated by mdust\" -o $outputDir/$destVCFFile 2> $logDir/$data.annotate_regions.log");
makeStep($tgt, $dep, @cmd);

#index sites file
$tgt = "$logDir/$destVCFFile.$indexExt.OK";
$dep = "$logDir/$destVCFFile.OK";
@cmd = ("$vt index $outputDir/$destVCFFile 2> $logDir/$data.sites.index.log");
makeStep($tgt, $dep, @cmd);

#############################
#trf/VNTRseek reference panel
#############################
$refFile = "/net/fantasia/home/atks/programs/VNTRseek-1.08/resource/hg19/t26__.seq";
$data = "trf.vntrseek";
$destVCFFile = "$data.sites.$ext";

#sort and bgzip
$tgt = "$logDir/$destVCFFile.OK";
$dep = "$refFASTAFile $logDir/$destBEDFile.OK";
@cmd = ("$binDir/convert_vntrseek_resource_2_vcf $refFile | $vt sort - -o +| $vt annotate_regions + -b $dustBEDFile -t DUST -d \"Low complexity sequence annotated by mdust\" -o $outputDir/$destVCFFile 2> $logDir/$data.annotate_regions.log");
makeStep($tgt, $dep, @cmd);

#index sites file
$tgt = "$logDir/$destVCFFile.$indexExt.OK";
$dep = "$logDir/$destVCFFile.OK";
@cmd = ("$vt index $outputDir/$destVCFFile 2> $logDir/$data.sites.index.log");
makeStep($tgt, $dep, @cmd);

#make bed file
$srcVCFFile = "$data.sites.$ext";
$destBEDFile = "$data.bed.gz";
$tgt = "$logDir/$destBEDFile.OK";
$dep = "$logDir/$srcVCFFile.OK";
@cmd = ("$binDir/convert_vcf_2_bed $outputDir/$srcVCFFile -o $outputDir/$destBEDFile");
makeStep($tgt, $dep, @cmd);

#index bed file
$srcBEDFile = "$data.bed.gz";
$tgt = "$logDir/$destBEDFile.tbi.OK";
$dep = "$logDir/$srcBEDFile.OK";
@cmd = ("tabix $outputDir/$srcBEDFile");
makeStep($tgt, $dep, @cmd);

################
#Broad OMNI chip
################
$srcVCFFile = "/net/fantasia/home/atks/ref/1000g/working/20131122_broad_omni/Omni25_genotypes_2141_samples.b37.v2.vcf.gz";
$data = "1000G.omni.chip";
$destVCFFile = "$data.snps.indels.genotypes.$ext";
$destSitesVCFFile = "$data.snps.indels.sites.$ext";

#remove unecessary fields, normalize variants and removing duplicates
$tgt = "$logDir/$destVCFFile.OK";
$dep = "$srcVCFFile";
@cmd = ("$binDir/clean_omni_chip $srcVCFFile " . 
      "| $vt normalize - -r $refFASTAFile 2> $logDir/$data.normalize.log " . 
      "| $vt uniq - 2> $logDir/$data.uniq.log " . 
      "| $vt annotate_regions - -b $dustBEDFile -t DUST -d \"Low complexity sequence annotated by mdust\" -o $outputDir/$destVCFFile 2> $logDir/$data.annotate_regions.log");
makeStep($tgt, $dep, @cmd);

#index file
$tgt = "$logDir/$destVCFFile.$indexExt.OK";
$dep = "$logDir/$destVCFFile.OK";
@cmd = ("$vt index $outputDir/$destVCFFile 2> $logDir/$data.index.log");
makeStep($tgt, $dep, @cmd);

#extract sites
$tgt = "$logDir/$destSitesVCFFile.OK";
$dep = "$logDir/$destVCFFile.OK";
@cmd = ("$vt view -h -s -p $outputDir/$destVCFFile -o $outputDir/$destSitesVCFFile 2> $logDir/$data.sites.log");
makeStep($tgt, $dep, @cmd);

#index sites file
$tgt = "$logDir/$destSitesVCFFile.$indexExt.OK";
$dep = "$logDir/$destSitesVCFFile.OK";
@cmd = ("$vt index $outputDir/$destSitesVCFFile 2> $logDir/$data.sites.index.log");
makeStep($tgt, $dep, @cmd);

#summarize file
$tgt = "$logDir/$data.peek.log.OK";
$dep = "$logDir/$destSitesVCFFile.OK";
@cmd = ("$vt peek $outputDir/$destSitesVCFFile 2> $logDir/$data.peek.log");
makeStep($tgt, $dep, @cmd);

####################
#Broad Knowledgebase
####################
$srcVCFFile = "/net/fantasia/home/atks/ref/gatk/2.8_b37/NA12878.knowledgebase.snapshot.20131119.b37.vcf.gz";
$data = "NA12878.broad.kb";
$destVCFFile = "$data.snps.indels.complex.genotypes.$ext";
$destSitesVCFFile = "$data.snps.indels.complex.sites.$ext";

#remove unecessary fields, normalize variants and removing duplicates
$tgt = "$logDir/$destVCFFile.OK";
$dep = "$srcVCFFile";
@cmd = ("zcat $srcVCFFile | $vt normalize - -r $refFASTAFile 2> $logDir/$data.normalize.log | $vt uniq - 2> $logDir/$data.uniq.log | $vt annotate_regions - -b $dustBEDFile -t DUST -d \"Low complexity sequence annotated by mdust\" -o $outputDir/$destVCFFile 2> $logDir/$data.annotate_regions.log");
makeStep($tgt, $dep, @cmd);

##index file
$tgt = "$logDir/$destVCFFile.$indexExt.OK";
$dep = "$logDir/$destVCFFile.OK";
@cmd = ("$vt index $outputDir/$destVCFFile 2> $logDir/$data.index.log");
makeStep($tgt, $dep, @cmd);

#summarize file
$tgt = "$logDir/$data.peek.log.OK";
$dep = "$logDir/$destVCFFile.OK";
@cmd = ("$vt peek $outputDir/$destVCFFile 2> $logDir/$data.peek.log");
makeStep($tgt, $dep, @cmd);

#extract sites
$tgt = "$logDir/$destSitesVCFFile.OK";
$dep = "$logDir/$destVCFFile.OK";
@cmd = ("$vt view -h -s -p $outputDir/$destVCFFile -o $outputDir/$destSitesVCFFile 2> $logDir/$data.sites.log");
makeStep($tgt, $dep, @cmd);

#index sites file
$tgt = "$logDir/$destSitesVCFFile.$indexExt.OK";
$dep = "$logDir/$destSitesVCFFile.OK";
@cmd = ("$vt index $outputDir/$destSitesVCFFile 2> $logDir/$data.sites.index.log");
makeStep($tgt, $dep, @cmd);

#####################
#Illumina Platinum v7
#####################
#ConfidentRegions.bed.gz  
$srcVCFFile = "/net/fantasia/home/atks/ref/platinum/v7/NA12878.vcf.gz";

$data = "NA12878.illumina.platinum.v7";
$destVCFFile = "$data.snps.indels.complex.genotypes.$ext";
$destSitesVCFFile = "$data.snps.indels.complex.sites.$ext";

#remove unecessary fields, normalize variants and removing duplicates
$tgt = "$logDir/$destVCFFile.OK";
$dep = "$srcVCFFile";
@cmd = ("$binDir/clean_illumina_platinum $srcVCFFile " . 
        "| $vt normalize - -r $refFASTAFile -o + 2> $logDir/$data.normalize.log " .
        "| $vt uniq - -o + 2> $logDir/$data.uniq.log " . 
        "| $vt annotate_regions - -b $dustBEDFile -t DUST -d \"Low complexity sequence annotated by mdust\" -o $outputDir/$destVCFFile 2> $logDir/$data.annotate_regions.log");
makeStep($tgt, $dep, @cmd);

##index file
$tgt = "$logDir/$destVCFFile.$indexExt.OK";
$dep = "$logDir/$destVCFFile.OK";
@cmd = ("$vt index $outputDir/$destVCFFile 2> $logDir/$data.index.log");
makeStep($tgt, $dep, @cmd);

#summarize file
$tgt = "$logDir/$data.peek.log.OK";
$dep = "$logDir/$destVCFFile.OK";
@cmd = ("$vt peek $outputDir/$destVCFFile 2> $logDir/$data.peek.log");
makeStep($tgt, $dep, @cmd);

#extract sites
$tgt = "$logDir/$destSitesVCFFile.OK";
$dep = "$logDir/$destVCFFile.OK";
@cmd = ("$vt view -h -s -p $outputDir/$destVCFFile -o $outputDir/$destSitesVCFFile 2> $logDir/$data.sites.log");
makeStep($tgt, $dep, @cmd);

#index sites file
$tgt = "$logDir/$destSitesVCFFile.$indexExt.OK";
$dep = "$logDir/$destSitesVCFFile.OK";
@cmd = ("$vt index $outputDir/$destSitesVCFFile 2> $logDir/$data.sites.index.log");
makeStep($tgt, $dep, @cmd);

#####################
#Illumina Platinum v8
#####################
$srcVCFFile = "/net/fantasia/home/atks/ref/platinum/v8/NA12878.vcf.gz";

$data = "NA12878.illumina.platinum.v8";
$destVCFFile = "$data.snps.indels.complex.genotypes.$ext";
$destSitesVCFFile = "$data.snps.indels.complex.sites.$ext";

#remove unecessary fields, normalize variants and removing duplicates
$tgt = "$logDir/$destVCFFile.OK";
$dep = "$srcVCFFile";
@cmd = ("$binDir/clean_illumina_platinum $srcVCFFile " . 
        "| $vt normalize - -r $refFASTAFile -o + 2> $logDir/$data.normalize.log " .
        "| $vt uniq - -o + 2> $logDir/$data.uniq.log " . 
        "| $vt annotate_regions - -b $dustBEDFile -t DUST -d \"Low complexity sequence annotated by mdust\" -o $outputDir/$destVCFFile 2> $logDir/$data.annotate_regions.log");
makeStep($tgt, $dep, @cmd);

##index file
$tgt = "$logDir/$destVCFFile.$indexExt.OK";
$dep = "$logDir/$destVCFFile.OK";
@cmd = ("$vt index $outputDir/$destVCFFile 2> $logDir/$data.index.log");
makeStep($tgt, $dep, @cmd);

#summarize file
$tgt = "$logDir/$data.peek.log.OK";
$dep = "$logDir/$destVCFFile.OK";
@cmd = ("$vt peek $outputDir/$destVCFFile 2> $logDir/$data.peek.log");
makeStep($tgt, $dep, @cmd);

#extract sites
$tgt = "$logDir/$destSitesVCFFile.OK";
$dep = "$logDir/$destVCFFile.OK";
@cmd = ("$vt view -h -s -p $outputDir/$destVCFFile -o $outputDir/$destSitesVCFFile 2> $logDir/$data.sites.log");
makeStep($tgt, $dep, @cmd);

#index sites file
$tgt = "$logDir/$destSitesVCFFile.$indexExt.OK";
$dep = "$logDir/$destSitesVCFFile.OK";
@cmd = ("$vt index $outputDir/$destSitesVCFFile 2> $logDir/$data.sites.index.log");
makeStep($tgt, $dep, @cmd);

##############################
#NIST Genome in a Bottle v2.19
##############################  
$srcVCFFile = "/net/fantasia/home/atks/ref/giab/v2.19/NISTIntegratedCalls_14datasets_131103_allcall_UGHapMerge_HetHomVarPASS_VQSRv2.19_2mindatasets_5minYesNoRatio_all_nouncert_excludesimplerep_excludesegdups_excludedecoy_excludeRepSeqSTRs_noCNVs.vcf.gz";

$data = "NA12878.nist.giab.v2.19";
$destVCFFile = "$data.snps.indels.complex.genotypes.$ext";
$destSitesVCFFile = "$data.snps.indels.complex.sites.$ext";

#add missing header INFO line, normalize variants and remove duplicates
$tgt = "$logDir/$destVCFFile.OK";
$dep = "$srcVCFFile";
@cmd = ("$binDir/clean_nist_giab_v2.19 $srcVCFFile " .
        "| $vt normalize - -r $refFASTAFile -o + 2> $logDir/$data.normalize.log " .
        "| $vt uniq - -o + 2> $logDir/$data.uniq.log " . 
        "| $vt annotate_regions - -b $dustBEDFile -t DUST -d \"Low complexity sequence annotated by mdust\" -o $outputDir/$destVCFFile 2> $logDir/$data.annotate_regions.log");
makeStep($tgt, $dep, @cmd);

#index file
$tgt = "$logDir/$destVCFFile.$indexExt.OK";
$dep = "$logDir/$destVCFFile.OK";
@cmd = ("$vt index $outputDir/$destVCFFile 2> $logDir/$data.index.log");
makeStep($tgt, $dep, @cmd);

#summarize file
$tgt = "$logDir/$data.peek.log.OK";
$dep = "$logDir/$destVCFFile.OK";
@cmd = ("$vt peek $outputDir/$destVCFFile 2> $logDir/$data.peek.log");
makeStep($tgt, $dep, @cmd);

##################
#DBSNP (from GATK)
##################
$srcVCFFile = "/net/fantasia/home/atks/ref/gatk/2.3_b37/dbsnp_138.b37.excluding_sites_after_129.vcf.gz";
$data = "dbSNP138";
$destVCFFile = "$data.snps.indels.complex.sites.$ext";

#remove unecessary fields, normalize variants and removing duplicates
$tgt = "$logDir/$destVCFFile.OK";
$dep = "$srcVCFFile";
@cmd = ("zcat $srcVCFFile | $vt normalize - -o + -r $refFASTAFile 2> $logDir/$data.normalize.log | $vt uniq + 2> $logDir/$data.uniq.log | $vt annotate_regions - -b $dustBEDFile -t DUST -d \"Low complexity sequence annotated by mdust\" -o $outputDir/$destVCFFile 2> $logDir/$data.annotate_regions.log");
makeStep($tgt, $dep, @cmd);

#index file
$tgt = "$logDir/$destVCFFile.$indexExt.OK";
$dep = "$logDir/$destVCFFile.OK";
@cmd = ("$vt index $outputDir/$destVCFFile 2> $logDir/$data.index.log");
makeStep($tgt, $dep, @cmd);

#summarize file
$tgt = "$logDir/$data.peek.log.OK";
$dep = "$logDir/$destVCFFile.OK";
@cmd = ("$vt peek $outputDir/$destVCFFile 2> $logDir/$data.peek.log");
makeStep($tgt, $dep, @cmd);

###########
#Mills Chip
###########
$srcVCFFile = "/net/fantasia/home/atks/ref/mills/src/mills_indel_genotypes_hg19.vcf";
$data = "mills.chip";
$destVCFFile = "$data.indels.genotypes.$ext";
$destSitesVCFFile = "$data.indels.sites.$ext";

#remove unecessary fields, normalize variants and removing duplicates
$tgt = "$logDir/$destVCFFile.OK";
$dep = "$srcVCFFile";
@cmd = ("$binDir/clean_mills_chip $srcVCFFile | $vt normalize - -r $refFASTAFile 2> $logDir/$data.normalize.log | $vt uniq - 2> $logDir/$data.uniq.log | $vt annotate_regions - -b $dustBEDFile -t DUST -d \"Low complexity sequence annotated by mdust\" -o $outputDir/$destVCFFile 2> $logDir/$data.annotate_regions.log");
makeStep($tgt, $dep, @cmd);

#index file
$tgt = "$logDir/$destVCFFile.$indexExt.OK";
$dep = "$logDir/$destVCFFile.OK";
@cmd = ("$vt index $outputDir/$destVCFFile 2> $logDir/$data.index.log");
makeStep($tgt, $dep, @cmd);

#summarize file
$tgt = "$logDir/$data.peek.log.OK";
$dep = "$logDir/$destVCFFile.OK";
@cmd = ("$vt peek $outputDir/$destVCFFile 2> $logDir/$data.peek.log");
makeStep($tgt, $dep, @cmd);

#extract sites
$tgt = "$logDir/$destSitesVCFFile.OK";
$dep = "$logDir/$destVCFFile.OK";
@cmd = ("$vt view -h -s -p $outputDir/$destVCFFile -o $outputDir/$destSitesVCFFile 2> $logDir/$data.sites.log");
makeStep($tgt, $dep, @cmd);

#index sites file
$tgt = "$logDir/$destSitesVCFFile.$indexExt.OK";
$dep = "$logDir/$destSitesVCFFile.OK";
@cmd = ("$vt index $outputDir/$destSitesVCFFile 2> $logDir/$data.sites.index.log");
makeStep($tgt, $dep, @cmd);

##################
#Mills (from GATK)
##################
$srcVCFFile = "/net/fantasia/home/atks/ref/gatk/2.3_b37/Mills_and_1000G_gold_standard.indels.b37.vcf.gz";
$data = "mills";
$destVCFFile = "$data.indels.sites.$ext";

#remove unecessary fields, normalize variants and removing duplicates
$tgt = "$logDir/$destVCFFile.OK";
$dep = "$srcVCFFile";
@cmd = ("$binDir/clean_mills_1000g_gatk $srcVCFFile | $vt normalize - -r $refFASTAFile 2> $logDir/$data.normalize.log | $vt uniq - 2> $logDir/$data.uniq.log | $vt annotate_regions - -b $dustBEDFile -t DUST -d \"Low complexity sequence annotated by mdust\" -o $outputDir/$destVCFFile 2> $logDir/$data.annotate_regions.log");
makeStep($tgt, $dep, @cmd);

#index file
$tgt = "$logDir/$destVCFFile.$indexExt.OK";
$dep = "$logDir/$destVCFFile.OK";
@cmd = ("$vt index $outputDir/$destVCFFile 2> $logDir/$data.index.log");
makeStep($tgt, $dep, @cmd);

#summarize file
$tgt = "$logDir/$data.peek.log.OK";
$dep = "$logDir/$destVCFFile.OK";
@cmd = ("$vt peek $outputDir/$destVCFFile 2> $logDir/$data.peek.log");
makeStep($tgt, $dep, @cmd);

################
#1000 Genomes v5
################
$srcVCFFile = "/net/1000g/1000g/release/20130502/supporting/bcf_files/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.bcf";
$data = "1000G.v5";
$destVCFFile = "$data.snps.indels.complex.svs.sites.$ext";

#remove unecessary fields, normalize variants and removing duplicates
$tgt = "$logDir/$destVCFFile.OK";
$dep = "$srcVCFFile";
@cmd = ("$vt view -h $srcVCFFile -s | $vt normalize - -r $refFASTAFile 2> $logDir/$data.normalize.log | $vt uniq - 2> $logDir/$data.uniq.log | $vt annotate_regions - -b $dustBEDFile -t DUST -d \"Low complexity sequence annotated by mdust\" -o $outputDir/$destVCFFile 2> $logDir/$data.annotate_regions.log");
makeStep($tgt,$dep,@cmd);

#index file
$tgt = "$logDir/$destVCFFile.$indexExt.OK";
$dep = "$logDir/$destVCFFile.OK";
@cmd = ("$vt index $outputDir/$destVCFFile 2> $logDir/$data.index.log");
makeStep($tgt,$dep,@cmd);

#summarize file
$tgt = "$logDir/$data.peek.log.OK";
$dep = "$logDir/$destVCFFile.OK";
@cmd = ("$vt peek $outputDir/$destVCFFile 2> $logDir/$data.peek.log");
makeStep($tgt,$dep,@cmd);

###############
#UK10K 20140722
###############
$srcVCFFile = "/net/fantasia/home/atks/ref/uk10k/UK10K_COHORT.20140722.sites.vcf.gz";
$data = "UK10K.20140722";
$destVCFFile = "$data.snps.indels.sites.$ext";

#remove unecessary fields, normalize variants and removing duplicates
$tgt = "$logDir/$destVCFFile.OK";
$dep = "$srcVCFFile";
@cmd = ("$vt view -h $srcVCFFile -s | $vt normalize - -r $refFASTAFile 2> $logDir/$data.normalize.log | $vt uniq - 2> $logDir/$data.uniq.log | $vt annotate_regions - -b $dustBEDFile -t DUST -d \"Low complexity sequence annotated by mdust\" -o $outputDir/$destVCFFile 2> $logDir/$data.annotate_regions.log");
makeStep($tgt,$dep,@cmd);

#index file
$tgt = "$logDir/$destVCFFile.$indexExt.OK";
$dep = "$logDir/$destVCFFile.OK";
@cmd = ("$vt index $outputDir/$destVCFFile 2> $logDir/$data.index.log");
makeStep($tgt,$dep,@cmd);

#summarize file
$tgt = "$logDir/$data.peek.log.OK";
$dep = "$logDir/$destVCFFile.OK";
@cmd = ("$vt peek $outputDir/$destVCFFile 2> $logDir/$data.peek.log");
makeStep($tgt,$dep,@cmd);

######################
#Affymetrix exome chip
#####################
$srcVCFFile = "/net/fantasia/home/atks/ref/1000g/working/20120208_axiom_genotypes/src/ALL.wex.axiom.20120206.snps_and_indels.genotypes.vcf.gz";
$data = "affy.exome.chip";
$destVCFFile = "$data.snps.indels.complex.genotypes.$ext";
$destSitesVCFFile = "$data.snps.indels.complex.sites.$ext";

#remove unecessary fields, normalize variants and removing duplicates
$tgt = "$logDir/$destVCFFile.OK";
$dep = "$srcVCFFile";
@cmd = ("$binDir/clean_affy_exome_chip $srcVCFFile | $vt sort - | $vt normalize - -r $refFASTAFile 2> $logDir/$data.normalize.log | $vt uniq - 2> $logDir/$data.uniq.log | $vt annotate_regions - -b $dustBEDFile -t DUST -d \"Low complexity sequence annotated by mdust\" -o $outputDir/$destVCFFile 2> $logDir/$data.annotate_regions.log");
makeStep($tgt, $dep, @cmd);

#index file
$tgt = "$logDir/$destVCFFile.$indexExt.OK";
$dep = "$logDir/$destVCFFile.OK";
@cmd = ("$vt index $outputDir/$destVCFFile 2> $logDir/$data.index.log");
makeStep($tgt, $dep, @cmd);

#summarize file
$tgt = "$logDir/$data.peek.log.OK";
$dep = "$logDir/$destVCFFile.OK";
@cmd = ("$vt peek $outputDir/$destVCFFile 2> $logDir/$data.peek.log");
makeStep($tgt, $dep, @cmd);

#extract sites
$tgt = "$logDir/$destSitesVCFFile.OK";
$dep = "$logDir/$destVCFFile.OK";
@cmd = ("$vt view -h -s -p $outputDir/$destVCFFile -o $outputDir/$destSitesVCFFile 2> $logDir/$data.sites.log");
makeStep($tgt, $dep, @cmd);

#index sites file
$tgt = "$logDir/$destSitesVCFFile.$indexExt.OK";
$dep = "$logDir/$destSitesVCFFile.OK";
@cmd = ("$vt index $outputDir/$destSitesVCFFile 2> $logDir/$data.sites.index.log");
makeStep($tgt, $dep, @cmd);

###################
#CODIS STR data set
###################
$srcVCFFile = "/net/fantasia/home/atks/ref/codis/codis.sites.vcf";
$data = "codis";
$destVCFFile = "$data.strs.sites.$ext";

#annotate variants
$tgt = "$logDir/$destVCFFile.OK";
$dep = "$srcVCFFile";
@cmd = ("$vt annotate_regions $srcVCFFile -b $dustBEDFile -t DUST -d \"Low complexity sequence annotated by mdust\" -o $outputDir/$destVCFFile 2> $logDir/$data.annotate_regions.log");
makeStep($tgt, $dep, @cmd);

#index sites file
$tgt = "$logDir/$destVCFFile.$indexExt.OK";
$dep = "$logDir/$destVCFFile.OK";
@cmd = ("$vt index $outputDir/$destVCFFile 2> $logDir/$data.sites.index.log");
makeStep($tgt, $dep, @cmd);

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