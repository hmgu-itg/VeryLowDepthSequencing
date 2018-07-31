#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Cwd 'abs_path';
use File::Path qw/make_path/;

our $LSF_COMMAND="";
our $JAVA="/software/jre1.7.0_21/bin/java";
our $JAVARGS="-XX:ParallelGCThreads=1 -Xmx6500m";
our $GATK_VERSION="3.1.1";
our $REFERENCE="/lustre/scratch113/teams/zeggini/ReferenceData/ref_human/hs37d5.fa";

our $GATK_CONSTANT_ARGS="--analysis_type VariantRecalibrator -U LENIENT_VCF_PROCESSING --mode SNP -l INFO --use_annotation HaplotypeScore --use_annotation BaseQRankSum --use_annotation MQRankSum --use_annotation ReadPosRankSum --use_annotation FS --ignore_filter LowQual --use_annotation InbreedingCoeff --target_titv 2.15 --TStranche 100.0 --TStranche 99.0 --TStranche 98.0 --TStranche 97.0 --TStranche 96.0 --TStranche 95.0 --TStranche 94.0 --TStranche 93.0 --TStranche 92.0 --TStranche 91.0 --TStranche 90.0 --TStranche 89.0 --TStranche 88.0 --TStranche 87.0 --TStranche 86.0 --TStranche 85.0 --TStranche 84.0 --TStranche 83.0 --TStranche 82.0 --TStranche 81.0 --TStranche 80.0 --TStranche 79.0 --TStranche 78.0 --TStranche 77.0 --TStranche 76.0 --TStranche 75.0 --TStranche 74.0 --TStranche 73.0 --TStranche 72.0 --TStranche 71.0 --TStranche 70.0 --TStranche 69.0 --TStranche 68.0 --TStranche 67.0 --TStranche 66.0 --TStranche 65.0 --TStranche 64.0 --TStranche 63.0 --TStranche 62.0 --TStranche 61.0 --TStranche 60.0 --reference_sequence $REFERENCE -dt NONE -L 11";


our $GATK_JAR=getGATK();
our $DBSNP_PATH=getdbSNPFile();
our $TABIX=getTabix();
our $OUT="/lustre/scratch113/projects/helic/Arthur/VARECH/test";
our @blod_range=(-2.0, -5.0);
our @posg_range=(3, 5, 8, 10, 12, 15);
our @negg_range=(2, 4);
our @numv_range=(1000, 5000, 10000, 12000, 15000, 20000);
our @dbSnpPrior_range=(2, 5);

main();


sub main{
info("Found dbSNP file at $DBSNP_PATH");
err("$ARGV[0] does not exist.") if ( ! -e $ARGV[0]);
info("Found input file at $ARGV[0]");
info("Writing redux file...");
my $dbsnp_redux_path=generateReduxVCF(130, 0, 0);
info("Generating outdir structure...");
generateOutDir();
info("Generating commands...");
my @commands=generateCommands($dbsnp_redux_path);
}


sub generateCommands{
    my ($dbsnp_redux_path)=@_;
    foreach my $blod (@blod_range){
	foreach my $posg (@posg_range){
	    foreach my $negg (@negg_range){
		foreach my $numv (@numv_range){
		    foreach my $dbSnpPrior (@dbSnpPrior_range){
			my $point="$blod.$posg.$negg.$numv.$dbSnpPrior";
			make_path("$OUT/grid/points/$point");
			my $RESOURCE="-resource:hapmap,known=false,training=true,truth=true,prior=15.0 /lustre/scratch113/teams/zeggini/ReferenceData/gatk-bundle/bundle/2.5/b37/hapmap_3.3.b37.vcf -resource:omni,known=false,training=true,truth=false,prior=12.0 /lustre/scratch113/teams/zeggini/ReferenceData/gatk-bundle/bundle/2.5/b37/1000G_omni2.5.b37.vcf -resource:1000g,known=false,training=true,truth=false,prior=10.0 /lustre/scratch113/teams/zeggini/ReferenceData/gatk-bundle/bundle/2.5/b37/1000G_phase1.snps.high_confidence.b37.vcf -resource:dbsnp,known=true,training=false,truth=false,prior=$dbSnpPrior $dbsnp_redux_path --recal_file /dev/null --tranches_file $OUT/grid/points/$point/$point.tranches --rscript_file $OUT/grid/points/$point/$point.R --input $ARGV[0]";
			my $CMD="$LSF_COMMAND $JAVA $JAVARGS -jar $GATK_JAR $GATK_CONSTANT_ARGS --badLodCutoff $blod --maxGaussians $posg --maxNegativeGaussians $negg --minNumBadVariants $numv $RESOURCE";

			print "$CMD\n";
		    }
		}
	    }
	}
    }

}

sub generateOutDir{
use File::Path qw/make_path/;
my $d="$OUT/grid";
if (-e $d) {
    err("the directory $d already exists.");
}
else {
    make_path("$d");
    make_path("$d/commands");
    make_path("$d/points");
    make_path("$d/err");
    make_path("$d/out");
}

}

sub getINFOfieldValue {
    my ($vcfstring, $fieldname)=@_;

    ## Remember to change the info field if different!
    my $INFO_FIELD=7;
    $vcfstring=~/$fieldname=(.*?);/;
    return $1;
    

}

sub getGATK{

    my $path_to_dir=dirname(abs_path($0));
    my $ret="$path_to_dir/gatk/GenomeAnalysisTK-$GATK_VERSION.jar";
    err("Could not find path $ret.") if( ! -e $ret );
    return $ret;

}

sub getTabix{

    my $path_to_dir=dirname(abs_path($0));
    my $ret="$path_to_dir/tabix/tabix";
    err("Could not find path $ret.") if( ! -e $ret );
    return $ret;

}


sub getdbSNPFile{

    my $path_to_dir=dirname(abs_path($0));
    my @files_in_dbsnp_dir=split(/\s/, `ls $path_to_dir/dbSNP`);
    if (scalar(@files_in_dbsnp_dir)>1){
	err("The $path_to_dir has more than one file in its dbSNP folder. This is wrong. That directory should contain only one gzipped VCF of the latest dbSNP release. Please correct and rerun.");
    }else{
	return "$path_to_dir/dbSNP/$files_in_dbsnp_dir[0]";
    }

}

sub generateReduxVCF{
    my ($build, $common, $kg)=@_;
    my $condition="0";    
    if($build){$condition.=' || getINFOfieldValue($_, "dbSNPBuildID")<=$build';}
    if($common){$condition.=' || $_=~/COMMON=1/'}
    if($kg){$condition.=' || $_=~/KGPhase1/'}
    open(IN, "zless $DBSNP_PATH | ") or err("$DBSNP_PATH is not a file.");
    if (! -e dirname(dirname($DBSNP_PATH))."/redux/$build.vcf.gz"){
	open (OUT, "| bgzip > ".dirname(dirname($DBSNP_PATH))."/redux/$build.vcf.gz") or err("Unknown error while trying to write file $build.vcf.gz");
	while(<IN>){
	    print OUT $_ if($_=~/^#/ or eval($condition));
	}
	my $INDEX_CMD="$TABIX ".dirname(dirname($DBSNP_PATH))."/redux/$build.vcf.gz";
	`$INDEX_CMD`;
    }else{ info("File has already been generated in previous run.");}
    return(dirname(dirname($DBSNP_PATH))."/redux/$build.vcf.gz");

}

sub intro{

}

sub info{
    my ($msg)=@_;
    print STDERR "INFO:\t$msg\n";
}

sub err{
    my ($msg)=@_;
    die("ERROR:\t$msg\n");
}


sub test{
 ## Worthless function containing various tests

open(IN, 'zgrep -v \'\#\' /lustre/scratch113/projects/helic/Arthur/VARECH/dbSNP/dbsnp_138.b37.vcf.gz |') or die "fail";
while(<IN>){
    chomp;
    print getINFOfieldValue($_, "dbSNPBuildID"), "\n";

}

}
