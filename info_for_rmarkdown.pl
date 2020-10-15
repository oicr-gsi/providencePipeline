use strict;
use JSON::PP;
use Data::Dumper;
use File::Basename;

my $trimming       = $ARGV[0];
my $samstats       = $ARGV[1];
my $vcf            = $ARGV[2];
my $blast          = $ARGV[3];
my $refseq         = $ARGV[4];
my $insertsizedist = $ARGV[5];
my $orf            = $ARGV[6];
my $mut            = $ARGV[7];
my %info;

#adpater trimming stats
my @trimStats = `cat $trimming`;

foreach my $line (@trimStats) {
    if ( $line =~ /^Input:/ ) {
        my @data     = split( "\t", $line );
        my $rawReads = $data[1];
        $rawReads =~ s/reads//;
        chomp $rawReads;
        $info{"Stats"}{RawReads} = $rawReads;
    }
    if ( $line =~ /^QTrimmed:/ ) {
        my @data                = split( "\t", $line );
        my $qTrimmedReads = $data[1];
        my $qTrimmedBases = $data[2];
        $qTrimmedReads =~ s/reads//;
        $qTrimmedBases =~ s/bases//;
        chomp $qTrimmedReads;
        chomp $qTrimmedBases;
        $info{"Stats"}{QualityTrimmedReads} = $qTrimmedReads;
        $info{"Stats"}{QualityTrimmedBases} = $qTrimmedBases;
    }
    if ( $line =~ /^KTrimmed:/ ) {
        my @data                = split( "\t", $line );
        my $kTrimmedReads = $data[1];
        my $kTrimmedBases = $data[2];
        $kTrimmedReads =~ s/reads//;
        $kTrimmedBases =~ s/bases//;
        chomp $kTrimmedReads;
        chomp $kTrimmedBases;
        $info{"Stats"}{AdapterTrimmedReads} = $kTrimmedReads;
        $info{"Stats"}{AdapterTrimmedBases} = $kTrimmedBases;
    }
		if ( $line =~ /^Trimmed by overlap: / ) {
        my @data                = split( "\t", $line );
        my $tTrimmedReads = $data[1];
        my $tTrimmedBases = $data[2];
        $tTrimmedReads =~ s/reads//;
        $tTrimmedBases =~ s/bases//;
        chomp $tTrimmedReads;
        chomp $tTrimmedBases;
        $info{"Stats"}{TrimmedByOverlapReads} = $tTrimmedReads;
        $info{"Stats"}{TrimmedByOverlapBases} = $tTrimmedBases;
    }
    if ( $line =~ /^Total Removed:/ ) {
        my @data                = split( "\t", $line );
        my $totalTrimmedReads = $data[1];
        my $totalTrimmedBases = $data[2];
        $totalTrimmedReads =~ s/reads//;
        $totalTrimmedBases =~ s/bases//;
        chomp $totalTrimmedReads;
        chomp $totalTrimmedBases;
        $info{"Stats"}{TotalReadsRemoved} = $totalTrimmedReads;
        $info{"Stats"}{TotalBasesRemoved} = $totalTrimmedBases;
    }
}

#samstats

my @samstats = `cat $samstats|grep ^SN`;
my %h;
for my $rec (@samstats) {
    my ( $key, $val ) = $rec =~ /SN\s+(\S.*):\s+([\d\.]+)/;
    $h{$key} = $val;
}

my @trim = (
    "raw total sequences",
    "reads mapped",
    "average length",
    "insert size average",
    "insert size standard deviation"
);
my @sam = map { $h{$_} || 0 } @trim;

$info{"Stats"}{PassedReads} = $sam[0];
$info{"Stats"}{MappedReads} = $sam[1];
my $percent = ( $sam[1] / $sam[0] ) * 100;
my $rounded = sprintf( "%.1f", $percent );
$rounded                       = $rounded . "%";
$info{"Stats"}{MappedPercent}  = $rounded;
$info{"Stats"}{MeanReadLength} = $sam[2];
$info{"Stats"}{MeanInsertSize} = $sam[3];
$info{"Stats"}{InsertSizeSD}    = $sam[4];

#depth and error percent calculation
my @data = `cat $vcf`;
chomp @data;
my $output = 'depth.txt';
open( FH, '>', $output ) or die $!;
print FH "Position\tDepth\tRefCount\tErrorCount\tErrorPercent\n";
my $error_percent_total = 0;
my $count               = 0;

foreach my $line (@data) {
    chomp $line;
    if ( $line !~ /^#/ ) {
        my @info = split( "\t", $line );
        print FH "$info[1]\t";
        if ( $info[7] =~ /DP4\=(.*);MQ.*/ ) {
            my @depth = split( ",", $1 );
            my $ref   = $depth[0] + $depth[1];
            my $alt   = $depth[2] + $depth[3];
            my $dp    = $ref + $alt;
            print FH "$dp\t$ref\t$alt\t";
            my $alt_percent = 0;
            if ( $alt != 0 ) {
                $alt_percent = ( $alt / $dp ) * 100;
            }
            $error_percent_total = $error_percent_total + $alt_percent;
            printf FH ( "%.2f", $alt_percent );
            print FH "\n";
        }
    }
    $count++;
}
close FH;
my $error_mean = $error_percent_total / $count;
my $round      = sprintf( "%.2f", $error_mean );
$info{"ErrorPercent"} = $round;

my @depth = `cat depth.txt`;
open( my $fh2, '>', 'error_gt_1.0.txt' );
print $fh2 "Position\tDepth\tFiltered\tError\tError_Percent\n";

my ( $zero, $point02, $point03, $point04, $point05, $one, $gtone ) = 0;
foreach my $line (@depth) {
    if ( $line !~ /^Pos/ ) {
        my @info  = split( "\t", $line );
        my $error = $info[4];
        if ( $error > 1.00 ) {
            print $fh2 $line;
            $gtone++;
        }
        elsif ( $error == 0 ) {
            $zero++;
        }
        elsif ( ( $error > 0 ) && ( $error <= 0.2 ) ) {
            $point02++;
        }
        elsif ( ( $error > 0.2 ) && ( $error <= 0.3 ) ) {
            $point03++;
        }
        elsif ( ( $error > 0.3 ) && ( $error <= 0.4 ) ) {
            $point04++;
        }
        elsif ( ( $error > 0.4 ) && ( $error <= 0.5 ) ) {
            $point05++;
        }
        elsif ( ( $error > 0.5 ) && ( $error <= 1.0 ) ) {
            $one++;
        }
    }

}
my $all = $zero + $point02 + $point03 + $point04 + $point05 + $one + $gtone;

my $percent_error;
$percent_error            = calc_pert( $all, $zero );
$info{"Error"}{'0'}       = "$zero ($percent_error%)";
$percent_error            = calc_pert( $all, $point02 );
$info{"Error"}{'0-0.2'}   = "$point02 ($percent_error%)";
$percent_error            = calc_pert( $all, $point03 );
$info{"Error"}{'0.2-0.3'} = "$point03 ($percent_error%)";
$percent_error            = calc_pert( $all, $point04 );
$info{"Error"}{"0.3-0.4"} = "$point04 ($percent_error%)";
$percent_error            = calc_pert( $all, $point05 );
$info{"Error"}{"0.4-0.5"} = "$point05 ($percent_error%)";
$percent_error            = calc_pert( $all, $one );
$info{"Error"}{"0.5-1"}   = "$one ($percent_error%)";
$percent_error            = calc_pert( $all, $gtone );
$info{"Error"}{">1"}      = "$gtone ($percent_error%)";

system("ln -s $blast blast.txt");
my @blast = `cat $blast`;
foreach my $line (@blast) {
    if ( $line =~ /Identities.*?(\d+%).*?/ ) {
        my $identity = $1;
        $info{"Consensus"} = $identity;
        last;
    }
}

my @mutdata =`cat $mut |grep -v \#`;
if (@mutdata ==0){
   $info{"Variants"} = "No variants were detected";
   $info{"ORFs"}     = "No changes in ORFs were detected";
}
else {
   $info{"Variants"} = "Variants were detected";
   $info{"ORFs"}     = "Changes in ORFs were detected"; 
}

my $seq = `cat $refseq | grep -v ">"`;

my $Acount  = () = $seq =~ /A/gi;
my $Tcount  = () = $seq =~ /T/gi;
my $Ccount  = () = $seq =~ /C/gi;
my $Gcount  = () = $seq =~ /G/gi;
my $GCcount = $Gcount + $Ccount;
my $total   = $Acount + $Tcount + $Ccount + $Gcount;
$round              = calc_pert( $total, $Acount );
$info{"Seq"}{"A"}   = "$Acount ($round%)";
$round              = calc_pert( $total, $Gcount );
$info{"Seq"}{"G"}   = "$Gcount ($round%)";
$round              = calc_pert( $total, $Ccount );
$info{"Seq"}{"C"}   = "$Ccount ($round%)";
$round              = calc_pert( $total, $Tcount );
$info{"Seq"}{"T"}   = "$Tcount ($round%)";
$info{"TotalBases"} = "$total";

my @insertsizedata = `cat $insertsizedist`;
my $counter        = 0;
foreach my $line (@insertsizedata) {
    if ( $line =~ /^insert_size/ ) {
        my $ioutput = 'insert_size.txt';
        open( IFH, '>', $ioutput ) or die $!;
        print IFH "Insert_Size\tRead_Count\n";
        for ( my $i = $counter + 1 ; $i < @insertsizedata ; $i++ ) {
            print IFH "$insertsizedata[$i]";
        }
    }
    $counter++;
}

my @orfdata = `cat $orf`;
system("ln -s $orf orf.txt");
my $orfout = 'orf_summary.txt';
open( OFH, '>', $orfout ) or die $!;
print OFH "ORF\tStart\tStop\tOrientation\tLength\n";
foreach my $line (@orfdata) {
    chomp $line;
    my $orf;
    if ( $line =~ /^\>/ ) {
        my @info = split( " ", $line );
        if ( $info[0] =~ /.*_(\d+)$/ ) {
            $orf = $1;
        }
        $info[1] =~ s/\[//g;
        $info[3] =~ s/\]//g;
        my $len      = $info[3] - $info[1];
        my $absolute = abs($len);
        my $aa       = ( $absolute + 1 ) / 3;
        print OFH "$orf\t$info[1]\t$info[3]\t";
        if ( $info[4] =~ /REVERSE/ ) {
            print OFH "REVERSE STRAND\t";
        }
        else {
            print OFH "FORWARD STRAND\t";
        }

        print OFH "$aa\n";
    }
}

#print Dumper \%info;
my $json = encode_json( \%info );
print $json;

sub calc_pert {
    my $all   = shift;
    my $count = shift;
    my $round;
    if ( $all != 0 ) {
        my $percent = ( $count / $all ) * 100;
        $round = sprintf( "%.2f", $percent );
    }
    else {
        $round = "0.00";
    }
    return $round;
}

