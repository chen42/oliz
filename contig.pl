#!/usr/bin/perl

# requirements: cap3 in /opt/cap3/
# contact Xiaoqiu Huang at huang@mtu.edu for cap3

# usage: perl contig.pl DirName
# DirName is the directory where extracted unigene clusters are.
# Contact: Hao Chen, Ph.D. (hchen@utmem.edu)

 
open (CONF, ".oliz") || print "Please read the README file before you proceed!\n";
@conf=<CONF>;
close (CONF);


foreach (@conf) {
	chomp ($_);
	$conf.=$_." ";
}
%conf=split (" ", $conf);

$DirName = $ARGV[0];

if (!defined $DirName) {
	print "Must provide one argument! \n";
	exit;
}

open (ARY, "$DirName.sort");
@array=<ARY>;

print "Number of files to be processed: ",  ($#array +1) ;

for ($i=0; $i < $#array + 1; $i++) {
        #===== move accession number forward so they show up in cap3 output file.
        chop($array[$i]);
        $array[$i]=$conf{species}.".".$array[$i];
        $input = $DirName . "/". $array[$i];
        $output = ">" . $input . ".tr";
        open(UNIN, $input) || die "Can't open $input.";
        open(UNIOUT, $output) || die "Can't open $output.";
        while ($line = <UNIN>) {
                if ($line =~ />/ ) {
                        $Left = index ($line, "gb=");
                        $Right = index ($line, "gi=");
                        $accession = substr ($line, $Left+3, $Right-$Left-5);
                        $line =~ s/>/>$accession /;
                }
                print UNIOUT "$line";
        }
        close UNIN; 
	close UNIOUT;
        #==== run cap3
        $temp1= $array[$i] . ".tr";
        $temp2= $array[$i] . ".out";
        
	
	system ("$conf{\"cap3\"} $DirName/$temp1 > $DirName/$temp2");
	
	
	print ".";
	system ("rm $DirName/$array[$i]");	
	system ("cp $DirName/$temp1 $DirName/$array[$i]");
	
}
system ("rm $DirName/*.tr");
system ("rm $DirName/*.singlets");
system ("rm $DirName/*.links");
system ("rm $DirName/*.ace");
system ("rm $DirName/*.info");
system ("rm $DirName/*.qual");

print "\nStep two done.\nStep three, run: perl utr.pl $DirName\n";

