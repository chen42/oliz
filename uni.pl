#!/usr/bin/perl

# this script extracts UniGene clusters from unigene flat files

# usage: perl uni.pl Input_FileName 

# Input_file is a text file, with one unigene ID 
# (with or without the "Rn") per line.
# Duplicates are OK and are reported in *.missing
# file, together with the Not Found
# unigenes (very likely they are 'retired').

# contact: Hao Chen, Ph.D. hchen@utmem.edu


use Bio::SeqIO; 

# uncomment the following line if the input file is in dos format.

#system ("dos2unix $ARGV[0]");

#creating the configuration file .oliz, which can be manually edited.
print "\n\n===============\nWelcome to Oliz.\n===============\n";

open (CONF, ".oliz") || print "Please read the README file before you proceed!\n";
@conf=<CONF>;
close (CONF);


foreach (@conf) {
	chomp ($_);
	$conf.=$_." ";
}
%conf=split (" ", $conf);

open (CONF, ">>.oliz");


$file="species";
if (! $conf{$file}) {
	print "What is the two letter UniGene code for the species you are working on? e.g. 'Hs' for homo sapien\nYour answer:";
	$conf{"species"}=ucfirst(<STDIN>);
	$conf{"species"}=~s/\.//g;
	chop ($conf{"species"});
	print "\t>>Oliz: I'll remember that. If you want to change that, edit the first line in file .oliz\n";
	print CONF "species $conf{\"species\"}\n";
}
	
$conf{"species"} = $conf{"species"} if  ($conf{"species"} eq "");


$file=$file_name=$conf{"species"}.".seq.all";
&check_file if (! $conf{$file}); 	
$file=$file_name=$conf{"species"}.".data";
&check_file if (! $conf{$file}); 
$file=$file_name="blastall";
&check_file if (! $conf{$file}); 
$file=$file_name="cap3";
&check_file if (! $conf{$file}); 
$file=$file_name="clustalw";
&check_file if (! $conf{$file}); 
$file=$file_name="prima";
&check_file if (! $conf{$file}); 
$file="";
$file_name="blast_database";
$more="Please include the database name\n";
&check_file if (! $conf{$file_name}); 

$file="";
$file_name="a_no_EST_gi_list";
$more="Please include the name of the file. Details on obtaining this list is in README.\n";

&check_file if (! $conf{$file_name}); 


sub check_file {
	my ($found, $dir, $dir_file);
	while (!$found) {
		print "\nHave you obtained $file_name? If not, please press 'CTRL-C' to exit. \nIf so, which directory is it in? $more\nyour answer:";
		
		$dir=<STDIN>;
		chop ($dir);
		 ($dir .= "/") if ($dir!~ /\/$/);
		chop ($dir) if ($file eq "");
		print "$dir, $file\n";
		
		 $dir_file=$dir.$file;
		
		if (-e $dir_file) {
			$found=1;
			print "\t>> Oliz: I'll remember that.\n";
			open (CONF, ">>.oliz");
			print CONF "$file_name $dir_file\n";
			close (CONF);
			$conf{$file_name}=$dir_file;
		} else {
			print "\t>> Oliz: Could not find $dir_file!! please try again.\n";
		}
	}
}



$inputfile= $ARGV[0];

if (length($inputfile) <=6) {
	print "please rename your input file so that it contains more than 6 characters.\n";
	exit;
}


$DirName=substr($inputfile, 0,5);
$DirName=~s/ //g;

open (TEST, "$DirName");
if (!-e TEST) {
	system ("mkdir $DirName");
} else {
	system ("rm $DirName/$conf{\"species\"}.*");
}
open(INPUT, $inputfile) || die "$inputfile.";

$data_=$conf{"species"}.".seq.all";
 
open (DATA, "$conf{$data_}") || die "cant' find $conf{\"species\"}.seq.all";

$notify= ">". $inputfile.".missing";
open (NOTIFY, $notify) || die "can't open file *missing";
$output= ">". $DirName. ".sort";
open (SORT, $output);
$output=">". $DirName. ".singlest";
open (SINGLE, $output);
print SINGLE "The following unigene clusters contain only one EST. Not realiable enough to be included for further analysis\n";

@array0= <INPUT>;

#=== delete "$conf{"species"}.", space, and empty lines
foreach $line (@array0) {
        chop($line);
        $line =~ s/$conf{"species"}|\s|\.//ig;
        if ($line ne "") {
                $array[$i]=$line;
                $i++;
        } else {
		$empt++;
	}
}

#=== delete duplicates and sort
@uniq{@array}=1;
@array= keys %uniq;
@array= sort{$a<=>$b}@array;
$dupli=$#array0-$#array-$empt;

#=====retrieve Unigene Name from $conf{"species"}.data file $name{$uid}
print "\n";
print ">>reading UniGene names... ";

$data_=$conf{"species"}.".data";

open (NAME, "$conf{$data_}") || die "can't open $conf{\"species\"}.data";

$output=">".$DirName.".name";
open (NM, $output);
$j=$currentID=0;
while ($line = <NAME>) { 
	if ($found eq "yes") { $found ="print";}
        if (substr($line, 0, 2) eq "ID") {
                $currentID=substr($line, 15, length($line)-16);
                $found="yes";
                $searchID = $array[$j];
        }	
	
        if ($currentID == $searchID) {
                if ($found eq "print") {
                        $name{$currentID} = substr($line, 12);
                        print NM "$conf{\"species\"}.$searchID,$name{$currentID}";
			$j++; $found = "no";
                }
        }
        if ($currentID > $searchID) {$j++;}
        if ($j == $#array+1) {last;}
}

print "done.\n";



#===  extract unigene clusters
print ">>extracting unigene clusters... \n";
$currentID=$searchID=$prt="";
$last= $array[$#array];
$i= 0;
while ($line = <DATA>) {
        $search = $array [$i];
	
        $searchID = substr($search, 0, length($search));
        if ($line =~ /$conf{"species"}./) {
                $Left = index ($line, "$conf{species}.");
                $Right = index ($line, "len=");
                $currentID = substr ($line, $Left+3, $Right-$Left-5);
        }
        if ($currentID == $searchID) {
                $output = ">>" . $DirName . "/$conf{species}." . $search;
                open(OUTPUT, $output) || die "Can't open $output.";
                $prt = "true";
                if ($line =~ /Set/) {
                        print "\t found $conf{species}.$currentID, last one is $conf{species}.$last\n";
			$found[$i]=$currentID;
                        $prt = "";
			close OUTPUT;
                        $i++;
                }
        }	
	#print "$currentID, $searchID, \n";
        if ($currentID > $searchID ) {
                if ($prt eq "") {
                        print NOTIFY "$conf{species}.$searchID not found\n";
			splice (@array, $i, 1);
			$ntf++;
		}
                
        }
        if ($prt) {
		$line =~ s/cds=UNKNOWN//;
		print OUTPUT "$line";
        }
        if ($i > $#array) {
                close (NOTIFY);
                last;
        }
}


#=== reads $conf{"species"}.unigene seq sets and do ->revcom if clone_end= 3';

print ">>processing EST data";

foreach $uid (@found) {
	
    print ".";
   
    $input= $DirName . "/$conf{species}." . $uid;
   
    # reads $conf{"species"}.unigene seq sets and do ->revcom if the clone_end= 3';

    $in= Bio::SeqIO->new('-file'=> $input ,
                        '-format'=> 'Fasta');
    $output= ">".$input.".tmp";
    open (OUT, $output);
    $num=0; 
    while ($seqobj= $in->next_seq()) {
    	$num++;
        $desc= $seqobj->desc();
        $id= $seqobj->id();
        $regular="true"	if ($desc !~ /clone_end/);
	if ($desc=~ /clone_end=3/) {
                $desc=~ s/3'/5'/g;
                $len= $seqobj->length();
                $tail= $seqobj->seq();
                $tail= substr($tail, ($len-50));
                if ($tail !~ /AAAAAAAAAAAAAAA/) {
                        $seqobj= $seqobj->revcom;
                }
        }
        $sequence= $seqobj->seq();
        $len = $seqobj->length();
        print OUT ">$id $desc\n";
      	$p=0;
	while ($p < $len ) {
               	$print = substr ($sequence, $p, 60);
               	$p +=60;
           	print OUT "$print\n";
        }
    }
    if ($num==1) {
        if ($name{$uid} !~ /EST/) { 	# cap3 only works when multiple seq exit,
		$p=0;  			# this duplicates single entry genes but not ESTs	 					
		print OUT ">$id\n";
		print SORT $uid , "\n";		
		while ($p < $len ) {
               		$print = substr ($sequence, $p, 60);
               		$p+=60;
           		print OUT "$print\n";			
        	}
	} else {
		print SINGLE "$conf{species}.$uid, $name{$uid}\n";
		$sngl++;
	}
    } else {
    	print SORT $uid, "\n";
    }
    unlink ($input);
    system ("cp $input.tmp $input");
    $input=$input.".tmp";
    unlink ($input);
}
 
#=== print warnings 
print "\n";
print " :={ Warning: \n\t$dupli duplicates found in input file.\n" if ($dupli !=0);
if (defined $ntf) {
	print " :={ Warning: \n\t$ntf unigenes were not found (maybe retired). Please see file '$inputfile.missing' for details\n";
} else {
	$notify= $inputfile.".missing";
	unlink ("$notify");
}
if (defined $sngl) {
	print " :={ Warning: \n\t$sngl unigene clusters contain only one EST, see file '$DirName.singlest' for the details.\n";
} else {
	unlink ("$DirName.singlest");
}	
print "Step one done.\nStep two, run: perl contig.pl $DirName\n";
