#!/usr/bin/perl

# usage: perl utr.pl $DirName

#=== INPUT
use Bio::Seq
$DirName = $ARGV[0];

if (!defined $DirName) {
print "Must provide one argument! \n";
exit;
}


open (ARY, "$DirName.sort")|| die "can't locate input file";
@array=<ARY>;

#=== import hand edited unigene, orient unknown accesion number $refseq{$uid};
$input=$DirName . ".ref";
open (REF, $input);
while ($line = <REF>) {
        chop ($line);
        $line =~ s/ //g;
        if ($line ne ""){       
        	$uid=substr($line, 0, index($line, ","));
        	$ref=substr($line, index($line, ",")+1);
        	$refseq{$uid}=$ref;
        	#print "$uid,$ref\n";
        }
}
close (REF);

#=== retrieve Unigene Name from uni.pl output $name{$uid}
$input=$DirName.".name";
open (NAME, $input) || die "$DirName.name";
while ($line=<NAME>) {
	chop ($line);
	if ($line ne ""){	
		$uid=substr($line, 0, index($line, ","));
		$nm=substr($line, index($line, ",")+1);
		$name{$uid}=$nm;
	}
}


print "Generating html files";

#=== file handles
$output = ">". $DirName . ".utr";
open (ARRAYSEQ, $output) || die "arrayseq";
$output = ">". $DirName . ".orient.html";
open (UNKNOWN, $output) || die "orient";
print UNKNOWN "<html><pre>the CDS/Orientation of these contigs are unknown, please edit *.ref file to provide a refseq. \n";
print UNKNOWN "Use refseq accession to indicate 5'-3' and use '-' (without the quotation) to indicate 3'-5'. \n";
print UNKNOWN "No space between Rn number and accession number. One pair per line. For example: \n====\nRn.8174,+\nRn.11929,-\n====\n";
print UNKNOWN "ececute the utr.pl script again incorporates the information into the results\n\n";
open (POLYA, ">$DirName.polya");

opendir (CUR, ".") || die "can't open dir"; 
@files= readdir(CUR);
foreach $file (@files) {
	($skip ="true")	if ($file eq "$DirName.ref")
}



$output = ">". $DirName . ".utr.html" ;
open (UTROUT, $output) ||  die "utr.html";
print UTROUT "<pre>\n";
$output = ">". $DirName. ".nocontig" || die "nocontig";
open (NOCONTIG, $output);
print NOCONTIG "CAP3 did not find a contig for these unigenes. \n";

#===  generates master html file.
$output = ">" . $ARGV[0] . ".html";
open(OUTPUT, $output) or die "Can't open $output.";
print OUTPUT "<html>\n<body link=#000080 vlink=#000080 alink=#FF0000>\n";
print OUTPUT "<center><h2> Data set: $ARGV[0]</h2>\n";
print OUTPUT "<h3>Number of Unigenes in this page: ", ($#array), "</h3><p>\n";


print OUTPUT "<h3><a href= $DirName.utr.html> List of 3'UTR</a></h3>\n";
print OUTPUT "<h3><a href= $DirName.utr> Parsed target 3'UTR sequences </a></h3>\n";

print OUTPUT "<h3><a href= $DirName.orient>Unknow contig orientation</a> </b><br><font size=-2>provided for easier homologue search</font></h3> \n";

print OUTPUT "<h3><a href= $DirName.search_log> Detailed log of specificity screen</a></h3>\n";
print OUTPUT "<h3><a href= $DirName.brief> Brief version of the blast search log </a></h3>\n";
print OUTPUT "<h3><a href= $DirName.final> Final oligo that passed specificity screen </a></h3>\n";






print OUTPUT "</center>Note: 3'UTR is <font color=red>colored</font> in the assembly view<hr>\n";
	
foreach $uid (@array) {
        print ".";
	chop ($uid);
	$uid = "Rn.".$uid;

	#=== continue to print master html file
        print OUTPUT "<a href=http://www.ncbi.nlm.nih.gov/UniGene/clust.cgi?ORG=Rn&CID=$uid&MGC=0>$uid<\/a> $name{$uid}<br>\n";
        system ("cp $DirName/$uid $DirName/$uid.txt");
        print OUTPUT "<a href = ./$DirName/$uid.txt>Sequences </a> |\n";
        print OUTPUT "<a href = ./$DirName/$uid.html>Assembly </a> |\n";
        print OUTPUT "<a href = ./$DirName/$uid.contig.html > Contig </a> |\n";
      
        #===== initiation
	$strand="unknown"; # contig orientation
	$polyA=$utrseq=$contig=$cntglen=$utrlen=$position=$contigseq[1]="";

        #===== formats *.out file into html
	$align = "./$DirName/$uid".".out";
        $outhtml = ">" . "./$DirName/$uid" . ".html";
        open (CAPOUT, $align) || die "Can't open Cap3 output file $align\n";
        open (OUTHTML, $outhtml) || die "Can't open file for output $outhtml\n";
        print OUTHTML "<html>\n<body>\n<font face = courier>\n<pre>\n";
        
	#==== gets theses: $cdend{$acce}, $len{$acce}, $refseq{$uid};
	$input=$DirName."/".$uid;
        open (RN, $input) || die "can't open RN.";
        while ($line = <RN>) {
                if ($line =~ />/) {
                        $Left = index ($line, "gb=");
                        $Right = index ($line, "gi=");
                        $acce = substr ($line, $Left+3, $Right-$Left-5);
                        $len=substr($line, index($line, "len=")+4);
                        $len=substr($len, 0, length($len)-1);
                        $len{$acce}=$len;
                        if (($line =~ /cds=\(/ ) && ($line !~ /cds=UNKNOWN/) ){
                                $Left = index($line, "cds=");
                                $Right= index($line, "gb=");
                                $cdstop= substr($line, $Left+5, $Right-$Left-8);
                                $cdstop= substr($cdstop, index($cdstop, ",")+1);
                                $cdend{$acce}=$cdstop-1;  # the -1 is added to deal with the situ. where the seq. ends at the cds.
                                $refseq{$uid}=$acce unless (defined($refseq{$uid}));
                                $refseq{$uid}=$acce if ($acce =~ /NM_/);
                                #print "$uid,  $refseq{$uid}, $cdend{$acce}\n";
                        }
                }
        }
        close (RN);

        #===== html formating
        while ($line = <CAPOUT>) {
                #======= base counter
                $acc= substr($line, 0, 20); # acc number
                $acc= substr($line, 0, index($acc, " ")-1);
                $forcount=substr($line, 22);
                $forcount=~ s/_| |-//g;
                $forc = substr($forcount, 0, index($forcount, " "));
                $c= length($forc);
                $counter{$acc} +=$c;
                $c=0;
                if ($line =~ /Contig/) {
                        $counter{"consensu"}=0;
                        $utr{"consensu"}="no";
                        $s= substr($line, index($line, "Contig")+7, 1); #number of contig
                }
		
                #==== decide if strand orientation is correct		
		# print "$uid, $acc, $refseq{$uid}\n";
                if ($acc eq $refseq{$uid}) {
                		$position =  $s;
                        if (substr($line, length($acc), 1) eq "-") {$strand="rev.com"; }
                        if (substr($line, length($acc), 1) eq "+") {$strand="correct"; }
                }

		
		#======= color code the end of the CDS   
                if ($strand eq "correct") {
                        if ($utr{$acc} eq "yes") { # full line in color
                        	$str1 = substr($line, 0, 22); $str2= substr($line, 22, length($line)-23);
                        	$line = $str1 . "<font color=red>" . $str2 . "<\/font>\n";
                        } else {
                                if ($cdend{$acc} ne "") {
                                	$start= $counter{$acc} - $cdend{$acc};
					#print "start=$start\n";
                                	if (($start <= 60) && ($start > 0)) {
	                                	$str0 = substr($line, 0, 22);
                                            	$str1 = substr($forcount, 0, 60-$start);
                                            	$str2 = substr($forcount, 60-$start, $start );
                                            	$line = $str0.$str1."<font color=red>".$str2."<\/font>\n";
                                            	$utr{$acc}="yes"; #color subsequent lines
						
                                            	if ($acc eq $refseq{$uid}){
                                            		$conse{$uid}=$counter{"consensu"}+60-$start;
							#print "$cdend{$acc} $conse{$uid} ", $counter{"consensu"}, " $start\n";
                                            	}
                                        }
                                }
                        }
			if ($conse{$uid} > $counter{"consesu"}) {
				$conse {$uid} = $cdend {$refseq{$uid}};
			}
                }
			
	
                #====== color coding the end of CDS, assembly is reversed orientation
                if ($strand eq "rev.com") {
                        if ($cdend{$acc} ne "") {
                                $front= $len{$acc} - $cdend{$acc};                                
                                if ($counter{$acc} < $front) {
                                                $str1 = substr($line, 0, 22); $str2= substr($line, 22, length($line)-23);
                                                $line = $str1 . "<font color=blue>" . $str2 . "<\/font>\n";
                                } else {
					$temp = $counter{$acc}-$front;
                                        if (($temp <= 60)&& ($temp > 0)) {
                                        	if ($counter{$acc}>60) {
         #print "->acc", $acc, "counter", $counter{$acc}, "front", $front, "temp", $temp,"\n";
                                                    	#print $forcount, "\n";
                                                         $str0 = substr($line, 0, 22);
                                                         $str1 = substr($forcount, 0, 60-$temp);
                                                         $str2 = substr($forcount, 60-$temp, $temp+1);
                                                         $line = $str0."<font color=blue>". $str1. "<\/font>" .$str2;
                                                         if ($acc eq $refseq{$uid}) {
                                                              $conse{$uid}=$counter{"consensu"}+60-$temp;
 
                                                         }
                                                 }   
                                        }
                                }
                        }
                }

                #=== print one line of the assembly
                $newline = substr($line, 0, length($line)-1);
                print OUTHTML "$newline";

                #===print base count;
                if ((substr ($newline, 0, 2) ne "  ") && ($counter{$acc} != 0)) {
                                if ($strand eq "rev.com"){
                                		$counterev=$len{$acc}-$counter{$acc};
                                		print OUTHTML "    $counterev of $len{$acc}";
                                } else {
                                        print OUTHTML "    $counter{$acc}";
                                }
                } 
                print OUTHTML "\n";
         
        } # end of while loop
	
        close CAPOUT;
        close OUTHTML;


        #======= generate contig.html
        $output = ">".$DirName ."/" . $uid . ".contig.html";
        open (CONTIGOUT, $output) || die "can't open contig output file";
        
        #===== head of contig.html
        print CONTIGOUT "<html><pre>\n";
        print CONTIGOUT "$uid $name{$uid}";        
        if ($strand eq "unknown") {
                print CONTIGOUT "Contig strand: <font color=red>Unknown </font><p>\n";
        } else {
		print CONTIGOUT "Contig strand: 5' - 3' <p>\n";
	}
  		
	#=== read contig sequences	
	$input = $DirName."/". $uid . ".tr.cap.contigs";
        open (CONTIGIN, $input) || die "can't open contig input file";
	while ($line = <CONTIGIN>) {
                chop ($line);
                if ($line =~ /Contig/) {
                        $cn=substr($line, 7, 1); #contig number
				$contigseq[$cn]="";
		} else {
				$contigseq[$cn] .= $line;
		}
	}
	
	#=== process manual input 
	if ($refseq{$uid} eq "+") {
		$strand="correct";
		$contig=$contigseq[1];
		
		 
	}
	
	if ($refseq{$uid} eq "-") {
		$strand="rev.com";
		$revcom=$contigseq[1];
		$contig=&revcom;
	}
	
	#print $uid, "position, $position, $contigseq[1]\n";
	
	
	
	#=== print seq to contig.html
	for ($x=1; $x<=$cn; $x++) {
		$leng=length ($contigseq[$x]);
		print CONTIGOUT "\n\n>Contig$x,";		
		print CONTIGOUT " length=$leng ";
		 
		if ($x == $position) {
			$contig = $contigseq [$x]; # utr extraction using THE contig contains refseq 
			#====  correct contig orientation if rev.com.
        		if ($strand eq "rev.com") { 
             			$revcom=$contig;
				$contig=&revcom;
			}
			#=== excise trailing seq after polyA tail
			$contig=~ s/a/A/g;
			if ($contig =~ /AAAAAAAAAAAAAAAA/) {
               		 	$A16 = index ($contig, "AAAAAAAAAAAAAAAA");
                		$contig = substr ($contig, 0, $A16+16);
        		}
						
			$cntglen=length($contig);
			print CONTIGOUT "<font color = blue> used for UTR extraction  </font><br>" ;			
		} 
		
		$p=0;
		while ($p < $leng ) {
                	$print = substr ($contigseq[$x], $p, 60);
                	$p+=60;
                	print CONTIGOUT "\n$print  ";
			if (length($print) == 60) {
				print CONTIGOUT $p;
			} else {
				print CONTIGOUT ($p-60+length($print));
			}			
                }        	
        }
	close (CONTIGOUT);
        close (CONTIGIN);

        #======== calculate 3'UTR 
	
	if ($strand eq "rev.com") { 
	 	$conse{$uid}=$cntglen-$conse{$uid}; #position cds ends		 
      	}
     	$conse{$uid}=0 if (!defined $cdend{$refseq{$uid}} );
        $utrseq=substr($contig, $conse{$uid});
	
	
	#=== check polyA tail	
        if (($utrseq =~ /ATTAAA|AATAAA|AATTAA|AATAAT|CATAAA|AGTAAA/) && ($utrseq =~ /AAAAAAAAAA/i)){
                $polyA = "true";
		print POLYA "$uid,polyA\n";
        } else {
		print POLYA "$uid,\n";
	}
	
	#===== guess orientation 
        if ($strand eq "unknown") {
        	if (($polyA) && ($utrseq =~ /AAAAAAAAAA/i)) {
                	$strand= "correct";
        	}
		if (substr($utrseq, 0, 15) =~ /TTTTTTTTTT/i) {
			if (substr($utrseq, 0, 50) =~ /TTTACT|TTTATG|ATTATT|TTAATT|TTTATT|TTTAAT|/){
				$strand="rev.com";
				$revcom=$utrseq;
				$utrseq=&revcom;	
				 
			}
		}
	}		
	
	
	$utrlen=length($utrseq);
	

	#=== print utr seq 		
	
	print UTROUT "<p><b>unigene:</b><a href=http://www.ncbi.nlm.nih.gov/UniGene/clust.cgi?ORG=Rn&CID=$uid&MGC=0>$uid<\/a> $name{$uid}\n";
	print UTROUT "<b>ref_seq</b>:<a href=http://www.ncbi.nlm.nih.gov/htbin-post/Entrez/query?db=n&form=6&uid=$refseq{$uid}&dopt=g>$refseq{$uid}<\/a>";
	
	(print UTROUT "\t") if (length($refseq{$uid})<=8);
	$refutr= $len{$refseq{$uid}}-$cdend{$refseq{$uid}};
	
	print UTROUT " \tlength=$len{$refseq{$uid}} \tcds_end=", (1+$cdend{$refseq{$uid}}), "\t3'utr_length=$refutr<br>";
	print UTROUT "<b>contig</b>: <a href=$DirName/$uid.contig.html>sequence<\/a> \tlength=$cntglen ";
	
	if ($conse{$uid}==0) {
		print UTROUT "<font color = red>CDS unknown, full seq. </font>";
	} else {
        	print UTROUT "\tcds_end=", (1+$conse{$uid});   #####
		print UTROUT "\t3'utr_length=$utrlen ";
	}
print UTROUT "has_polyA " if ($polyA);
	print UTROUT "strand=$strand <br>";
	print UTROUT "<b>assembly</b>:<a href=./$DirName/$uid.html>cap3</a><br>";
	$p=0;
        while ($p < $utrlen) {                        
		$print = substr ($utrseq, $p, 60);
                $p+=60;
                print UTROUT "\n$print  ";
		if (length($print) == 60) {
			print UTROUT $p;
		} else {
			print UTROUT ($p-60+length($print));
		}
        }
	
        #=== finish output html master file.
        if ($polyA) {
        	$Acount++;
                print OUTPUT "polyA: yes | ";
        } else {
	        print OUTPUT "polyA: no  |";
        }
        print OUTPUT " Contig strand: $strand<p><p>\n";
     
        #== Generate seq file for EMBOSS prima. 
        if ($strand ne "unknown") {
                print ARRAYSEQ ">$uid\n";
                if ($polyA eq "") {
                        #print ARRAYSEQ " no polyA sig,";
                        if (($utrlen >= 100) && ($utrlen <= 400)) {
                                #print ARRAYSEQ "3'UTR is $utrlen bp, use 3'UTR";
                        }
			if ($utrlen > 400) {
                            #print ARRAYSEQ " 3'UTR is $utrlen bp, use last 400bp";
                            $utrseq= substr($utrseq, $utrlen-400);
			       
                        }
                        if ($utrlen < 100) {
                            #print ARRAYSEQ " 3'UTR is $utrlen bp, use last 200bp of contig";
                            $utrseq= substr($contig, length($contig)-200); 
			                           
                        }
			if (length($contig)<200) { # last 200bp
			    $utrseq= $contig;
			    
			}
                        
                } else {
                        #print ARRAYSEQ " has polyA sig, ";
                        if (($utrlen >= 100) && ($utrlen<=500) ){
                         	#print ARRAYSEQ " 3'UTR is $utrlen bp, use 3'UTR";
                               
                        }			
			if ($utrlen >500) {
                                #print ARRAYSEQ " 3'UTR is $utrlen bp, use last 500bp";
                                $utrseq= substr($utrseq, $utrlen-500);
				
                                                        }
                        if ($utrlen <100) {
                                #print ARRAYSEQ " 3'UTR 3'UTR is $utrlen bp, use last 500 bp of contig";
                                $utrseq= substr($contig, length($contig)-500);
				
                        }
                        
                }
                
		 
		print ARRAYSEQ "$utrseq\n";
        } else {  # strand="unknown"
		if (length($contig) >0) { 
			if ($skip eq "") {
				$output=">>". $DirName . ".ref";
				open (REFOUT, $output) || die "refout";
				print REFOUT "$uid,\n";
			}
			
			print UNKNOWN "><a href=http://www.ncbi.nlm.nih.gov/HomoloGene/homolquery.cgi?TEXT=$uid>$uid</a> $name{$uid}\n";
                	$p=0;
			$unkn++;
			while ($p < length($contig)) {
                  	      	$print = substr ($contig, $p, 60);
                     	   	$p+=60;
                      	  	print UNKNOWN "$print\n";
			}
			
		} else {
			print NOCONTIG "$uid $name[$i]\n";
			$sngl="true";
		}
        }
         #print "$uid \t$refseq{$uid} \t strand->", $strand, "\n";
}
print OUTPUT "<hr>\n";
print OUTPUT "$Acount unigene have polyA signal <p>\n";
print OUTPUT "</html>\n";

if ($unkn ne "") {
	print " :<{ Warning: the orientation of  ->$unkn<- contigs are unknown, see $DirName.orient file for details.\n";
} else {
	unlink ("$DirName.orient");
}
if ($sngl) {	
	print " :+{ Warning: Some unigene clusters did not produce contig, see $DirName.nocontig for a list.\n";
} else {
	unlink ("$DirName.nocontig");
}
print "\nStep Three done\nStep Four: perl uniq.pl $DirName\n";

sub revcom {
	my ($seqobj, $seq);
	$seqobj = Bio::Seq->new(-seq => $revcom,
		         -id         => $id,
		         -desc       => $desc,
		         );

        $seq=$seqobj->revcom->seq();
}


