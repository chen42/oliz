#!/usr/bin/perl
 
use Bio::Seq;
use Bio::DB::GenBank;

open (CONF, ".oliz") || print "Please read the README file before you proceed!\n";
@conf=<CONF>;
close (CONF);

foreach (@conf) {
	chomp ($_);
	$conf.=$_." ";
}
%conf=split (" ", $conf);
$fastacmd=$conf{blastall};
$fastacmd=~ s/blastall/fastacmd/;

#=== run EMBOSS prima

$DirName=$ARGV[0];
$file=$DirName.".utr";
open (INPUT, $file);
$output=">$DirName.oligo.csv";
open (OLI, "$output");
$prima_in =">prima.seq";

while (<INPUT>) {	 
	#print $l;
	if (/^>/) {		 
		chomp($_);
		$rn=substr($_,1);
		open (PRIMA_IN, "$prima_in");
		print PRIMA_IN $_, "\n";
	} else {
		print PRIMA_IN $_;
		close PRIMA_IN;
		  
		system ("$conf{prima} -sequence prima.seq -targetrange N -outf prima.out -minprimerlen 50  -maxprimerlen 50  -minprimertm 70 -maxprimertm 80 -minprodlen 1 -maxprodlen 1000 ");
		$input="prima.out";
		open (PRIMA_OUT, "$input");
		while (<PRIMA_OUT>) {
			if (/Forward|Reverse/) {
				$oli50=substr($_, index($_,"5'")+3, 50);
				print OLI "$rn,$oli50\n";
			}
		}	
	}
}
 
unlink ("prima.seq");	
unlink ("prima.out");	

#=== creating directories
$dir_b = $DirName."_blast";
open (TEST, "$dir_b");
if(!-e TEST) {
	system ("mkdir $dir_b");
}
close TEST;
$dir_c = $DirName."_clustalw";
open (TEST, "$dir_c");
if(!-e TEST) {
	system ("mkdir $dir_c");
}
$log=$DirName. ".search_log";
$final=$DirName. ".final";
$input= ">". $log;
open (LOG, $input) || die "log file";
$input=">". $final;
open (FINAL, $input) || die "final output";
$input=">". $DirName . ".brief";
open (BRIEF, $input) || die "brief";
print BRIEF "This is a brief version of the search log, listed EST accession is the one\n";
print BRIEF "that disqualified all the oligos in the first round of search, the listed 100% accession\n";
print BRIEF "matches with the oligo 100% but overall similarity with the contig is less than 90%\n";
print BRIEF "UID,REFSEQ,50mer,EST,100%\n";

#========= put the oligos into an array

$file="$DirName.oligo.csv";
open (INPUT, $file);
while ($line =<INPUT>){
        chop ($line);
        $c_uid= substr($line, 0, index($line, ","));
        $oli= substr($line, index($line, ",")+1);
        if ($c_uid eq $uid[$a-1]) {
                $b++;
        } else {
                $b=0;
                $uid[$a]=$c_uid;
                $a++;
        }
        $oligo{$uid[$a-1], $b}=$oli;
}

#===== start the search
for ($i=0; $i<$a; $i++) {
        $j=0;
        print "\n$uid[$i]\n";
        print LOG "\n-->>$uid[$i]\n";

        #=== generate a file contains all accession# from a unigene
        $temp = $DirName. "/". $uid[$i];
        open (ACC, $temp) || die "can't open unigene source $uid[$i]";
        $uni_acc = ">". $temp.".accs";
        open (UNI, $uni_acc) || die "cant' open temp-uni-acc";
        $uni_acc = substr($uni_acc, 1);
        while ($line = <ACC>) {
                if ($line =~ /gb=/) {
                        $Left = index ($line, "gb=");
                        $Right = index ($line, "gi=");
                        $acce = substr ($line, $Left+3, $Right-$Left-5);
                        print UNI "$acce "; $s++;
                        if ($s==6) {
                                $s=0; print UNI "\n";
                        }
                }
        }
        close (UNI);

        #=== process one cluster
        print "\n oligo ";
        while ($oligo{$uid[$i], $j} ne "") {
                print " $j ";
		$countself=0 if ($noest eq "");	
		$satisfied = "yes";
                if ( $oligo{$uid[$i], $j} =~ /gggg/ig) {
                        print LOG "\n  oligo_$j has more than 4 continuous G, discard.";
                        $satisfied = "no";
                        goto LABEL;
                }

                #==== I have to repeat it to make it work during the second no est search, strange
                if ( $oligo{$uid[$i], $j} =~ /gggg/ig) {
                        print LOG "\n>>oligo_$j has more than 4 continuous G, discard.";
                        $satisfied = "no";
                        goto LABEL;
                }

                #======= stand alone blast search
                open (BLASTIN, ">blast.tmp") || die "blastin";
                print BLASTIN $oligo{$uid[$i], $j};
                close (BLASTIN); 
                $boutput = $dir_b. "/" . $uid[$i]. "_oli_". $j. "_blast.html";				       
		if ($noest eq "noest") {
			#print "$conf{a_no_EST_gil_list}\n";
        		system ("$conf{blastall} -p blastn -d $conf{blast_database} -i blast.tmp -T T -F F -l $conf{a_no_EST_gi_list} -o $boutput  >/dev/null");
		} else {
			#print "normal\n";
        		system ("$conf{blastall} -p blastn -d $conf{blast_database} -i blast.tmp -T T -F F  -o $boutput  >/dev/null");
		}

                #======== check blast result
                open (BRESULT, $boutput) || die "blastout" ;
               	$send=$sstart=$qend=$qstart=0;
                $acc=$start="";
                while ($line = <BRESULT>) {
                    if ($line =~ /No hits found/) { #blast failure
				print LOG "oligo $j, no hits found\n";				
				print "\n";			
				goto LABEL;
	            }
		    if ($line =~ /a name/) {
                    	    $revcom="";
                    	    $acc= substr($line,  index($line, "\|")+1);
                    	    $acc= substr($acc, 0, index($acc, "|"));
                    	    if ($acc =~ /\.1/) {
		    		    $acc= substr ($acc, 0, index($acc, ".1")); 
		    	    }
		    	    $start="true";
		    	    
                    }
		    if ($start) {
			if ($line =~ /Identities/) {
                        	$match=substr($line, 14, index($line, "\/")-14);
                        	$total = substr($line, index($line, "\/")+1, index($line, " (")-17);
                        }
                        if ($line =~ /Query/) {
                          	$line=~ s/: |a|g|c|t//ig;
                        	$qstart=substr($line, index($line, "Query")+5, index($line, " ")-5);
                        	$qend=substr($line, rindex($line, " ")+1);
                        	chop($qend);
                        }
                        if ($line =~ /Sbjct/) {
                        	$line=~ s/: |a|g|c|t//ig;
                        	$sstart=substr($line, index($line, "Sbjct")+4, index($line, " ")-3);
                         	$send=substr($line, rindex($line, " ")+1);
                        	chop($send);
                        }
                        $EST="true" if ($line=~ /cDNA clone/);
                        $EST="true" if ($line=~ /EST/);			
                        $fiveprime="true" if ($line=~ /5'/);
                        $threeprime=$revcom="true" if ($line=~ /3'/);			
		
                        if  (($line =~ /<\/PRE>/) && ($acc ne ""))  {
 				#print "acc=$acc, EST=$EST match=$match, total=$total, qry_start=$qstart, qry_end=$qend,  sbj_start=$sstart, sbj_end=$send.\n";
                                #====== check strand
                                if ($EST) {
                                        if (($fiveprime) && ($send > $sstart)) { $strand = "true"; }
                                        if (($threeprime) && ($send < $sstart)) { $strand = "true"; }
                                } else {
                                        if ($send > $sstart) { $strand = "true"; }
                                }
                                $EST=$fiveprime=$threeprime="";
                                
                                if ($strand)  {  # correct strand
                                        $strand="";
                                        #== check if this is from the same unigene set
                                        open (UNI, $uni_acc) || die "can't open temp-uni for read";
                                        $uni="other";
			                while ($line = <UNI>) {
                                               if ($line =~ $acc) {
                                                        $uni="own";
							$countself++ if ($noest eq "");							
                                                        #print "$acc self \n";
                                                        #print LOG "|$acc self ";
                                                        ($refseq{$uid[$i]}=$acc) if ($refseq{$uid[$i]} !~ /NM_/) ;
                                                        last;
                                                }
                                        }
                                        if ($uni eq "other") { #need further check
                                                if ($acc ne "") {
							system ("$fastacmd -d $conf{blast_database} -s $acc >seq.fasta") ;
							open (SEQ, "seq.fasta");
                                                        $sequence = "";
                                                        while ($line = <SEQ>) {
								chop ($line);
								if ($line =~/>/) {
                                                                        $desc = $line;
                                                                } else {
                                                                        $sequence .= $line;
                                                                        $sequence =~ s/ //g;
                                                                }
                                                        }
                                                        close (SEQ);
                                                        $seq = Bio::Seq->new(-seq => $sequence,
		         						-id         => $acc,
		        						-primary_id => $acc,
		         						-desc       => $desc,
								        );
                                                        if (!defined $seq) { # just in case, get seq from NCBI
                                                                & gb_net;								
                                                                $seq=$seq_net;
							}
							$seqlen=$seq->length();
							if ($send < $sstart)  {
                                                		$sstart += ($qstart-1);
                                                		$send = $send - (50-$qend);
                                               		 	$send=1 if ($send<=0);
                                                		$sstart=$seqlen if ($sstart > $seqlen);
                                                		#print "< $acc, send=$send, sstart=$sstart, len=$seqlen:\n";
								$fifty = $seq->trunc($send, $sstart); #fifty is the extended 50mer from the seq found by blast
	                                        		$fifty = $fifty->revcom->seq();
                                        		} else {
                                                		$sstart= $sstart-$qstart+1;
                                                		$send=$send+50-$qend;
                                                		$sstart=1 if ($sstart<=0);
                                                		$send=$seqlen if ($send > $seqlen);						
                                                		$fifty=$seq->trunc($sstart, $send)->seq();
                                        		}

        						#== run clustalw
        						open (CLUSTAL, ">clustaltemp") || die "can't open clustal file";
        						print CLUSTAL ">50mer\n$oligo{$uid[$i], $j}\n>$acc\n$fifty\n";
        						close CLUSTAL;
        						$input= $dir_c. "/". $uid[$i]. ".". $acc. ".aln";
        						system ("/opt/clustal/clustalw /ALGN /INFILE=clustaltemp /GAPOPEN=100 /OUTFILE=$input >/dev/null" );
       						 	open (COUNT, $input) || die "can't open clustal file";
        						while ($line =<COUNT>) {
                						if ($line =~ /\*/) {
                        						$count = ($line =~ tr/\*//);
                        						$stretch=$line;
                						}
        						}
        						$percent = $count/50*100;
							print LOG "\n  oligo_$j: $acc $match/$total, 50mer=$percent% ";
        						if ($percent >75 ) {
                						if ($percent ==100)  {
                        						& check_uni;
                        						#print LOG "\t $uid[$i] does not contain \$acc, need check\n";
                        						$self= "$acc" if ($pct < 90);
                						} else {
                        						$satisfied="no";
                       					 		print LOG "\t->$percent % similar to $acc, discard. ";
                       							last;
                						}
        						} else {
                						if ($percent > 50 ) {			
                    							if ($stretch =~ /\Q***************/) {
                                					$satisfied = "no";
                                					print LOG " >15 cont. bp, discard. ";
                                					last;
                       							}
                						}
        						}
                                                } #if acc=""
	                       		} # if uni
                               } # if strand
                        } # if <pre>
                        if ($line =~ /Number of letters/) {last; }
		} #if start
                } # while loop
		
                if ($satisfied eq "yes") {
                        print BRIEF "$uid[$i],$refseq{$uid[$i]},$oligo{$uid[$i],$j},$disqualify,$self\n";
			print FINAL "$uid[$i],$oligo{$uid[$i],$j}\n";
                        $count_dis++ if ($disqualify ne "");
                        #$count_self++ if ($self ne "");
                        $self=$disqualify="";
                        #print LOG "\n $uid[$i] oligo $j is qualified, it hits $countself sequences in unigene cluster\n";
                        last;
                }
LABEL:          $j++;
	        print LOG "\n";
        } # while loop $j
        if ($satisfied eq "no") {
                if ($noest eq "") {
                        print LOG "\n\n-->blasting again, excluding EST\n";
                        $i=$i-1;
                        $noest = "noest";
                        $disqualify = $acc;
                } else {
                        
			print BRIEF "$uid[$i],$refseq{$uid[$i]},,$acc,$self\n";
                        $disqualify=$self=""; #######
			$count_self++ if ($self ne "");
                        $count_no++;
                        $noest = "";
                        print LOG "\n \tXX: $uid[$i] has no satisfactory oligo \n";
                }
        } else {
                $noest = "";
        }
} # for $i

unlink ("clustaltemp");
unlink ("clustaltemp.dnd");
unlink ("chkuni");
unlink ("chkuni.dnd");
unlink ("blast.tmp");
unlink ("seq.fasta");
close INPUT;
#print LOG "\n statistics:\n total input unigene \t$i;\n no oligo qualify: \t$count_no;\n";
 
print "\n\nJob completed!\nThe qualified oligos are in $DirName.final\nSee $DirName.search_log and $DirName.brief for more details\n\n";
#===
sub gb_net {
	my ($gb);
        print "retrive seq from genbank\n";
        $gb = new Bio::DB::GenBank();
        GENBANK: $seq_net = $gb->get_Seq_by_acc($acc);
        if (! defined ($seq_net)) {
              print "--> failed to retrieve data from genbank, redo\n";
              goto GENBANK;
        }

}

#===== sub
sub check_uni {
        my ($input, $seq, $rc, $x);
        $redo=0;
        system ("$fastacmd -d $conf{blast_database} -s $acc >seq.fasta") ;
                open (SEQ, "seq.fasta");
                $sequence = "";
                while ($line = <SEQ>) {
                        chop ($line);
                        if ($line =~/>/) {
                                $desc = $line;
                        } else {
                                $sequence .= $line;
                                $sequence =~ s/ //g;
                        }
                }
                $seq = Bio::Seq->new(   -seq        => $sequence,
                                        -id         => $acc,
                                        -primary_id => $acc,
                                        -desc       => $desc,
                                        );
        START:
                $pct=$cnt=$lng=$countline=0;
                open (CHKUNI, ">chkuni");
                print CHKUNI  ">$acc\n", $seq->seq(), "\n";
                $input = $DirName . "/$uid[$i]" . ".tr.cap.contigs";
                open (CNTG, $input) || die "sub chek_uni";
                while ($l=<CNTG>) {
                        last if ($l =~ /Contig2/);
                        print CHKUNI $l;
                }
                close CHKUNI;
                $suboutput=$DirName."_clustalw/". $uid[$i]. ".oli_$j.$acc".".self";
                system ("/opt/clustalx1.81/clustalw /ALGN /INFILE=chkuni /GAPOPEN=150 /OUTFILE=$suboutput >/dev/null" );

                open (RESU, $suboutput) || die "chkuni";

                while ($line =<RESU>) {
                        if ($line =~ /\*/) {
                                chop ($line);
                                $l1=substr($line, index($line, "*")) if ($countline ==0);
                                $l2=substr($line, index($line, "*"));
                                $l2=substr($l2, 0, rindex($line, "*")+1);
                                $countline++ ;
                                $cnt += ($line =~ tr/\*//);
                        }
                }
                $lng = length($l1) + 60*($countline-2) + length($l2);
                $pct = $cnt/$lng*100;

                if ($pct > 90) {
                        $input = ">>". $DirName. "/" .$uid[$i] . ".accs";
                        open (APEND, $input) || die "accs";
                        print APEND  " $acc";
                        close (APEND);
                        unlink $suboutput;
                }

                if (($pct < 50) && ($redo == 0)) {
                        $seq  = $seq ->revcom;
                        $redo++;
                        goto START;
                        unlink $suboutput;
                }
                $pct = substr ($pct, 0, 4);
                #print "\n\t >>>Exact match, but $acc is not in $uid[$i], overall $pct% similar to contig\n\n";
                print LOG "\t=>Exact match, but $acc is not in $uid[$i], overall it's $pct% similar to contig";

}

