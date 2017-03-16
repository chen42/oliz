#/usr/bin/bash
                                             
 if (test $# = 0)                                              
      then echo 'usage: ./onestep.sh genelist_file '            
      exit                                                     
 fi                                                            
  
perl uni.pl $1 
echo '>>assembling contigs..'
perl contig.pl ${1:0:5}
echo '>>parsing UTRs..'
perl utr.pl ${1:0:5}
echo '>>finding uniq oligo sequences..'
perl uniq.pl ${1:0:5}

