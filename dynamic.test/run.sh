# Safe test - should be always working
#../t_coffee/src/t_coffee -reg -seq ./seqs.fa -nseq 2 -reg_tree ./tree.dnd -method famsa_msa
export method_4_CLTCOFFEE=TMalign_pair
export template_file_4_CLTCOFFEE="template_list.txt"
export VERBOSE_4_DYNAMIC=1
#export DUMP_ALN_BUCKETS=1
#export DUMP_SEQ_BUCKETS_ONLY=1

# Only uncomment to recomplile t_coffee
cd ../t_coffee/
rm -r ./src
mkdir src
cp ./src_safe/makefile ./src/
cd src
make t_coffee
cd ../../dynamic.test

../t_coffee/src/t_coffee -reg -reg_method dynamic_msa \
          -seq /home/luisasantus/Desktop/tcoffee/dynamic.test/seqs.fa \
          -reg_tree /home/luisasantus/Desktop/tcoffee/dynamic.test/tree.dnd \
          -reg_nseq 2 \
          -dynamic 1 \
          -dynamic_config /home/luisasantus/Desktop/tcoffee/dynamic.test/dynamic.config \
          -output fasta_aln \
          -reg_homoplasy \
          -thread 2 \
          -outfile /home/luisasantus/Desktop/tcoffee/dynamic.test/seqs_dynamic.aln

###-method TMalign \
