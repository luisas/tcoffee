
# recompile
cd /home/luisasantus/Desktop/crg_cluster/projects/tcoffee/t_coffee/src/
rm tree_util.o_
g++ -ggdb -w  -MM -I. -I/home/luisasantus/Desktop/crg_cluster/projects/tcoffee/t_coffee/src/../../lib /home/luisasantus/Desktop/crg_cluster/projects/tcoffee/t_coffee/src/../../lib/io_lib/tree_util.h -MF tree_util.d
g++ -ggdb -w  -I. -I/home/luisasantus/Desktop/crg_cluster/projects/tcoffee/t_coffee/src/../../lib -c /home/luisasantus/Desktop/crg_cluster/projects/tcoffee/t_coffee/src/../../lib/io_lib/tree_util.c -o tree_util.o_
g++ -ggdb -w  -lm -o t_coffee util_constraints_list.o_ util_job_handling.o_ util_dps.o_ util_domain_constraints_list.o_ util_analyse_constraints_list.o_ util_aln_analyze.o_ aln_convertion_util.o_ util_declare.o_ hsearch.o_ random.o_ util_make_tree.o_ util.o_ reformat_struc.o_ reformat.o_ aln_compare.o_ io_func.o_ pb_util_read_sequence.o_ pb_util_read_seq_util.o_ tree_util.o_ util_graph_maln.o_ util_dp_clean_maln.o_ util_dp_ssec_pwaln.o_ util_dp_sim.o_ util_dp_mm_nw.o_ util_dp_gotoh_nw.o_ util_dp_suboptimal_nw.o_ util_dp_cdna_fasta_nw.o_ util_dp_generic_fasta_nw.o_ util_dp_fasta_nw.o_ util_dp_fasta_sw.o_ util_dp_gotoh_sw.o_ util_dp_est.o_ util_domain_dp_drivers.o_ util_dp_drivers.o_ util_domain_dp.o_ CUSTOM_evaluate_for_struc.o_ evaluate_for_struc.o_ evaluate_for_domain.o_ evaluate_dirichlet.o_ evaluate.o_ phylo3d.o_ showpair.o_ fsa_dp.o_ pavie_dp.o_ dev1.o_ fastal.o_ parttree.o_ tree.o_ diagonal.o_ fastal_opt_parsing.o_ scoring.o_ iteration.o_ Stack.o_ Vector.o_ classes.o_ km_util.o_ kmeans.o_ km_coffee.o_ t_coffee.o_

/home/luisasantus/Desktop/crg_cluster/projects/tcoffee/t_coffee/src/t_coffee -reg \
                -seq /home/luisasantus/Desktop/crg_cluster/projects/tcoffee/merge_test/data/tiny.fa \
                -reg_nseq 3 \
                -outfile /home/luisasantus/Desktop/crg_cluster/projects/tcoffee/merge_test/data_regressive_standard/tiny_test.aln \
                -outtree /home/luisasantus/Desktop/crg_cluster/projects/tcoffee/merge_test/data_regressive_standard/tiny_test.mbed 2> log.txt


# Command line to generate
#/home/luisasantus/Desktop/crg_cluster/projects/tcoffee/t_coffee/src/t_coffee -reg \
#                -seq /home/luisasantus/Desktop/crg_cluster/projects/tcoffee/merge_test/data/tiny.fa \
#                -reg_nseq 2 \
#                -outfile /home/luisasantus/Desktop/crg_cluster/projects/tcoffee/merge_test/data_regressive_standard/tiny_test.aln \
#                -outtree /home/luisasantus/Desktop/crg_cluster/projects/tcoffee/merge_test/data_regressive_standard/tiny_test.mbed
