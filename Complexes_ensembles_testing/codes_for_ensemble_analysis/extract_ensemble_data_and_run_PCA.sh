#!/bin/bash

touch ensemble_RMSDs.txt
while read folder; do
	cp ../$folder/ensemble.txt .
	while read p; do
		struct=${p% *}
		chains=${p#* }
		cp ../$folder/$struct".pdb" .
	done <ensemble.txt

	head -n 1 ensemble.txt > ref.txt
	while read p; do
		ref_pdb=${p% *}
		ref_chains=${p#* }
	done < ref.txt
	
	echo "Analyzing $folder..."
	./write_CA_P $ref_pdb $ref_chains

	grep " CA " $ref_pdb"_"$ref_chains"_CA_P.pdb" > CA.pdb
	num_CA=$(wc -l < CA.pdb)
	rm CA.pdb
	n_CA_P=$(wc -l < $ref_pdb"_"$ref_chains"_CA_P.pdb")

	# Aligning the ensembles using only CA-P and checking structure consistency
	while read p; do
		struct=${p% *}
		chains=${p#* }
		./write_CA_P $struct $chains
		yes 0 | gmx confrms -f1 $ref_pdb"_"$ref_chains"_CA_P.pdb" -f2 $struct"_"$chains"_CA_P.pdb" -one -o temp_file.pdb
		awk 'substr($0,13,4)==" CA " || substr($0,13,4)==" P  "' temp_file.pdb > $struct"_"$chains"_CA_P_ali.pdb"
		rm temp_file.pdb
		n_CA_P_other=$(wc -l < $struct"_"$chains"_CA_P_ali.pdb")
		if test $n_CA_P_other -ne $n_CA_P
		then
			echo $struct"_"$chains "has a different number of CA/P atoms - check it!"
			exit 0
		fi
		./add_chain_label_ali $struct"_"$chains"_CA_P" $struct"_"$chains"_CA_P_ali"
	done <ensemble.txt

	# Compute maximum RMSDs
	max_RMSD=0.0
	grep -F -v "$ref_pdb $ref_chains" ensemble.txt > ensemble_reduced.txt
	num_remaining=$(grep -cve '^[[:space:]]*$' ensemble_reduced.txt)
	pdb_ali=$ref_pdb
	chains_ali=$ref_chains

	while [[ $num_remaining != 0 ]]; do
	    while read p; do
	        struct=${p% *}
	        chains=${p#* }

	        yes 0 | gmx confrms \
	            -f1 "${pdb_ali}_${chains_ali}_CA_P_ali_labeled.pdb" \
	            -f2 "${struct}_${chains}_CA_P_ali_labeled.pdb" > temp.txt

	        RMSD_nm=$(awk '/Root mean square deviation/{print $9}' temp.txt)
	        RMSD_A=$(echo "$RMSD_nm * 10" | bc)

	        rm -f temp.txt fit.pdb

	        if (( $(echo "$RMSD_A > $max_RMSD" | bc -l) )); then
	            max_RMSD=$RMSD_A
	        fi

	    done < <(cat ensemble_reduced.txt)
	    read p < ensemble_reduced.txt
	    pdb_ali=${p% *}
	    chains_ali=${p#* }
	    grep -F -v "$pdb_ali $chains_ali" ensemble_reduced.txt > ensemble_reduced_tmp.txt
	    mv ensemble_reduced_tmp.txt ensemble_reduced.txt
	    num_remaining=$(grep -cve '^[[:space:]]*$' ensemble_reduced.txt)
	done

	# Run PCA
	n_exp_conf=$(wc -l < ensemble.txt)
	echo "Total number of conformations in the ensemble: $n_exp_conf"
	echo " "
	echo "Running PCA ... "
	./run_pca_symm_red $n_CA_P $n_exp_conf "A"
	echo " "
	echo "PCA completed!"
	echo " "
	mv pc_variances.txt ../$folder
	mv PCs.txt ../$folder
	mv exp_ensemble_proj_PC.txt ../$folder

	# Now alignign the ensembles using CA-P-C1-C2 and saving CA-P-C1-C2 aligned PDBs in the folder
	./write_CA_P_C2_C1 $ref_pdb $ref_chains
	gmx make_ndx -f $ref_pdb"_"$ref_chains"_CA_P_C1_C2.pdb" -o index_ref.ndx << EOF
a CA | a P 
q
EOF
	grep "\[" index_ref.ndx > temp.txt
	num_list_ref=$(wc -l <temp.txt)
	rm temp.txt
	echo $num_list_ref
	let num_list_ref=num_list_ref-1
	while read p; do
		struct=${p% *}
		chains=${p#* }
		./write_CA_P_C2_C1 $struct $chains
		gmx make_ndx -f $struct"_"$chains"_CA_P_C1_C2.pdb" -o index_tar.ndx << EOF
a CA | a P 
q
EOF
		grep "\[" index_tar.ndx > temp.txt
		num_list_tar=$(wc -l <temp.txt)
		rm temp.txt
		echo $num_list_tar
		let num_list_tar=num_list_tar-1
		echo "$num_list_ref $num_list_tar" | gmx confrms -f1 $ref_pdb"_"$ref_chains"_CA_P_C1_C2.pdb" -f2 $struct"_"$chains"_CA_P_C1_C2.pdb" -n1 index_ref.ndx -n2 index_tar.ndx -one -o temp_file.pdb
		rm index_tar.ndx
		awk "substr(\$0,13,4)==\" CA \" || substr(\$0,13,4)==\" P  \" || substr(\$0,13,4)==\" C1'\" || substr(\$0,13,4)==\" C2 \"" temp_file.pdb > $struct"_"$chains"_CA_P_C1_C2_ali.pdb"
		rm temp_file.pdb
		./add_chain_label_ali $struct"_"$chains"_CA_P_C1_C2" $struct"_"$chains"_CA_P_C1_C2_ali"
	done <ensemble.txt
	rm index_ref.ndx

	mv *_CA_P_C1_C2_ali_labeled.pdb ../$folder

	rm ensemble.txt ensemble_reduced.txt ref.txt
	rm *.pdb

	echo "$folder , REF: $ref_pdb - $ref_chains, maxRMSD: $max_RMSD, num_CA: $num_CA, num_CA_P: $n_CA_P" >> ensemble_RMSDs.txt
done < folders.txt