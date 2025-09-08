for F in "MYBrelated" "Homeobox" "C2C2dof" "WRKY" "TCP" "bZIP" "MYB" "MADS" "bHLH" "Trihelix" "NAC" "G2like" "HSF" "AP2EREBP" "C2C2gata"                                                                                                               "MYBrelated" "Homeobox" "C2C2dof" "WRKY" "TCP" "bZIP" "MYB" "MADS" "bHLH"                                               "Trihelix" "NAC" "G2like" "HSF" "AP2EREBP" "C2C2gata"
do
  for O in "Athaliana" "Bnapus"
  do
    cd ../results/nemo/${O}/masked_graphpart/${F}_insertion/
    rm -f TF_names.lst
    for FILE in *.tsv
    do
      basename $FILE _medium_exp.tsv >> TF_names.lst
    done
    cd -
  done
done
