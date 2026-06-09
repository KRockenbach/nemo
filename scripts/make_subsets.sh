# exclude .npz, .png, .h5, .html
for file in $(find ../results/* | grep -v "\.npz" | grep -v "\.png" | grep ".*\....$")
do
  new=$(echo $file | sed 's/results/result_subset/g')
  mkdir -p $(dirname $new)
  cp -r $file $new
done


#for dir in $(ls -d ../data/*/parent_data/TF* ../data/motifs ../data/plot_data)
#do
#  new=$(echo $dir | sed 's/data/data_subset/')
#  mkdir -p $(dirname $new)
#  cp -r $dir $new
#done
