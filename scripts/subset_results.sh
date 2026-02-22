# exclude .npz, .png, .h5, .html
for file in $(find ../results/* | grep -v "\.npz" | grep -v "\.png" | grep ".*\....$")
do
  new=$(echo $file | sed 's/results/result_subset/g')
  mkdir -p $(dirname $new)
  cp -r $file $new
done
