input=($(echo $@))
clusters=${input[0]}
files=$(echo ${input[@]} | tr ' ' '\n' | grep "[a-z]*[0-9].clusters" | paste -sd " ")
#echo -e "$files\n$clusters"
name_fmt=$(realpath ${input[@]} | egrep -o '\w+[12]\.' | sed 's/\.//') &&
#echo $name_fmt
echo chr start end $name_fmt | tr ' ' '\t' &&
bedtools annotate -counts -i <(cat $clusters | cut -f 1,2,3) -files $files
