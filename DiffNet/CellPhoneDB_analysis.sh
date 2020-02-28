# cellphonedb method statistical_analysis --help

inpath=$1
outpath=$2
iter=$3

###############
infile=expression_matrix.txt
inmeta=meta_data.txt

cellphonedb method statistical_analysis $inpath/$inmeta $inpath/$infile --counts-data hgnc_symbol --threads 60 --output-path $outpath --iterations $iter







