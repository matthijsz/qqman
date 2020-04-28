# qqman
My QQ/Manhattan-plotter written in Python. Depending on the size of the input data this is about 2-3x faster than doing the equivalent in R (bigger input files result in bigger improvements).<br>
It can be used directly in your own Python script, or it can be run from the command line. <br>

# Command-line usage
To see the help page:<br>
python qqman.py -h<br>
You’ll probably end up running something like this:<br>
`python qqman.py --tsv my_sumstats.txt.gz --qq_name MyQQ.png –man_name MyManhattan.png --col_chr chromosome --col_bp basepair --col_rsid snp`<br>
<br>
There’s an additional ‘options’ argument which I’ve designed to maintain full functionality from the command line, without adding too much rarely used arguments to the list. The input is a JSON-style character string (single quotes on the outside, double quotes on the inside) with your own options (see python qqman.py –adv_help). So for example, people often hate the colors I use, so to change the colors of the manhattan to blue and red:<br>
`python qqman.py --tsv my_sumstats.txt.gz --qq_name MyQQ.png --man_name MyManhattan.png --options '{"qq":{}, "man": {"pointcolor": ["red", "blue"]}}'`<br>
If you want to know which colors you can use: https://matplotlib.org/tutorials/colors/colors.html

