# qqman
My QQ/Manhattan-plotter written in Python.

# Command-line usage
To see the help page:
python qqman.py -h
You’ll probably end up running something like this:
python qqman.py --tsv my_sumstats.txt.gz --qq_name MyQQ.png –man_name MyManhattan.png --col_chr chromosome --col_bp basepair --col_rsid snp

There’s an additional ‘options’ argument which I’ve designed to maintain full functionality from the command line, without adding too much rarely used arguments to the list. The input is a JSON-style character string (single quotes on the outside, double quotes on the inside) with your own options (see python qqman.py –adv_help). So for example, people often hate the colors I use, so to change the colors of the manhattan to blue and red:
python qqman.py --tsv my_sumstats.txt.gz --qq_name MyQQ.png --man_name MyManhattan.png --options '{"qq":{}, "man": {"pointcolor": ["red", "blue"]}}'
If you want to know which colors you can use: https://matplotlib.org/tutorials/colors/colors.html

