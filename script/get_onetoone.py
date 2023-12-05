import sys, re
import pandas as pd
from argparse import ArgumentParser


def get_option():
	argparser = ArgumentParser()
	argparser.add_argument('ortho', help='One-to-one ortholog file to be sanitized')
	argparser.add_argument('--sanitize', action='store_true', help='Sanitizing ortholog table')
	argparser.add_argument('--manualcuration', type=str,
						default='',
						help='Specify the file containing a list of manually curated one-to-one orthologs')
	return argparser.parse_args()

def sanitize(df_org):

	df0 = df_org[~df_org.duplicated()]
	colnames = df0.columns.values

	for c in colnames:
		df_tmp = df0[~df0.duplicated(subset=c, keep=False)]
		df0 = df_tmp

	return(df0)

def addwhite(df_org, in_white):
	df_white = pd.read_csv(in_white, sep='\t')
	df00 = pd.concat([df_org, df_white], axis=0, ignore_index=True)

	return(df00)


if __name__ == '__main__':
	args = get_option()

	df_in = pd.read_csv(args.ortho, sep='\t')
	colnames = df_in.columns.values
	colnames_new = [s.replace('gene stable ID', 'gene ID') for s in colnames]
	df_in.columns = colnames_new

	df_hml = df_in.filter(like='homology type', axis=1)
	df1 = df_in.filter(regex='gene ID|gene name', axis=1).fillna(method='ffill', axis=1).loc[(df_hml == 'ortholog_one2one').all(axis=1), :]

	if args.manualcuration != '':
		df1 = addwhite(df1, args.manualcuration)

	if args.sanitize:
		df1 = sanitize(df1)

	df1.to_csv(sys.stdout, header=True, index=False, sep='\t')


