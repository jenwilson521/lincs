# Written to look at the crizotinib pertubation data
# downloaded by Emily and saved in the shared drive
# 10-11-16, JLW

import csv, itertools, pickle
from collections import defaultdict
import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm 
import numpy as np

# directories and file names for downloaded files
share_dir = '/share/PI/rbaltman/erflynn/l1000/data/'
col_file = share_dir+'geo_col_metadata_l3.csv'
row_file = share_dir+'geo_row_metadata_l3.csv'
# row headers
row_headers_raw='ind,id,pr_analyte_id,pr_analyte_num,pr_bset_id,pr_gene_id,pr_gene_symbol,pr_gene_title,pr_is_bing,pr_is_inf,pr_is_lmark,pr_lua_id,pr_model_id,pr_pool_id'
col_headers_raw='ind,id,CL_Center_Specific_ID,RN_Target_Gene_ID,SM_Center_Compound_ID,SM_Dose,SM_Dose_Unit,SM_LINCS_ID,SM_Name,SM_Pert_Type,SM_Time,SM_Time_Unit,det_plate,det_well'
rw_hdrs=row_headers_raw.split(',')
cl_hdrs=col_headers_raw.split(',')

# create dictionary readers and parse data for 978 measured genes
# start with CRKL, get all redundant shRNA treatments
rr = csv.DictReader(open(row_file,'r'))
lmark_rows=[]
for row in rr:
	if row['pr_is_lmark']=='Y':
		lmark_rows.append(row)
row_ids=[row['id'] for row in lmark_rows]

cr=csv.DictReader(open(col_file,'r'))
criz_col=[]
for row in cr:
	if row['SM_Pert_Type']=='trt_cp' and row['SM_Name']=='crizotinib':
		criz_col.append(row)

column_ids=[row['id'] for row in criz_col]

# import l1k tools to read the enormous data object
import cmap.io.gct as gct
import cmap.io.plategrp as grp

path_to_gctx = share_dir+'GSE70138_Broad_LINCS_Level3_INF_mlr12k_n115209x22268_2015-12-31.gctx'
GCTObject = gct.GCT(path_to_gctx)
GCTObject.read(cid=column_ids,rid=row_ids)

GCTObject.write('../data/criz_lmark_matrix.gctx')

# Manipulate data further
data_array = GCTObject.matrix
symbols = GCTObject.get_row_meta('pr_gene_symbol')
cell_lines = GCTObject.get_column_meta('CL_Center_Specific_ID')
time_points = GCTObject.get_column_meta('SM_Time')
doses = GCTObject.get_column_meta('SM_Dose')
dose_unit = GCTObject.get_column_meta('SM_Dose_Unit')

# print data information
print('extracted crizotinib pertubation data')
print('cell lines: ',set(cell_lines))
print('time points: ',set(time_points))

# plot some example genes to show trends in raw data
gene_data=zip(symbols,data_array)
unique_CL = ['HEPG2', 'MDAMB231', 'A375', 'PC3', 'MCF7', 'HA1E', 'MCF10A', 'HT29', 'BT20', 'HS578T', 'HCC515', 'SKBR3', 'A549']
ass_fea = zip(cell_lines,time_points,doses) # assembled features
### color=cm.rainbow(np.linspace(0,1,len(unique_CL)))
colors = itertools.cycle(['r','b','g'])
# Create scatter plots looking at individual genes across cell lines
print('starting plots')
for (gene,row_data) in gene_data[0:10]:
	fig,ax=plt.subplots()
	### for (CL,clr) in zip(unique_CL,color):
	for CL in unique_CL:
		x = unique_CL.index(CL)
		col_inds = [i for (i,(cl,tp,ds)) in enumerate(ass_fea) if cl==CL and tp=='24' and ds=='1.11']
		yd = [row_data[a] for a in col_inds]
		xd = [x]*len(col_inds)
		ax.scatter(xd,yd,color=next(colors)) #c=clr)

	ax.set_xticks(np.arange(-1,14,1.0))
	ax.set_xticklabels(['']+unique_CL+[''],rotation=45)
	ax.set_title(gene+' exp across cell lines')
	ax.set_xlabel('cell lines')
	ax.set_ylabel('landmark gene expression')
	ax.set_ylim([0,15])
	ax.set_xlim([-1, 14])
	plt.savefig('../plots/'+gene+'_24hr_1.11uM.png',format='png')

# write to output for use with mutual information heuristic
print('summarizing data by cell line')
missing_data = [] # keep track of genes missing the 3 hr time point
#st = sorted(list(set(time_points)),reverse=True)
sd = sorted(list(set(doses)))
for CL in unique_CL:
	print('cell line ',CL)
	#outf=open('../data/'+CL+'_1.11uM.txt','w')
	outf=open('../data/'+CL+'_24hr.txt','w')
	header=['gene']+[str(a)+'uM' for a in sd]+['\n'] 
	outf.write('\t'.join(header))
	for (gene,row_data) in gene_data:
		d_data = []
		for dose in sd:
			col_inds = [i for (i,(cl,tp,ds)) in enumerate(ass_fea) if cl==CL and tp=='24' and ds==dose]
			yd = [row_data[a] for a in col_inds]
			d_data.append(np.median(yd))
			if len(yd)==0:
				missing_data.append((CL,gene,dose))
		output_data = [gene]+[str(a) for a in d_data]+['\n']
		outf.write('\t'.join(output_data))
	outf.close()
pickle.dump(missing_data,open('../data/doses_missing.pkl','w')) 


