# Written to look at the shRNA pertubation data
# downloaded by Emily and saved in the shared drive
# 10-11-16, JLW

import csv

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
CRKL_col=[]
for row in cr:
	if row['SM_Pert_Type']=='trt_xpr' and row['SM_Name']=='CRKL':
		CRKL_col.append(row)

column_ids=[row['id'] for row in CRKL_col]

# import l1k tools to read the enormous data object
import cmap.io.gct as gct
import cmap.io.plategrp as grp

path_to_gctx = share_dir+'GSE70138_Broad_LINCS_Level3_INF_mlr12k_n115209x22268_2015-12-31.gctx'
GCTObject = gct.GCT(path_to_gctx)
GCTObject.read(cid=column_ids,rid=row_ids)

GCTObject.write('../data/CRKL_lmark_matrix.gctx')

# Manipulate data further
data_array = GCTObject.matrix
symbols = GCTObject.get_row_meta('pr_gene_symbol')
cell_lines=GCTObject.get_column_meta('CL_Center_Specific_ID')
