import pandas as pd

new_data = pd.read_excel('/Users/u2176312/OneDrive - University of Warwick/CSP/AllelePops/FilteredDataAllele.xlsx',
                         sheet_name='FilteredAlleles', index_col=0)

listallelesnetpanmhc = pd.read_table(
    '/Users/u2176312/OneDrive - University of Warwick/CSP/ListMHC_humans_netmhcpan.txt',
    header=None, sep='\t')

All_HLAs = new_data.index.tolist()
temp_hlas = ['HLA-' + str(x) for x in All_HLAs]
netmhcalleles = list(listallelesnetpanmhc[0])
temp_list = [x.replace(' ', '') for x in netmhcalleles]
common_mhc = list(set(temp_hlas) & set(temp_list))
print('" "'.join(set(common_mhc[:])))
