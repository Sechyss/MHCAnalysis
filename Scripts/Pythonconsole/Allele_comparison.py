import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn3, venn3_circles, venn2, venn2_circles
from venn import venn

fullproteinrun = pd.read_table('/Users/u2176312/OneDrive - University of Warwick/CSP/KenyaPopulation_mhc_i_pf3d7.txt')
filtered = fullproteinrun[fullproteinrun['rank'] <= 1]
filtered_region = set([filtered.loc[x]['allele'] for x, row in filtered.iterrows() if filtered.loc[x]['start'] >= 283
                       and filtered.loc[x]['end'] <= 342])
fullalleles = set(filtered['allele'])
region308342 = pd.read_table('/Users/u2176312/OneDrive - University of Warwick/'
                             'CSP/CSP_SNP_region/Kmer_CSP_predictionbinding_region_308-342.txt')
filtered308342 = region308342[region308342['rank'] <= 1]
region308 = set(filtered308342['allele'])
region283307 = pd.read_table('/Users/u2176312/OneDrive - University of Warwick/'
                             'CSP/CSP_SNP_region/Kmer_CSP_predictionbinding_region_283-307.txt')
filtered283307 = region283307[region283307['rank'] <= 1]
region283 = set(filtered283307['allele'])

v = venn3([filtered_region, region283, region308], ('ReferencePf3D7', 'Region283-307 (SNPKmers)',
                                                    'Region308-342 (SNPKmers)'))
c = venn3_circles([filtered_region, region283, region308], linestyle='dotted', linewidth=0.4)
plt.title('Matching HLAs for CSP (283-342 aa)')
plt.tight_layout()
#plt.savefig('/Users/u2176312/OneDrive - University of Warwick/CSP/Figures/VennKenya283-342region.pdf', dpi=300)
plt.show()

# %% Plot the matching alleles from Finland and Kenya

HLAsFinland = {'HLA-A*30:02', 'HLA-C*03:03', 'HLA-B*27:05', 'HLA-B*39:06', 'HLA-C*01:02', 'HLA-C*07:01', 'HLA-C*07:02',
               'HLA-B*35:01', 'HLA-C*05:01', 'HLA-C*07:04', 'HLA-B*44:02', 'HLA-B*55:01', 'HLA-A*23:01', 'HLA-B*27:02',
               'HLA-A*01:01', 'HLA-A*25:01', 'HLA-C*12:02', 'HLA-B*40:02', 'HLA-C*16:01', 'HLA-C*04:01', 'HLA-C*06:02',
               'HLA-A*68:06', 'HLA-B*41:02', 'HLA-C*14:02', 'HLA-C*17:01', 'HLA-B*47:01', 'HLA-B*40:01', 'HLA-B*44:03',
               'HLA-B*13:02', 'HLA-A*26:01', 'HLA-A*30:01', 'HLA-C*12:03', 'HLA-C*04:07', 'HLA-B*18:01', 'HLA-B*56:01',
               'HLA-A*11:01', 'HLA-B*08:01', 'HLA-A*02:01', 'HLA-B*41:01', 'HLA-B*39:03', 'HLA-A*24:02', 'HLA-C*15:02',
               'HLA-B*27:04', 'HLA-B*57:01', 'HLA-B*39:01', 'HLA-A*03:01', 'HLA-A*29:01', 'HLA-A*32:01'}

HLAsKenya = {'HLA-B*07:02', 'HLA-A*01:01', 'HLA-A*23:05', 'HLA-B*51:12', 'HLA-A*80:01', 'HLA-A*02:14', 'HLA-B*56:01',
             'HLA-A*24:06', 'HLA-C*07:27', 'HLA-C*03:05', 'HLA-A*68:01', 'HLA-B*18:01', 'HLA-C*04:01', 'HLA-C*08:01',
             'HLA-B*14:03', 'HLA-B*14:01', 'HLA-B*55:01', 'HLA-A*03:01', 'HLA-C*05:01', 'HLA-B*82:02', 'HLA-A*30:01',
             'HLA-A*29:03', 'HLA-A*74:01', 'HLA-A*02:05', 'HLA-B*53:02', 'HLA-C*01:02', 'HLA-C*15:02', 'HLA-B*13:01',
             'HLA-B*40:12', 'HLA-C*07:14', 'HLA-A*26:09', 'HLA-B*78:01', 'HLA-B*15:01', 'HLA-A*66:01', 'HLA-B*14:02',
             'HLA-C*03:10', 'HLA-A*23:01', 'HLA-A*32:01', 'HLA-A*26:01', 'HLA-A*26:05', 'HLA-B*38:01', 'HLA-B*51:01',
             'HLA-A*69:01', 'HLA-C*16:04', 'HLA-B*82:01', 'HLA-B*35:01', 'HLA-A*01:06', 'HLA-A*02:01', 'HLA-A*01:09',
             'HLA-B*40:16', 'HLA-B*35:03', 'HLA-A*24:02', 'HLA-A*02:06', 'HLA-B*07:07', 'HLA-A*43:01', 'HLA-C*16:01',
             'HLA-B*15:03', 'HLA-B*15:31', 'HLA-C*03:08', 'HLA-A*03:02', 'HLA-C*04:07', 'HLA-A*31:03', 'HLA-C*12:02',
             'HLA-A*33:03', 'HLA-C*08:02', 'HLA-C*16:02', 'HLA-A*66:02', 'HLA-B*44:03', 'HLA-B*67:01', 'HLA-B*52:01',
             'HLA-B*53:01', 'HLA-B*44:15', 'HLA-A*02:24', 'HLA-A*24:13', 'HLA-B*07:17', 'HLA-B*57:03', 'HLA-C*06:11',
             'HLA-C*14:03', 'HLA-B*08:01', 'HLA-B*35:34', 'HLA-B*44:07', 'HLA-A*02:02', 'HLA-B*37:01', 'HLA-A*34:01',
             'HLA-B*47:02', 'HLA-A*02:34', 'HLA-C*15:05', 'HLA-B*15:10', 'HLA-C*02:02', 'HLA-B*44:26', 'HLA-A*36:01',
             'HLA-C*08:04', 'HLA-B*39:01', 'HLA-C*06:02', 'HLA-B*49:01', 'HLA-B*41:01', 'HLA-B*27:05', 'HLA-B*15:17',
             'HLA-B*42:01', 'HLA-A*02:36', 'HLA-B*81:01', 'HLA-A*30:06', 'HLA-B*07:05', 'HLA-B*40:02', 'HLA-B*73:01',
             'HLA-B*15:37', 'HLA-C*03:04', 'HLA-A*24:23', 'HLA-B*45:02', 'HLA-B*50:01', 'HLA-A*02:04', 'HLA-B*15:67',
             'HLA-B*41:02', 'HLA-B*15:16', 'HLA-C*04:04', 'HLA-A*01:02', 'HLA-A*02:25', 'HLA-A*30:03', 'HLA-C*07:05',
             'HLA-C*14:02', 'HLA-B*27:03', 'HLA-B*39:03', 'HLA-A*02:17', 'HLA-C*12:03', 'HLA-B*35:08', 'HLA-B*58:02',
             'HLA-B*15:18', 'HLA-A*30:02', 'HLA-A*11:01', 'HLA-B*57:01', 'HLA-B*58:01', 'HLA-A*24:03', 'HLA-A*34:02',
             'HLA-B*41:03', 'HLA-A*31:04', 'HLA-B*35:02', 'HLA-C*07:04', 'HLA-A*26:12', 'HLA-B*42:02', 'HLA-B*47:03',
             'HLA-B*51:08', 'HLA-C*03:02', 'HLA-A*02:26', 'HLA-B*40:01', 'HLA-A*74:03', 'HLA-B*53:05', 'HLA-A*29:02',
             'HLA-A*25:01', 'HLA-A*68:02', 'HLA-C*03:03', 'HLA-A*02:11', 'HLA-B*13:02', 'HLA-B*57:02', 'HLA-C*07:01',
             'HLA-B*18:07', 'HLA-B*44:02', 'HLA-B*18:03', 'HLA-A*29:01', 'HLA-A*26:08', 'HLA-C*17:01', 'HLA-B*47:01',
             'HLA-B*15:55', 'HLA-A*01:03', 'HLA-B*15:09', 'HLA-C*18:01', 'HLA-B*44:05', 'HLA-B*13:03', 'HLA-B*39:10',
             'HLA-B*45:01', 'HLA-A*33:01', 'HLA-C*02:10', 'HLA-A*29:15', 'HLA-A*30:04', 'HLA-C*15:03', 'HLA-C*07:02'}

HLAsPapuaNewGuinea = {'HLA-A*02:01', 'HLA-A*02:01', 'HLA-A*02:01', 'HLA-A*02:01', 'HLA-A*02:06', 'HLA-A*02:06',
                      'HLA-A*02:06', 'HLA-A*02:06', 'HLA-A*03:01', 'HLA-A*03:01', 'HLA-A*03:01', 'HLA-A*03:01',
                      'HLA-A*11:01', 'HLA-A*11:01', 'HLA-A*11:01', 'HLA-A*11:01', 'HLA-A*24:02', 'HLA-A*24:02',
                      'HLA-A*24:02', 'HLA-A*24:02', 'HLA-A*24:07', 'HLA-A*24:07', 'HLA-A*24:07', 'HLA-A*24:07',
                      'HLA-A*25:01', 'HLA-A*25:01', 'HLA-A*25:01', 'HLA-A*25:01', 'HLA-A*26:01', 'HLA-A*26:01',
                      'HLA-A*26:01', 'HLA-A*26:01', 'HLA-A*32:01', 'HLA-A*32:01', 'HLA-A*32:01', 'HLA-A*32:01',
                      'HLA-A*34:01', 'HLA-A*34:01', 'HLA-A*34:01', 'HLA-A*34:01', 'HLA-A*68:02', 'HLA-A*68:02',
                      'HLA-A*68:02', 'HLA-A*68:02', 'HLA-B*07:02', 'HLA-B*07:02', 'HLA-B*07:02', 'HLA-B*07:02',
                      'HLA-B*13:01', 'HLA-B*13:01', 'HLA-B*13:01', 'HLA-B*13:01', 'HLA-B*14:01', 'HLA-B*14:01',
                      'HLA-B*14:01', 'HLA-B*14:01', 'HLA-B*15:01', 'HLA-B*15:01', 'HLA-B*15:01', 'HLA-B*15:01',
                      'HLA-B*15:05', 'HLA-B*15:05', 'HLA-B*15:05', 'HLA-B*15:05', 'HLA-B*15:06', 'HLA-B*15:06',
                      'HLA-B*15:06', 'HLA-B*15:06', 'HLA-B*15:21', 'HLA-B*15:21', 'HLA-B*15:21', 'HLA-B*15:21',
                      'HLA-B*15:25', 'HLA-B*15:25', 'HLA-B*15:25', 'HLA-B*15:25', 'HLA-B*15:36', 'HLA-B*15:36',
                      'HLA-B*15:36', 'HLA-B*15:36', 'HLA-B*18:01', 'HLA-B*18:01', 'HLA-B*18:01', 'HLA-B*18:01',
                      'HLA-B*27:04', 'HLA-B*27:04', 'HLA-B*27:04', 'HLA-B*27:04', 'HLA-B*27:06', 'HLA-B*27:06',
                      'HLA-B*27:06', 'HLA-B*27:06', 'HLA-B*35:05', 'HLA-B*35:05', 'HLA-B*35:05', 'HLA-B*35:05',
                      'HLA-B*39:01', 'HLA-B*39:01', 'HLA-B*39:01', 'HLA-B*39:01', 'HLA-B*39:03', 'HLA-B*39:03',
                      'HLA-B*39:03', 'HLA-B*39:03', 'HLA-B*39:06', 'HLA-B*39:06', 'HLA-B*39:06', 'HLA-B*39:06',
                      'HLA-B*40:01', 'HLA-B*40:01', 'HLA-B*40:01', 'HLA-B*40:01', 'HLA-B*40:02', 'HLA-B*40:02',
                      'HLA-B*40:02', 'HLA-B*40:02', 'HLA-B*40:10', 'HLA-B*40:10', 'HLA-B*40:10', 'HLA-B*40:10',
                      'HLA-B*44:02', 'HLA-B*44:02', 'HLA-B*44:02', 'HLA-B*44:02', 'HLA-B*48:01', 'HLA-B*48:01',
                      'HLA-B*48:01', 'HLA-B*48:01', 'HLA-B*51:01', 'HLA-B*51:01', 'HLA-B*51:01', 'HLA-B*51:01',
                      'HLA-B*53:01', 'HLA-B*53:01', 'HLA-B*53:01', 'HLA-B*53:01', 'HLA-B*55:02', 'HLA-B*55:02',
                      'HLA-B*55:02', 'HLA-B*55:02', 'HLA-B*56:01', 'HLA-B*56:01', 'HLA-B*56:01', 'HLA-B*56:01',
                      'HLA-B*56:02', 'HLA-B*56:02', 'HLA-B*56:02', 'HLA-B*56:02', 'HLA-C*01:02', 'HLA-C*01:02',
                      'HLA-C*01:02', 'HLA-C*01:02', 'HLA-C*02:02', 'HLA-C*02:02', 'HLA-C*02:02', 'HLA-C*02:02',
                      'HLA-C*03:03', 'HLA-C*03:03', 'HLA-C*03:03', 'HLA-C*03:03', 'HLA-C*03:04', 'HLA-C*03:04',
                      'HLA-C*03:04', 'HLA-C*03:04', 'HLA-C*04:01', 'HLA-C*04:01', 'HLA-C*04:01', 'HLA-C*04:01',
                      'HLA-C*04:03', 'HLA-C*04:03', 'HLA-C*04:03', 'HLA-C*04:03', 'HLA-C*07:01', 'HLA-C*07:01',
                      'HLA-C*07:01', 'HLA-C*07:01', 'HLA-C*07:02', 'HLA-C*07:02', 'HLA-C*07:02', 'HLA-C*07:02',
                      'HLA-C*08:01', 'HLA-C*08:01', 'HLA-C*08:01', 'HLA-C*08:01', 'HLA-C*08:02', 'HLA-C*08:02',
                      'HLA-C*08:02', 'HLA-C*08:02', 'HLA-C*12:02', 'HLA-C*12:02', 'HLA-C*12:02', 'HLA-C*12:02',
                      'HLA-C*12:03', 'HLA-C*12:03', 'HLA-C*12:03', 'HLA-C*12:03', 'HLA-C*15:02', 'HLA-C*15:02',
                      'HLA-C*15:02', 'HLA-C*15:02'}

y = venn3([HLAsKenya, HLAsFinland, HLAsPapuaNewGuinea], ('Kenya populations', 'Finland populations',
                                                         'Papua New Guinea populations'),
          set_colors=('green', 'gold', 'red'))
z = venn3_circles([HLAsKenya, HLAsFinland, HLAsPapuaNewGuinea], linestyle='dotted')

plt.tight_layout()

#plt.savefig('/Users/u2176312/OneDrive - University of Warwick/CSP/Figures/VennKenyaFinlandPapuaNewGuineacomparison.pdf',
#            dpi=300)
plt.show()

# %%

sets = {
    'Reference Pf3D7': filtered_region,
    'Region283-307 (SNPKmers)': region283,
    'Region308-342 (SNPKmers)': region308,
}

fig, ax = plt.subplots(figsize=(10, 10))
venn(sets, ax=ax, legend_loc='lower left')

plt.show()
