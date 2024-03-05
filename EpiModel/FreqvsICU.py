import pandas as pd
from MHCPipeline.Utils import flatten_array
import matplotlib.pyplot as plt

excel_path = ('/Users/u2176312/OneDrive - University of '
              'Warwick/SARS_COV_2/gisaid_variants_statistics_2024_02_19_1835.xlsx')

excel_file = pd.ExcelFile(excel_path)

table_Japan = {}
table_UK = {}
columns = []

for sheet_name in excel_file.sheet_names:
    if sheet_name == 'README':
        # Handle README sheet content if needed
        pass
    else:
        test = pd.read_excel(excel_file, sheet_name)
        columns = test.columns[2:]  # Assuming columns start from index 2
        japan_index = test.loc[test['Unnamed: 0'] == 'Japan'].index
        uk_index = test.loc[test['Unnamed: 0'] == 'United Kingdom'].index

        # Extract data for Japan and UK, potentially adjust indexing as needed
        japan_data = flatten_array(test.loc[japan_index + 1, columns].values.tolist())
        uk_data = flatten_array(test.loc[uk_index + 1, columns].values.tolist())

        table_Japan[sheet_name] = japan_data
        table_UK[sheet_name] = uk_data

# Create DataFrames, potentially handle errors in `pd.to_datetime()`
table_Japan = pd.DataFrame.from_dict(table_Japan, orient='index', columns=columns).transpose()

table_UK = pd.DataFrame.from_dict(table_UK, orient='index', columns=columns).transpose()
#writer = pd.ExcelWriter('/Users/u2176312/OneDrive - University of '
#                        'Warwick/SARS_COV_2/Dataset_Japan_UK.xlsx', engine='openpyxl')
#table_Japan.to_excel(writer, sheet_name='Japan')
#table_UK.to_excel(writer, sheet_name='United Kingdom')
#writer.close()

#new_Japan_df = pd.DataFrame(index=table_Japan.index, columns=table_Japan.columns)

new_Japan_df = table_Japan.div(table_Japan.sum(axis=1), axis=0)
new
new_Japan_df.plot( kind='area', use_index=True, xticks=new_Japan_df.index, colormap='Paired')
plt.show()
