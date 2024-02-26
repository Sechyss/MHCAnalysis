import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from MHCPipeline.Utils import flatten_array

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
table_Japan = pd.DataFrame.from_dict(table_Japan)
table_UK = pd.DataFrame.from_dict(table_UK)
table_UK['Dates'] = pd.to_datetime(columns, format='%Y-%m-%d', errors='coerce', dayfirst=True)
table_Japan['Dates'] = pd.to_datetime(columns, format='%Y-%m-%d', errors='coerce', dayfirst=True)

#%% Plotting the results

sns.kdeplot(data=table_Japan, multiple='fill', common_norm=False)
plt.show()
