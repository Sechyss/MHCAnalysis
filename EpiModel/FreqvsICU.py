import numpy as np
import pandas as pd
from MHCPipeline.Utils import flatten_array
import matplotlib.pyplot as plt


def main():
    excel_path = ('/Users/u2176312/OneDrive - University of '
                  'Warwick/SARS_COV_2/gisaid_variants_statistics_2024_02_19_1835.xlsx')

    excel_file = pd.ExcelFile(excel_path)

    table_Japan = {}
    table_Japan_2 = {}
    table_Japan_3 = {}
    table_UK = {}
    table_UK_2 = {}
    table_UK_3 = {}
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
            japan_data_1 = flatten_array(test.loc[japan_index + 1, columns].values.tolist())
            japan_data_2 = flatten_array(test.loc[japan_index, columns].values.tolist())
            japan_data = np.array(japan_data_2) / np.array(japan_data_1)
            uk_data_1 = flatten_array(test.loc[uk_index + 1, columns].values.tolist())
            uk_data_2 = flatten_array(test.loc[uk_index, columns].values.tolist())
            uk_data = np.array(uk_data_2) / np.array(uk_data_1)

            table_Japan[sheet_name] = japan_data_2
            table_Japan_2[sheet_name] = japan_data
            table_Japan_3[sheet_name] = japan_data_1
            table_UK[sheet_name] = uk_data_2
            table_UK_2[sheet_name] = uk_data
            table_UK_3[sheet_name] = uk_data_1

    # Create DataFrames
    table_Japan = pd.DataFrame.from_dict(table_Japan, orient='index', columns=columns).transpose()
    table_Japan_2 = pd.DataFrame.from_dict(table_Japan_2, orient='index', columns=columns).transpose()
    table_Japan_3 = pd.DataFrame.from_dict(table_Japan_3, orient='index', columns=columns).transpose()
    table_UK = pd.DataFrame.from_dict(table_UK, orient='index', columns=columns).transpose()
    table_UK_2 = pd.DataFrame.from_dict(table_UK_2, orient='index', columns=columns).transpose()
    table_UK_3 = pd.DataFrame.from_dict(table_UK_3, orient='index', columns=columns).transpose()

    writer = pd.ExcelWriter('/Users/u2176312/OneDrive - University of '
                            'Warwick/SARS_COV_2/Covid_strains_counts.xlsx', engine='openpyxl')
    table_Japan.to_excel(writer, sheet_name='Raw Counts Japan')
    table_Japan_3.to_excel(writer, sheet_name='Total Counts Japan')
    table_Japan_2.to_excel(writer, sheet_name='Proportion Counts Japan')

    table_UK.to_excel(writer, sheet_name='Raw Counts UK')
    table_UK_3.to_excel(writer, sheet_name='Total Counts UK')
    table_UK_2.to_excel(writer, sheet_name='Proportion Counts UK')

    writer.close()

    # Ensure the date column is a proper DateTime format before plotting
    new_Japan_df = table_Japan.div(table_Japan.sum(axis=1), axis=0)
    new_Japan_df.reset_index(inplace=True)
    new_Japan_df['index'] = pd.to_datetime(new_Japan_df['index'])  # Assuming 'index' is the date column
    new_Japan_df.rename(columns={'index': 'Date'}, inplace=True)

    new_UK_df = table_UK.div(table_UK.sum(axis=1), axis=0)
    new_UK_df.reset_index(inplace=True)
    new_UK_df['index'] = pd.to_datetime(new_UK_df['index'])
    new_UK_df.rename(columns={'index': 'Date'}, inplace=True)

    # %% Addition of hospitalization

    dataset = pd.read_csv('/Users/u2176312/OneDrive - University of '
                          'Warwick/SARS_COV_2/covid-hospitalizations.csv')
    Japan_hospitalizations = dataset[dataset['entity'] == 'Japan']
    Japan_hospitalizations['date'] = pd.to_datetime(Japan_hospitalizations['date'])
    Japan_hospitalizations = Japan_hospitalizations[
        Japan_hospitalizations['indicator'] == ('Daily hospital occupancy per '
                                                'million')]

    UK_hospitalizations = dataset[dataset['entity'] == 'United Kingdom']
    UK_hospitalizations['date'] = pd.to_datetime(UK_hospitalizations['date'])
    UK_hospitalizations = UK_hospitalizations[
        UK_hospitalizations['indicator'] == 'Daily hospital occupancy per million']
    # %% Combine the two figures
    fig, ax1 = plt.subplots(figsize=(12, 6))

    # Plot the new_Japan_df data on ax1
    new_Japan_df.plot(kind='area', x='Date', colormap='Paired', ax=ax1, linewidth=0)

    # Create the secondary axis after plotting on ax1
    ax2 = ax1.twinx()

    # Plot the hospitalizations on ax2
    Japan_hospitalizations.plot(kind='line', x='date', y='value', ax=ax2, legend=False, color='black')

    # Customize the plot for clarity
    ax1.set_xlabel('Date')
    ax1.set_ylabel('Variant proportions')
    ax2.set_ylabel('ICU Hospitalizations per million')
    ax1.legend(loc='upper left', prop={'size': 10})

    plt.tight_layout()  # Adjust layout for better spacing
    plt.savefig('/Users/u2176312/OneDrive - University of Warwick/SARS_COV_2/ICUvsDensity_Japan.pdf')

    fig, ax1 = plt.subplots(figsize=(12, 6))

    # Plot the new_Japan_df data on ax1
    new_UK_df.plot(kind='area', x='Date', colormap='Paired', ax=ax1, linewidth=0)

    # Create the secondary axis after plotting on ax1
    ax2 = ax1.twinx()

    # Plot the hospitalizations on ax2
    UK_hospitalizations.plot(kind='line', x='date', y='value', ax=ax2, legend=False, color='black')

    # Customize the plot for clarity
    ax1.set_xlabel('Date')
    ax1.set_ylabel('Variant proportions')
    ax2.set_ylabel('ICU Hospitalizations per million')
    ax1.legend(loc='upper left', prop={'size': 10})

    plt.tight_layout()  # Adjust layout for better spacing
    plt.savefig('/Users/u2176312/OneDrive - University of Warwick/SARS_COV_2/ICUvsDensity_UK.pdf')

    # %% Combine the two figures

    table_UK.reset_index(inplace=True)
    table_UK['index'] = pd.to_datetime(table_UK['index'])  # Assuming 'index' is the date column
    table_UK.rename(columns={'index': 'Date'}, inplace=True)
    fig, ax1 = plt.subplots(figsize=(12, 6))

    # Plot the new_Japan_df data on ax1
    table_UK.plot(kind='area', x='Date', colormap='Paired', ax=ax1, linewidth=0)

    # Create the secondary axis after plotting on ax1
    ax2 = ax1.twinx()

    # Plot the hospitalizations on ax2
    UK_hospitalizations.plot(kind='line', x='date', y='value', ax=ax2, legend=False, color='black')

    # Customize the plot for clarity
    ax1.set_xlabel('Date')
    ax1.set_ylabel('Variant counts')
    ax2.set_ylabel('ICU Hospitalizations per million')
    ax1.legend(loc='upper left', prop={'size': 10})

    plt.tight_layout()  # Adjust layout for better spacing
    plt.savefig('/Users/u2176312/OneDrive - University of Warwick/SARS_COV_2/ICUvsDensity_UK_2.pdf')

    table_Japan.reset_index(inplace=True)
    table_Japan['index'] = pd.to_datetime(table_Japan['index'])  # Assuming 'index' is the date column
    table_Japan.rename(columns={'index': 'Date'}, inplace=True)
    fig, ax1 = plt.subplots(figsize=(12, 6))

    # Plot the new_Japan_df data on ax1
    table_Japan.plot(kind='area', x='Date', colormap='Paired', ax=ax1, linewidth=0)

    # Create the secondary axis after plotting on ax1
    ax2 = ax1.twinx()

    # Plot the hospitalizations on ax2
    Japan_hospitalizations.plot(kind='line', x='date', y='value', ax=ax2, legend=False, color='black')

    # Customize the plot for clarity
    ax1.set_xlabel('Date')
    ax1.set_ylabel('Variant counts')
    ax2.set_ylabel('ICU Hospitalizations per million')
    ax1.legend(loc='upper left', prop={'size': 10})

    plt.tight_layout()  # Adjust layout for better spacing
    plt.savefig('/Users/u2176312/OneDrive - University of Warwick/SARS_COV_2/ICUvsDensity_Japan_2.pdf')


# %% Run main

if __name__ == '__main__':
    main()
