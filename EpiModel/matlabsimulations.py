import argparse
import re
import textwrap

import pandas as pd
import scipy.io as sc
from matplotlib import pyplot as plt

# =============================================================================

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description=textwrap.dedent('''\
Analysis of Matlab simulations results
------------------------------------------
Tasks added so far:
haplotable -- Takes the results from the HLA_haplo_tbl variable from the simulation and finds the haplotypes that 
              were not in the original population, returning the time point when they reach a given threshold of the 
              original population.
plottable  -- Takes the results from the HLA_haplo_tbl variable from the simulation and plots the time where a given 
              T-cell response probability reaches a given threshold of the original population.
stacked    -- Takes the results from the HLA_haplo_tbl variable from the simulation and plots a stacked bar chart of the 
              successful simulations that reached the threshold of the original population.
'''))
# Parse command line arguments
# -i FILE -t TASK -h HELP -o OUTPUT -th Threshold

parser.add_argument("-i", "--input", metavar='file.mat', dest="input", help="File to load", type=str)
parser.add_argument("-t", "--task", dest="task", help="Task to perform", type=str, choices=["haplotable",
                                                                                            "plottable",
                                                                                            "stacked"])
parser.add_argument("-o", "--output", metavar='file.tsv', dest="output", help="Output file to save", type=str)
parser.add_argument("-th", "--threshold", dest="threshold", help="Threshold for haplotable task", type=float)

args = parser.parse_args()

# =============================================================================


print('Opening your target file...' + str(args.input))


def analysis_hpl_table():
    target = sc.loadmat(args.input)
    tablehplt = pd.DataFrame(target["HLA_haplo_tbl"])
    limitpop = int(tablehplt.iloc[0].max())  # Select the maximum value in the  population
    thresholdpop = int(limitpop * args.threshold)  # Apply the threshold to the population
    newalleles = tablehplt.loc[:, (tablehplt.iloc[0] == 0)]
    year = 0
    for index, row in newalleles.iterrows():
        if newalleles.iloc[index].max() >= thresholdpop:  # If the value is greater or equal than the threshold breaks
            # the loop
            year = (index + 1) * 100
            break
        else:  # If the value never reaches the thresolhold population index becomes nan
            year = 'NaN'

    print('Saving the final data in ' + args.output)
    with open(args.output, 'a') as f:
        f.write(str(re.findall('Simulations/(.*)repeat', args.input)[0]) + '\t' +  # Probability of T-Cell response
                str(year) + '\t' +  # Years until threshold reached
                str(re.findall('repeat(.*).mat', args.input)[0]) + '\t' +
                str(args.threshold * 100) + '%' + '\n')  # Simulation repeat number


def plot_results_haplotable():
    print('Plotting the table of results from the HLA_haplo_tbl variable from the different simulations...')
    table = pd.read_table(args.input, names=['ProbabilityTcell', 'Years', 'Repetition', 'ThresholdPop'], header=None)
    probabilities = list(set(table['ProbabilityTcell']))
    probabilities.sort()
    average = []
    std = []
    for value in probabilities:
        test = table[table['ProbabilityTcell'] == value]
        average.append(test['Years'].mean())
        std.append(test['Years'].std())

    plt.figure(figsize=(10, 5))
    plt.rcParams["axes.labelweight"] = "bold"

    plt.errorbar(probabilities, average, yerr=std, capsize=2, capthick=1,
                 fmt='-o',
                 markerfacecolor='none',
                 color='black',
                 markeredgecolor='black',
                 )

    plt.xlabel('Probability of T-cell response')
    plt.ylabel('Years until threshold reached')
    plt.tight_layout()
    print('Saving the final plot in ' + args.output)
    plt.savefig(args.output, dpi=300)


def stacked_plot_results_haplotable():
    table = pd.read_table(args.input, names=['ProbabilityTcell', 'Years', 'Repetition', 'ThresholdPop'], header=None)
    probabilities = list(set(table['ProbabilityTcell']))
    probabilities.sort()
    temp = {}
    for values in probabilities:
        filteredtable = table[table['ProbabilityTcell'] == values]
        totalcases = len(filteredtable['Years'])
        yearcolumn = filteredtable['Years']
        casesbelow5000 = yearcolumn[yearcolumn < 5000].count() / totalcases
        casesover5000 = len(filteredtable[filteredtable['Years'].between(5000, 10000, inclusive='both')]) / totalcases
        casesover10000 = len(
            filteredtable[filteredtable['Years'].between(10000, 15000, inclusive='neither')]) / totalcases
        casesover15000 = yearcolumn[15000 <= yearcolumn].count() / totalcases
        temp.update({values: [casesbelow5000, casesover5000, casesover10000, casesover15000]})

    df_to_plot = pd.DataFrame.from_dict(temp, orient='index', columns=['Sim<5000y', 'Sim>5000y',
                                                                       'Sim>10000y', 'Sim>15000y'])

    # Make a coordinates by cycling through the colors you care about
    # to match the length of your data.
    colors_bar = ['#1E555C', '#E28413', '#000022', '#DE3C4B']
    df_to_plot.plot(kind='bar', stacked=True, mark_right=False, color=colors_bar)
    plt.xlabel('Probability of T-cell response', weight='bold')
    plt.ylabel('%simulations - threshold population', weight='bold')
    plt.legend(loc='best', prop={'size': 8})
    print('Saving the final plot in ' + args.output)
    plt.savefig(args.output, dpi=300)


def main():
    task = {'haplotable': analysis_hpl_table, 'plottable': plot_results_haplotable,
            'stacked': stacked_plot_results_haplotable}
    choosetask = str(args.task)
    task[choosetask]()


if __name__ == '__main__':
    main()
