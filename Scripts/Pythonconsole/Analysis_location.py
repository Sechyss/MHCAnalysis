import argparse
import pickle
import textwrap

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description=textwrap.dedent('''\
    Quick script to analyze if Kmers can bind to different regions compared to reference Pf3D7.
    ------------------------------------------
    '''))

parser.add_argument("-i", "--input", metavar='file.pickle', dest="dictionary", help="Dictionary with locations",
                    type=str)
args = parser.parse_args()

filepath = open(args.dictionary, 'rb')
dictionary_location = pickle.load(filepath)

for i in dictionary_location.keys():
    starting_pos = str(dictionary_location[i][0]).rsplit('_')[1]
    for y in dictionary_location[i]:
        if str(y).rsplit('_')[1] != starting_pos:
            print("".join(['Kmer ', str(i), ' binds to ', str(starting_pos), ' and ', str(y).rsplit('_')[1]]))
