import gzip

import seaborn as sns
from matplotlib import pyplot as plt


def get_vcf_names(vcf_path):
    with gzip.open(vcf_path, "rt") as ifile:
        for line in ifile:
            if line.startswith("#CHROM"):
                vcf_names = [x for x in line.split('\t')]
                break
    ifile.close()

    return vcf_names


def plot_windowed_variant_density(pos, qual, title=None, chromosome=None, y_title=None):
    # use window midpoints as x coordinate
    x = pos

    # compute variant density in each window
    y = qual

    # plot
    fig, ax = plt.subplots(figsize=(12, 6))
    plt.rcParams["axes.labelweight"] = "bold"
    sns.despine(ax=ax, offset=10)
    ax.plot(x, y)

    ax.set_xlabel(chromosome+' Gene position (bp)', weight='bold')
    ax.set_ylabel(y_title)
    if title:
        ax.set_title(title)
