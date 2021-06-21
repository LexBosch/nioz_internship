""" Author: Lex Bosch
    Date: 2021/03/11

    This script generates histograms with different options to generate them.
"""

import json
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats


def reject_outliers(data, m=6):
    return data[abs(data - np.mean(data)) < m * np.std(data)]


def read_contig_file(filename: str, div_by_len: bool = False):
    """ Reads and processes given contig file. Creates 2d list [[#_reads]]
        If div_by_len is set to True, #_reads is divided by the length of the contig
    
    :param filename: File to be processed. Json file with format {Contig_name:[sequence_length, #_mapped_reads]}
    :param div_by_len: Boolean, optional. If True, #_reads value will be divided by the sequence length
    :return: Returns list in format [[#_reads]]
    """
    with open(filename, "r") as open_json:
        complete_list = []
        json_dict = json.load(open_json)
        for contig in json_dict:
            list_per_contig = []
            for contig_name in contig.keys():
                if div_by_len:
                    list_per_contig.append((contig[contig_name][1]) / (contig[contig_name][0]))
                else:
                    list_per_contig.append(contig[contig_name][1])
            complete_list.append(list_per_contig)
    return complete_list


def create_hist(elem_list: list, save_path: str = None, alt_plot_text: str = None):
    """ Creates a series of histograms using the given list of elements

    :param elem_list: 2d list containing #_reads for each histogram
    :param save_path: string, optional. path to the folder in which to save the plots
    :param save_path:
    :return: None
    """
    for count, single in enumerate(elem_list):
        single_pd = np.log10(np.array(single)+1)
        n, x, _ = plt.hist(single_pd, bins=100, edgecolor="black", density=True)
        density = stats.gaussian_kde(single_pd)

        plt.plot(x, density(x))
        plt.xlabel("Log 10 of amount of mapped read-segments per contig")
        plt.ylabel("Density of amount of contigs with a\ngiven number of mapped reads")
        if alt_plot_text:
            plt.title(
                f"Histogram showing the distribution of mapped read \nsegments divided by the length in log10of the sequence per contig of all bins")
        else:
            plt.title('A histogram showing the distribution of mapped read segments per contig')

        if save_path is not None:
            plt.savefig('output/plots/{}/gen_histogram{}.png'.format(save_path, count + 1))
        plt.show()


# SET THESE VARIABLES!
contig_info_path = ""  # Path to contig info path




com_list = read_contig_file(contig_info_path, False)

less_list = []
for i in com_list:
    for elem in i:
        less_list.append(elem)


create_hist([reject_outliers(np.array(less_list))])