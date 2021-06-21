""" Author: Lex Bosch
    Date: 2021/03/23

    This script generates stacked bar pltos with different options to generate them.
"""


import json
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from os import listdir
from Bio import SeqIO
from statistics import median, mean, mode
from math import log


def read_bin_file(folder: str) -> list:
    """ Reads the binned files in the folder given

    :param folder: filepath which folder contains all the binning files
    :return: List containing each bin and their contigs
    """
    file_index = listdir(folder)
    full_list = []
    for file in file_index:
        partial_dict = {}
        cur_file = (folder + "/" + file)
        for record in SeqIO.parse(cur_file, "fasta"):
            partial_dict[record.description] = (str(record.seq))
        full_list.append(partial_dict)
    return full_list


def get_tpmfromfile(quant_filepath: str) -> pd.DataFrame:
    """ Reads the given file and places this in a pandas database

    :param quant_filepath: filepath to the quantification file
    :return: pandas dataframe with the quantification data
    """
    data_frame = pd.read_csv(quant_filepath, sep="\t")
    return data_frame


def read_contig_file(filename: str, div_by_len: bool = False) -> list:
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


def place_in_dict(data_list: list) -> dict:
    """ Generates a dictionary containing mean, median and sum of the given values

    :param data_list: List containing floats
    :return: dictionary containing the mean, median and sum of the given list
    """
    return {"mean": mean(data_list), "median": median(data_list), "sum": sum(data_list)}


def calculate_realtive(amounts_list: list) -> list:
    """ Returns the relative abundance values of the given list.
        (For instance, a list of [5, 5, 10] would become a list of [0.25, 0.25, 0.5])

    :param amounts_list: List containing numeric values
    :return: returns list with the relative abundances
    """
    return [i / sum(amounts_list) for i in amounts_list]


def create_bars(bar_value_dict: dict, savepath: str = None, plot_name: str = "gen_plot",
                sample_name: str = None) -> None:
    """ Creates stacked bar plots of the given data

    :param bar_value_dict: Dictionary containing data used to create stacked bar plots.
    :param savepath: Path to save the plots to. If not given, will not save the plots.
    :param plot_name: Name of the plot when saved. If not given, defaults to 'gen_plot'
    :param sample_name: Name of the sample of the plot. If not given, will not be implemented in the figure
    :return: Returns nothing, but creates and (optionally saves) plots
    """
    for pltnum, anlysis_type in enumerate(bar_value_dict):
        fig, ax = plt.subplots()
        layer_height = [0 for _ in bar_value_dict[anlysis_type].keys()]
        for sec_type in bar_value_dict[anlysis_type].keys():
            for sincount, _ in enumerate(bar_value_dict[anlysis_type][sec_type]):
                data_layer = [bar_value_dict[anlysis_type][i][sincount] for i in bar_value_dict[anlysis_type].keys()]
                labels = bar_value_dict[anlysis_type].keys()
                ax.bar(labels, data_layer, 0.35, bottom=layer_height, label=real_bins[sincount])
                temp_coupled_list = [layer_height, data_layer]
                layer_height = [sum(x) for x in zip(*temp_coupled_list)]
            if sample_name is None:
                plt.title("Relative abundance of {}".format(anlysis_type))
            else:
                plt.title("Relative abundance of {} of the {} sample".format(anlysis_type, sample_name))
            plt.ylabel("Relative abundance")
            plt.ylim(0, 1)
            ax.legend(loc="upper right", framealpha=0.4)
            if savepath is not None:
                plt.savefig('{}/{}{}.png'.format(
                    savepath, plot_name, pltnum + 1))
            plt.show()
            break


# SET THESE VARIABLES!
binning_folder = ""  # Set path to folder containing the bins
TPM_quant = ""  # Set path to file containing Salmon quantifiaction data
display_sample_name = None  # Edit this to the name of the sample

# Opening binning and quantification files
com_list_notdiv = read_bin_file(binning_folder)
quant_data = get_tpmfromfile(TPM_quant)
quant_data["TPM_log2"] = np.log2(quant_data["TPM"] + 1)


# Generating different ways of calculating the abundance.
# The keys are also used as labels for the bar plots.
contigs = {}
for num, contig in enumerate(com_list_notdiv):
    samples_in_bin = quant_data.loc[quant_data["Name"].isin(contig.keys())]
    contigs[num] = {}
    contigs[num]["TPM"] = place_in_dict(list(samples_in_bin["TPM"]))
    contigs[num]["TPM_log2"] = place_in_dict(list(samples_in_bin["TPM_log2"]))

# Calculate all sums to generate relative
sum_of_values = {}
for single in contigs[0].keys():
    sum_of_values[single] = {"mean": sum([contigs[i][single]["mean"] for i in contigs.keys()]),
                             "median": sum([contigs[i][single]["median"] for i in contigs.keys()]),
                             "sum": sum([contigs[i][single]["sum"] for i in contigs.keys()]),
                             }

# Get relative values of the lists
for single in contigs.keys():
    for count_type in contigs[single].keys():
        for type_an in contigs[single][count_type]:
            contigs[single][count_type][type_an] = contigs[single][count_type][type_an] / sum_of_values[count_type][type_an]

# Changing the form the dictionary to allow for the making of stacked bar plots
# new_dict = {TPM: {mean, median, sum}}
abundance_dict = {}
for single in contigs[0].keys():
    abundance_dict[single] = {}
    for an_type in contigs[0][single].keys():
        abundance_dict[single][an_type] = [contigs[i][single][an_type] for i in contigs.keys()]

# Adding the label of the true bins to list
real_bins = []
for single in com_list_notdiv:
    real_bins.append(mode([i[0:2] for i in single]))


create_bars(abundance_dict, sample_name=display_sample_name)

