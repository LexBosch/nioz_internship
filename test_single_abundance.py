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
    return {"mean": mean(data_list)}


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
                ax.bar(labels, data_layer, 0.35, bottom=layer_height)
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


def create_dict_contigs(name: str, binning_folder: str, quant_file: str) -> list:
    """ Prepares data for the creation of plots

    :param name: Name of the plot
    :param binning_folder: Folder containing the bins
    :param quant_file: Quantification file data
    :return: list containing sorted quantification data
    """
    bin_full = read_bin_file(binning_folder)
    quant_data = get_tpmfromfile(quant_file)
    quant_data["TPM_log2"] = np.log2(quant_data["TPM"] + 1)
    if name == "log_tpm":
        quant_data["TPM_log2"] = np.log2(quant_data["TPM"] + 1)
        output_list = (calculate_realtive(
            [mean(list(quant_data.loc[quant_data["Name"].isin(contig.keys())]["TPM_log2"])) for contig in bin_full]))
    elif name == "tpm":
        quant_data["TPM"] = (quant_data["TPM"] + 1)
        output_list = (calculate_realtive(
            [mean(list(quant_data.loc[quant_data["Name"].isin(contig.keys())]["TPM"])) for contig in bin_full]))
    output_list.sort()
    return output_list

# SET THESE VARIABLES
bins_folder = ""  # Path to folder containing bins
quant_path = ""  # Path to quant file created by samtools


full_dict = {}
full_dict["log_tpm_even"] = create_dict_contigs("log_tpm", bins_folder, quant_path)

plz = (full_dict.values())
plz2 = max([len(i) for i in full_dict.values()])

for i in full_dict.keys():
    while len(full_dict[i]) < plz2:
        full_dict[i].append(0)
bar_value_dict = {"TPM_log": full_dict}
create_bars(bar_value_dict)


