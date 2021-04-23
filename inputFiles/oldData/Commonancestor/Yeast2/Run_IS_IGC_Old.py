#! /usr/bin/python3

from IGCexpansion.CodonGeneconv import ReCodonGeneconv
import argparse, os
import numpy as np

def check_folder(folder_name):
    # if the folder doesn't exist, 
    # this function creats it by assuming the outer folder exists (potential bug...)
    if not os.path.isdir(folder_name):
        os.mkdir(folder_name)

if __name__ == '__main__':
    paralog = ['YAL056W', 'YOR371C']
    Force = None
    alignment_file = './YAL056W_YOR371C_input.fasta'
    newicktree = './YeastTestTree.newick'
    Force = None
    model = 'HKY'  # choose from 'HKY' and 'MG94'
    save_folder = './save/'
    check_folder(save_folder)
    save_name = save_folder + model + 'Yeastori.txt'

    summary_folder = './summary/'
    check_folder(summary_folder)

    test = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = model, Force = Force, clock = None, save_path = './save/', save_name = save_name,
                            tau=0.1, omega=0.1, kappa=1.0, inibl=0.1)
    test.get_mle()
    test.get_ExpectedNumGeneconv()

    test.get_individual_summary(summary_path='./save/')
