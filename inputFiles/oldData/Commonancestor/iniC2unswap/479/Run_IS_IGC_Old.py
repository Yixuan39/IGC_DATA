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
    paralog = ['01', '02']
    Force = None
    alignment_file = './Pillar479Swapped.fasta'
    newicktree = './Pillar479Swapped.newick'
    Force = None
    model = 'MG94'  # choose from 'HKY' and 'MG94'
    save_folder = './save/'
    check_folder(save_folder)
    save_name = save_folder + model + '479iniC2unswap.txt'

    summary_folder = './summary/'
    check_folder(summary_folder)

    test = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = model, Force = Force, clock = None, save_path = './save/', save_name = save_name,
                            tau=3.5, omega=0.01, kappa=3.0, inibl=0.1)
    test.get_mle()
    test.get_ExpectedNumGeneconv()

    test.get_individual_summary(summary_path='./save/')
