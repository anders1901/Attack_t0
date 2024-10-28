import os

import copy as py_copy

from datetime import timedelta

import re 

import numpy as np 

import subprocess

import Dilithium_functions as dilithium


def run_my_process(*args):
    arguments = list(args)
    process = subprocess.Popen(arguments, 
                               stdout=subprocess.PIPE,
                               universal_newlines=True)
    process.wait()
    while True:
        output = process.stdout.readlines()
        return output

    
def print_np_file_infos(npz_file_name):
    """
    This function reads an .npz file name prints its infos and return the corresponding np object 
                
    Parameters
    ----------
    npz_file_name                (str): Input lp file name.

    Returns
    ----------
    loaded_array (np object of arrays): Object containing the arrays in the npz_file_name
    """
    
    loaded_array = np.load(npz_file_name)
    
    print(f"loading: {npz_file_name}")
    for element in loaded_array.files:
        print(f"  >>> '{element}': {loaded_array[element].shape}")
    return loaded_array

    
def open_keys(nb_keys = 1, include_sk = True, keys_file_name = f"dilithium/ref/PQCsignKAT_Dilithium{dilithium.MODE}.rsp"):
    """
    This function reads nb_keys dilithium keys from the file keys_file_name. 
                
    Parameters
    ----------
    nb_keys           (int) : Number of keys to read from the keys file.
    include_sk        (bool): If the sk is included in the file 
    keys_file_name    (str) : Name of the file containing the keys to read. 

    Returns
    ----------
    dict_pks (dict[int]:str): At the index i contains the public key i. 
    dict_sks (dict[int]:str): At the index i contains the secret key i. (if include_sk=True)
    """
    if not isinstance(nb_keys, int):
        raise ValueError(f"nb_keys must be an int between 1 and 100, not {nb_keys}")
        
    if nb_keys < 0:
        raise ValueError(f"nb_keys must be an int between 1 and 100, not {nb_keys}")
        
    count_pattern = r"count\s*=\s*(\d+)"
    pk_pattern = r"pk\s*=\s*([\da-fA-F]*)"
    if include_sk:
        sk_pattern = r"sk\s*=\s*([\da-fA-F]*)"
        
    with open(f"{os.getcwd()}/{keys_file_name}", "r") as file:
        # We discard the first flag for the version of Dilithium as well as the first and last \n
        lines = "".join(file.readlines())
    
    counts = re.findall(count_pattern, lines)
    pks    = re.findall(pk_pattern, lines)
    if include_sk:
        sks = re.findall(sk_pattern, lines) 
        
    dict_pks = {int(counts[key_index]) : pks[key_index] for key_index in range(nb_keys)}
    if include_sk:
        dict_sks = {int(counts[key_index]) : sks[key_index] for key_index in range(nb_keys)}
        return dict_pks, dict_sks
    
    return dict_pks


def Antt2Aintt(A):
    """
    This function converts a dilithium matrix A in the NTT 
    domain to its representation in the normal domain. 
  
    Parameters
    ----------
    A (K-D array): Dilithium matrix A of shape K x L x N elements, in NTT domain.

    Returns
    ----------
    A (K-D array): Dilithium matrix A of shape K x L x N elements, in normal domain.
    """ 
    A_intt_ = []
    for i in range(dilithium.K):
        A_k = []
        for j in range(dilithium.L):
            a_test = py_copy.deepcopy(A[i][j])
            dilithium.poly_reduce(a_test)
            dilithium.invntt_frominvmont(a_test)
            A_k.append(a_test)
        A_intt_.append(A_k)
    return [[[dilithium.montgomery_reduce(a) for a in al] for al in ak]for ak in A_intt_]


def format_results(reults, include_sk = True):
    """
    This function formats the results of the uncompression of t0. 
  
    Parameters
    ----------
    results                    list(str): List containing the resultst for each radius
    include_sk                      bool: If set to true, suppose we know the real t0 and compare
    
    Returns
    ----------
    (hours, minutes, seconds) tuple(int): Time taken for the attack
    total_signs                      int: Number of signatures used 
    mkdwn_table                 Markdown: Table that sums up the attack
    """
    if include_sk:
        mkdwn_table = f'''| Round | Time         | C   |# Signatures | # Inequalities0  | # Inequalities1  |  Min $\infty$ norm | Max $\infty$ norm |  \n|---------|-----------------------|--------------|----------------------|----------------------|----------------------|----------------------|----------------------|  \n'''
    else:
        mkdwn_table = f'''| Round | Time         | C   |# Signatures | # Inequalities0  | # Inequalities1  |  \n|---------|-----------------------|--------------|----------------------|----------------------|----------------------|  \n'''

    total_time   = timedelta() 
    total_signs  = 0
    total_ineq   = 0
    total_ineq0  = 0
    total_ineq1  = 0
    
    for round_number, output_ in enumerate(reults):
        val_test = float(output_[1].split(" sec/")[0].split(": ")[1])
        hours_, remainder = divmod(val_test, 3600)
        minutes_, seconds_ = divmod(remainder, 60)  
        
        hours_ = int(hours_)
        minutes_ = int(minutes_)
        seconds_ = int(seconds_)
        
        total_time  += timedelta(hours = hours_, minutes = minutes_, seconds = seconds_)

        nb_ineq   = int(output_[7].split(": ")[1][:-1])
        total_ineq += nb_ineq
        
        nb_ineq0   = max(eval("[" + output_[1].split("/")[2][:-9] + "]"))
        total_ineq0 += nb_ineq0

        nb_ineq1   = max(eval("[" + output_[1].split("/")[3][:-9] + "]"))
        total_ineq1 += nb_ineq1

        if include_sk:
            min_error = float(output_[-2].split(": ")[1][:-1])
            max_error = float(output_[-1].split(": ")[1][:-1])
        
        signs_       = int(output_[1].split("/")[1].split(" ")[0])
        total_signs = max(total_signs, signs_)
        
        if include_sk:
            mkdwn_table += f'''|{round_number+1}     | {hours_}h{minutes_}m{seconds_}s  |   {float(output_[0].split(":")[1][1:])} |   {signs_}      |  {nb_ineq0}       |  {nb_ineq1}       | {min_error} | {max_error} |  \n'''
        else:
            mkdwn_table += f'''|{round_number+1}     | {hours_}h{minutes_}m{seconds_}s  |   {float(output_[0].split(":")[1][1:])} |   {signs_}      |  {nb_ineq0}       |  {nb_ineq1}       |  \n'''
    hours, remainder = divmod(total_time.seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    
    if include_sk:
        mkdwn_table += f'''|Total     | {hours}h{minutes}m{seconds}s  |   - |   {total_signs}      |  {total_ineq0}       |  {total_ineq1}       | - | - |  \n'''
    else:
        mkdwn_table += f'''|Total     | {hours}h{minutes}m{seconds}s  |   - |   {total_signs}      |  {total_ineq0}       |  {total_ineq1}       |  \n'''
    return (hours, minutes, seconds), total_signs, mkdwn_table