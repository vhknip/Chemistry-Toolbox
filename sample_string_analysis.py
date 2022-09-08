
import pandas as pd


def common_prefix(seq_samples: list) -> str:
    """function takes in list of strings and returns longest common prefix"""
    prefix = ''
    first_word = seq_samples[0]
    for char in range(len(first_word)):
        for sample in seq_samples:
            if char == len(sample) or sample[char] != first_word[char]:
                return prefix
        prefix += first_word[char]
    return prefix



        
def string_rep_convert(rep_string: str) -> list:
    """this function converts the string representaiton of a list to a list"""
    unbracketed = rep_string.strip('][')
    converted = unbracketed.split(',')
    return converted

def sample_group(df: pd.DataFrame(),
                 attr: str,
                 attr_groups: dict,
                 sample_col: str='Sample',) -> pd.DataFrame():
    """
    

    Parameters
    ----------
    df : pandas.DataFrame
        dataframe containing sample names and analysis information
    attr : str
        name of attribute of information extracted from sample name
    attr_groups: dict
        dictionary of all possible attributes and abbreviations or 
        representations of these attributes that can be found as a substring
        in the sample name
        Example: {'apple':['apple','app','apl'],'mango'['mango','mang']}
    sample_col: str
        Name of the sample column in the pandas dataframe.
        Values of this column are analyzed for substrings to assign attributes
        The default value of this function is 'Sample'
        
    

    Returns
    -------
    df :
        modification of original dataframe with grouping based on substring 
        analysis and assignment of values to specified sample attribute

    """
    
    df.insert(1, attr, ['CHECK LABEL'] * len(df))
    group_vals = list(attr_groups.values())
    group_labels = list(attr_groups.keys())
    for abbrevs_list in group_vals:
        for rows in range(len(df)):
            sample_name = df.loc[rows, sample_col]
            if any(abbrev.lower() in sample_name.lower() for abbrev in abbrevs_list):
                df.loc[rows,attr] = group_labels[group_vals.index(abbrevs_list)]
    return df