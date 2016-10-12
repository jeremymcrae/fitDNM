def get_rows_to_exclude(data):
    ''' identify exclusion rows
    
    Args:
        data: dataframe of either mutation rates per site, or severity scores
            per site.
    
    Returns:
        vector of booleans, where TRUE indicates a row to exclude
    '''
    
    # find rows with NA values so we can exclude them later
    bases = ['A', 'C', 'G', 'T']
    exclude = (data[bases] < 0) | data[bases].isnull()
    exclude = exclude.apply(any, axis=1)
    
    return exclude

def tidy_table(data, exclude):
    ''' prepare dataframe for analysis
    
    This changes some columns names to lower case, and excludes some predefined
    rows. We want to exclude rows without severity scores. This function is then
    applied for both the severity scores and mutation rates dataframes, so as to
    keep the dataframes consistent.
    
    Args:
        data: dataframe of either mutation rates per site, or severity scores
            per site.
        exclude: vector of booleans, where TRUE indicates a row to exclude
    
    Returns:
        dataframe with  columns renamed, and unecessary rows excluded
    '''
    
    cols = ['Gene', 'Chr', 'Ref', 'Pos']
    
    new_names = []
    for name in data.columns:
        if name in cols:
            name = name.lower()
        new_names.append(name)
    
    data.columns = new_names
    
    return data[~exclude]
