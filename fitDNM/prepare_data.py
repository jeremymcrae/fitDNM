#' identify exclusion rows
#'
#' @param data dataframe of either mutation rates per site, or severity scores
#'        per site.
#'
#' @export
#'
#' @return vector of booleans, where TRUE indicates a row to exclude
get_rows_to_exclude <- function(data) {
    
    # find rows with NA values so we can exclude them later
    bases = c('A', 'C', 'G', 'T')
    exclude = data[, bases] < 0 | is.na(data[, bases])
    exclude = apply(exclude, 1, any)
    
    return(exclude)
}

#' prepare dataframe for analysis
#'
#' This changes some columns names to lower case, and excludes some predefined
#' rows. We want to exclude rows without severity scores. This function is then
#' applied for both the severity scores and mutation rates dataframes, so as to
#' keep the dataframes consistent.
#'
#' @param data dataframe of either mutation rates per site, or severity scores
#'        per site.
#' @param exclude vector of booleans, where TRUE indicates a row to exclude
#'
#' @export
#'
#' @return dataframe with  columns renamed, and unecessary rows excluded
tidy_table <- function(data, exclude) {
    
    cols = c('Gene', 'Chr', 'Ref', 'Pos')
    
    new_names = c()
    for (name in names(data)) {
        if (name %in% cols) { name = tolower(name) }
        new_names = c(new_names, name)
    }
    
    names(data) = new_names
    
    return(data[!exclude, ])
}
