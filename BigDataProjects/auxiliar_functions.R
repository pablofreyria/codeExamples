#File to keep all functions that are used in repeatedly in the dataAnalysis file

library(dplyr)
library(reshape)
library(rlist)

# 
# Function to compute the number of claims by category of data_claims
num_claims_by_category <- function(df,category,drop_var = T){
  # Function to change dataframe `df` to a wide format on the category variable
  # Inputs:
  # -df: Data.frame with variables MemberID and Year whose unique combinations 
  #      are used as key to create new rows
  # -Category: Name of variable in Data frame whose unique values will provide 
  #           column names for the resulting dataframe. 
  # -drop_var : Logical to indicate if column with lowest variance is to be dropped to avoid colinearity
  
  # Outputs:
  # -cat_wide: Wide data frame format with the number of times category values appear
  #            for each MemberID and year
  
  f <- paste("MemberID + Year~",category,sep = "")
  
  cat_claims <- df %>%group_by_at(category) %>%count(MemberID, Year)
  cat_wide <- cast(cat_claims,as.formula(f),value = "n")
  cat_wide[is.na(cat_wide)] <- 0
  
  vars <- numeric(ncol(cat_wide)-2)
  for(j in c(3:ncol(cat_wide))){
    vars[(j-2)] <- var(cat_wide[,j],na.rm = T)
  }
  if(drop_var){
    print(colnames(cat_wide[-c(1,2)]))
    print(vars)
    ids <- which(vars==min(vars))[1]+2
    cat_wide <- cat_wide[,-ids]
  }
  
  return(cat_wide)
}

# Function to provide the mean summary for numerical variables specified in a name vector
aggregate_by_vec <- function(df,name_vec,fun = "mean"){
  # Function to create a summary of a numeric variable from repeated observations 
  # by unique MemberID and Year.
  # Inputs:
  # -df: Data.frame with variables MemberID and Year whose unique combinations 
  #      are used as key to create new rows
  # -name_vec: Name of variable in Data frame that is going to be sumarrized
  # -fun : Function to summarize the variable (e.g. mean, sum)
  
  # Outputs:
  # -data_sum: Data frame with the summarized value for memberID and Year
  
    
  
  f <- paste(name_vec[1],"~ MemberID + Year",sep = "")
  data_sum <- aggregate(as.formula(f), df, fun)
  
  if(length(name_vec)>1){
    for(i in c(2:length(name_vec))){
      
      f <- paste(name_vec[i],"~ MemberID + Year",sep = "")
      new_data <- aggregate(as.formula(f), df, fun)
      
      data_sum <- merge(data_sum,new_data,by = c("MemberID","Year"),all=T)
    }
  }
  
  return(data_sum)
}

# Function that, given a vector of column names, provides a list with wide tables
# with columns by each factor in variable name and value count of events
build_wide_format_table <- function(df,vec_categories,drop_var = T){
 # Function to return a list, each element is a wide format data frame as resulting 
 # from `num_claims_by_category`
 # Inputs:
 # -df: Data.frame with variables MemberID where variable names found in vec_categories 
 #      will be used to create a wide data frame format for each element.
 # -Vec_categories: Vector with variable names found in dataframe df to create 1
 #           wide data frame for each element in vector
 # Outputs:
 # -table_list: List with the generated tables
  
    table_list <- vector(mode = "list", length = length(vec_categories))
    i <- 1
    for(name in vec_categories){
      print_str <-sprintf("Run number %i,category %s",i,name)
      print(print_str)
      table_list[[i]] <- num_claims_by_category(df,name,drop_var)
      
      
      i <- i+1
    }
    names(table_list) <- vec_categories
    
    return(table_list)
  
}

#Function to change from list of data frames as produced by "build_wide_format_table"
# to a big data frame, merged by MemberID & Year
list_to_data_frame <- function(list_of_table){

  n <- length(list_of_table)
  df <- list_of_table[[1]]
  
  table_names <- names(list_of_table)
  colnames(df)[c(3:ncol(df))] <- paste(table_names[1],
                                       colnames(df)[c(3:ncol(df))],sep = ".")
  
  for (i in c(2:n)) {
    # add variable name to column
    n_cols <- ncol(list_of_table[[i]])
    colnames(list_of_table[[i]])[c(3:n_cols)] <- paste(table_names[i],
                                                  colnames(list_of_table[[i]])[c(3:n_cols)],
                                                  sep = ".")
    
    # create 1 gib data frame
    df <- merge(df,list_of_table[[i]],by=c("MemberID","Year"),all=T)
  }
  
  return(df)
  
}


# quick function for mse and rsme from a predicted vec and a data frame on "DIH"
get_rmse <- function(preds,df,varcol = "DIH"){
  
  mse <- mean((preds-df[,varcol])^2)
  rmse <- sqrt(mse)
  
  n_obs <- nrow(df)
  
  res <- data.frame("MSE"= mse,"RMSE" = rmse, "NObs" = n_obs)
  
  return(res)
}

# method to screen variables in linear models
screen_lm <- function(df,cutoff,varname = "DIH",FDR = T){
  # Returns linear model and vector of indices of significant variables
  # if FDR = T, the cutoff% of variavles with smallest p-values are selected
  # otherwise, variables with p-values less than cutoff% are selected
  f <- as.formula(paste(varname,"~."))
  
  lin_mod <- lm(f, data = df)
  anova_mod <- anova(lin_mod)
  var_pvals <- anova_mod$`Pr(>F)`[-ncol(df)] #remove analysis of residuals
  if(FDR){
    cutoff_val <- sort(var_pvals)[round(cutoff*ncol(data_train[,-1]),0)]
  }else{
    cutoff_val <- cutoff
  }
  
  var_ids <- (var_pvals<=cutoff_val)
  
  res <- list("Model" = lin_mod,"SelVars"=var_ids) 
  return(res)
}

# method to build the final tables with main and square terms after screening 
convert2screened <- function(df,id_lins,id_sqrs,id_rmve){
  
  new_data <- df
  new_data <- new_data[,id_lins]

  if(sum(id_sqrs)>0){
    sqr_data <- df[,-id_rmve]
    sqr_data <- sqr_data[,id_sqrs]^2
    colnames(sqr_data) <- paste(colnames(sqr_data),"^2",sep = "")
    
    new_data <- cbind.data.frame(new_data,sqr_data)  
  }
  
  
  return(new_data)
  
  
}
