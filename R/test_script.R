# Test R Script

run_tests <- function() {

database_loc <- normalizePath("data/vdem_data.sqlite")

select_vars <- c('e_migdppcln','year_factor','country_name')
varnames <- paste0('V',2:900)

output <- run_vdem(varnames = varnames,select_vars=select_vars,full_formula = as.formula(v2x_polyarchy~e_migdppcln + year_factor + country_name),
                  dbcon=database_loc,num_iters=5)
return(output)
}
