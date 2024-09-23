# title: "Helper functions"
# author: Guoqiang Zhang
# email: guoqiang.zhang@ki.se


##### Function to Calculate Proportions of Sexual Identities after Imputation #####

# among demographic subgroups
calc_prop_imp_subgroup <- function( implist, design, sexual_identities, demog_var, year ) {
  
  # calculate sample size
  sample_size <- lapply( implist, function( df ) {
    table( df[[ demog_var ]] ) 
  } )
  combined_freqs <- as.data.frame( Reduce( "+", sample_size ) / length( implist ) )
  colnames( combined_freqs ) <- c( "subgroup", "sample_size" )
  combined_freqs$sample_size <- prettyNum( 
    round( combined_freqs$sample_size, 0 ), big.mark = ",", preserve.width = "none" )
  names( combined_freqs )[ names( combined_freqs ) == "sample_size" ] <- paste0( "sample_size_", year )
  
  # calculate proportions
  combined_results_list <- list()
  
  for ( identity in sexual_identities ) {
    combined_summary <- summary( MIcombine( 
      with( design,
            svyby( formula = as.formula( paste0( "~ I( sexual_identity_", year, " == '", identity, "')" ) ),
                   by = as.formula( paste0( "~", demog_var ) ),
                   FUN = svyciprop,
                   method = "beta" ) ) ) )
    
    results_df <- rownames_to_column( combined_summary[ , c( "results", "(lower", "upper)" ) ], var = "subgroup" )
    colnames( results_df ) <- c( "subgroup", paste0( identity, "_point_estimate_", year ), paste0( identity, "_lower_ci_", year ), paste0( identity, "_upper_ci_", year ) )
    
    combined_results_list[[ identity ]] <- results_df
  }
  
  final_combined_df <- Reduce( function( x, y ) {
    merge( x, y, by = "subgroup" ) 
  }, 
  combined_results_list )
  
  final_combined_df <- merge( final_combined_df, combined_freqs, by = "subgroup" )
  
  return( final_combined_df )
}


# in Stockholm County
calc_prop_imp_overall <- function( design, sexual_identities, year ) {
  
  # calculate proportions
  combined_results_list <- list()
  
  for ( identity in sexual_identities ) {
    combined_summary <- summary( MIcombine( 
      with( design,
            svyciprop( formula = as.formula( paste0( "~ I( sexual_identity_", year, " == '", identity, "')" ) ),
                       method = "beta" ) ) ) )
    
    results_df <- rownames_to_column( combined_summary[ , c( "results", "(lower", "upper)" ) ], var = "subgroup" )
    colnames( results_df ) <- c( "subgroup", paste0( identity, "_point_estimate_", year ), paste0( identity, "_lower_ci_", year ), paste0( identity, "_upper_ci_", year ) )
    results_df[ 1, 1 ] <- "Stockholm County"
    
    combined_results_list[[ identity ]] <- results_df
  }
  
  final_combined_df <- Reduce( function( x, y ) {
    merge( x, y, by = "subgroup" ) 
  }, 
  combined_results_list )
  
  final_combined_df$sample_size <- prettyNum( 
    nrow( get( paste0( "d_", year ) ) ), big.mark = ",", preserve.width = "none" )
  names( final_combined_df )[ names( final_combined_df ) == "sample_size" ] <- paste0( "sample_size_", year )
  
  return( final_combined_df )
}


##### Function to Calculate Proportion of Change in Sexual Identity after Imputation #####

# in Stockholm County
calc_prop_fluidity_imp_overall <- function( design, year ) {
  
  combined_summary <- summary( MIcombine( 
    with( design,
          svyciprop( formula = ~ I( sexual_identity_fluidity_cat == "changed" ),
                     method = "beta" ) ) ) )
  
  results_df <- rownames_to_column( combined_summary[ , c( "results", "(lower", "upper)" ) ], var = "subgroup" )
  colnames( results_df ) <- c( "subgroup", paste0( "changed_point_estimate_", year ), paste0( "changed_lower_ci_", year ), paste0( "changed_upper_ci_", year ) )
  results_df[ 1, 1 ] <- "Stockholm County"
  
  results_df$sample_size <- prettyNum( 
    nrow( get( paste0( "d_", year ) ) ), big.mark = ",", preserve.width = "none" )
  names( results_df )[ names( results_df ) == "sample_size" ] <- paste0( "sample_size_", year )
  
  return( results_df )
}