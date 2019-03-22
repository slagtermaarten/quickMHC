# QuickMHC

Wrapper around NetMHCpan3.0 and a C++ implementation of our peptide similarity to self 
determination code for running on *nix systems.  Results are stored in a PostgreSQL 
database for future lookups.  The $USER is assumed to have read/write access to this 
database, named `binding_affinity`.

Binding affinity predictions are stored as such (example from the `A0201` table)

    peptide     peptide_score_log50k   affinity   percentile_rank
    ELVISLIVE   2.48e+04               6.48e-02   29
    SYFPEITHI   1.55e-01               9.32e+03   13

STS predictions are stored as such (example from the `STS_A1101_1.9` table)
 
    peptide     different_from_self
    SYFPEITHI   TRUE
    ELVISLIVE   FALSE


# Installation

* Install PostgreSQL on your system, initialize a database and couple a database user to 
  your Linux user account

* Configure `quickMHC` to know what user and database to use. You can do this by defining 
  any config file (found using the regex: `.*quickMHC.*\\.yaml`) in the following 
  directories (ordered in descending priority): ~/.config, the current working directory, the location of quickMHC's installation. 

# Usage

* Load quickMHC in R session
* Create quickMHC object and query peptides of interest
    

    ```
    ## For binding affinity predictions
    library(quickMHC)
    BA_predictor <- BindingPredictor$new(hla_allele = 'A0201')
    ## Peptide will be computed the first time around, will be looked up from database in 
    ## subsequent queries
    BA_predictor$query(c('SYFPEITHI'))
    # peptide     peptide_score_log50k   affinity   percentile_rank
    # SYFPEITHI   1.55e-01               9.32e+03   13

    ## For STS predictions
    STS_predictor <- STSPredictor$new(hla_allele = 'A0201')
    BA_predictor$query(c('SYFHPETHI'))
    # peptide     different_from_self
    # ELVISLIVE   FALSE
    ```


# TODO

* Make location of NetMHCpan configurable
* Add functionality to create self lists for STS
