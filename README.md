# QuickMHC

Wrapper around NetMHCpan3.0 and a C++ implementation of our peptide similarity to self 
determination code for running on *nix systems.  Results are stored in a PostgreSQL 
database for future lookups.  The $USER is assumed to have read/write access to this 
database, named `binding_affinity`.  I can make this configurable upon request.

INSERT USAGE EXAMPLES

Binding affinity predictions are stored as such (example from the `A0201` table)

    peptide     peptide_score_log50k   affinity   percentile_rank
    ELVISLIVE   2.48e+04               6.48e-02   29
    SYFPEITHI   1.55e-01               9.32e+03   13

STS predictions are stored as such (example from the `STS_A1101_1.9` table)
 
    peptide     different_from_self
    SYFPEITHI   TRUE
    ELVISLIVE   FALSE
