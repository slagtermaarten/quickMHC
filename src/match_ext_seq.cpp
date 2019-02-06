#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

bool match_seq_ext_single(const string &s1,
                          const string &s2,
                          const Rcpp::NumericMatrix &scorematrix,
                          const string &dimnames,
                          const bool debug = false) {
  int s1_L = s1.size(), s2_L = s2.size();

  // if (debug) {
  //   for (int i=0; i < scorematrix.nrow(); i++) {
  //     for (int j=0; j < scorematrix.ncol(); j++) {
  //       cout << scorematrix(i,j) << ' ';
  //     }
  //     cout << endl;
  //   }
  // }

  // Input checking
  // if (false && s1_L != s2_L) {
  //   throw range_error("Input peps of unequal lengths");
  // }
  // if (false && s1_L != 9) {
  //   throw range_error("Input pep 1 not 9 aas longs");
  // }

  Rcpp::NumericVector scorevec(s1_L);
  Rcpp::LogicalVector pmatches(s1_L);

  int r_idx, c_idx;
  for (int i=0; i<s1_L; i++) {
    r_idx = dimnames.find(s1[i]);
    c_idx = dimnames.find(s2[i]);
    scorevec[i] = scorematrix(r_idx, c_idx);
    pmatches[i] = (s1[i] == s2[i]);
    if (debug) {
      Rcout << r_idx << dimnames[r_idx] << ' ' << 
               c_idx <<  dimnames[c_idx] << ' ' << 
               // scorematrix[r_idx, c_idx] << ' ' << 
               scorevec[i] << ' ' << 
               pmatches[i] << endl;
    }
  }

  if (!pmatches(4)) return true;

  // Substitution matrix score
  int SM_score = 0;
  // Positional changes regardless of biophysical effect
  int N_muts = 0;
  // Sum positions 3 til 8
  for (int i=2; i<=7; i++) {
    // cout << i << endl;
    SM_score = SM_score + scorevec[i];
    N_muts = N_muts + abs(1-pmatches[i]);
  }
  if (N_muts >= 3 || SM_score <= 5) 
    return true;

  // Sum positions 6 til 8
  int N_muts_right = 0;
  for (int i=5; i<=7; i++) {
    N_muts_right = N_muts_right + abs(1-pmatches[i]);
  }
  if (N_muts >= 2 && N_muts_right >= 2) 
    return true;

  // Sum positions 3 and 4
  int N_muts_left = 0;
  for (int i=2; i<=3; i++) {
    N_muts_left = N_muts_left + abs(1-pmatches[i]);
  }
  if (N_muts_left == 2) 
    return true;

  return false;
}

// [[Rcpp::export]]
Rcpp::LogicalVector match_seq_ext_Rcpp(const Rcpp::CharacterVector &query_peps,
    const Rcpp::CharacterVector &ref_list,
    const Rcpp::NumericMatrix &scorematrix,
    const bool debug = false) {
  // Amount of query peps
  int N_peps = query_peps.size();
  // Amount of ref peptides
  int N_ref = ref_list.size();

  // Get colnames from matrix, assume they're equal to rownames
  // Convert from array to string
  Rcpp::List dimnames_both = scorematrix.attr("dimnames");
  vector< string > dimnames_vec = Rcpp::as< vector<string> >(dimnames_both[0]);
  string dimnames;
  for (int i = 0; i < dimnames_vec.size(); i++) {
    dimnames += dimnames_vec[i];
  }

  // Instantiate results vector different from self; not guilty until proven
  // otherwise
  Rcpp::LogicalVector DFS(N_peps);
  int i = 0;
  for (int j = 0; j < N_peps; j++) {
    string query_seq_str = as< string >(query_peps[j]);
    // query_seq_str <- query_peps[j];
    DFS(j) = true;
    i = 0;
    // Stop when a similar peptide has been encountered in the ref_list
    while (DFS(j) == true && i < N_ref) {
      string ref_seq = as< string >(ref_list[i]);
      // ref_seq = ref_list(i);
      DFS(j) = match_seq_ext_single(query_seq_str, ref_seq, scorematrix, 
          dimnames, debug);
      i++;
    }
  }
  return DFS;
}
