#include <Rcpp.h>
#include <string>
#include <map>
#include <regex>
#include <vector>
#include <algorithm>
#include <string>
using namespace Rcpp;

// [[Rcpp::plugins("cpp11")]]



// [[Rcpp::export]]
NumericVector get_features(std::string rna_seq, CharacterVector Struct_seq){

  // Features Vector
  NumericVector features(514);

  // Hash table for frequencies of tri-motifs
  std::map<std::string, double> probabilities;

  // Vector of motifs ordered according to what we want
  std::vector<std::string> sequences = {
    "aaa000","aaa001","aaa010","aaa011","aaa100","aaa101","aaa110","aaa111","aau000","aau001",
    "aau010","aau011","aau100","aau101","aau110","aau111","aac000","aac001","aac010","aac011",
    "aac100","aac101","aac110","aac111","aag000","aag001","aag010","aag011","aag100","aag101",
    "aag110","aag111","aua000","aua001","aua010","aua011","aua100","aua101","aua110","aua111",
    "auu000","auu001","auu010","auu011","auu100","auu101","auu110","auu111","auc000","auc001",
    "auc010","auc011","auc100","auc101","auc110","auc111","aug000","aug001","aug010","aug011",
    "aug100","aug101","aug110","aug111","aca000","aca001","aca010","aca011","aca100","aca101",
    "aca110","aca111","acu000","acu001","acu010","acu011","acu100","acu101","acu110","acu111",
    "acc000","acc001","acc010","acc011","acc100","acc101","acc110","acc111","acg000","acg001",
    "acg010","acg011","acg100","acg101","acg110","acg111","aga000","aga001","aga010","aga011",
    "aga100","aga101","aga110","aga111","agu000","agu001","agu010","agu011","agu100","agu101",
    "agu110","agu111","agc000","agc001","agc010","agc011","agc100","agc101","agc110","agc111",
    "agg000","agg001","agg010","agg011","agg100","agg101","agg110","agg111","uaa000","uaa001",
    "uaa010","uaa011","uaa100","uaa101","uaa110","uaa111","uau000","uau001","uau010","uau011",
    "uau100","uau101","uau110","uau111","uac000","uac001","uac010","uac011","uac100","uac101",
    "uac110","uac111","uag000","uag001","uag010","uag011","uag100","uag101","uag110","uag111",
    "uua000","uua001","uua010","uua011","uua100","uua101","uua110","uua111","uuu000","uuu001",
    "uuu010","uuu011","uuu100","uuu101","uuu110","uuu111","uuc000","uuc001","uuc010","uuc011",
    "uuc100","uuc101","uuc110","uuc111","uug000","uug001","uug010","uug011","uug100","uug101",
    "uug110","uug111","uca000","uca001","uca010","uca011","uca100","uca101","uca110","uca111",
    "ucu000","ucu001","ucu010","ucu011","ucu100","ucu101","ucu110","ucu111","ucc000","ucc001",
    "ucc010","ucc011","ucc100","ucc101","ucc110","ucc111","ucg000","ucg001","ucg010","ucg011",
    "ucg100","ucg101","ucg110","ucg111","uga000","uga001","uga010","uga011","uga100","uga101",
    "uga110","uga111","ugu000","ugu001","ugu010","ugu011","ugu100","ugu101","ugu110","ugu111",
    "ugc000","ugc001","ugc010","ugc011","ugc100","ugc101","ugc110","ugc111","ugg000","ugg001",
    "ugg010","ugg011","ugg100","ugg101","ugg110","ugg111","caa000","caa001","caa010","caa011",
    "caa100","caa101","caa110","caa111","cau000","cau001","cau010","cau011","cau100","cau101",
    "cau110","cau111","cac000","cac001","cac010","cac011","cac100","cac101","cac110","cac111",
    "cag000","cag001","cag010","cag011","cag100","cag101","cag110","cag111","cua000","cua001",
    "cua010","cua011","cua100","cua101","cua110","cua111","cuu000","cuu001","cuu010","cuu011",
    "cuu100","cuu101","cuu110","cuu111","cuc000","cuc001","cuc010","cuc011","cuc100","cuc101",
    "cuc110","cuc111","cug000","cug001","cug010","cug011","cug100","cug101","cug110","cug111",
    "cca000","cca001","cca010","cca011","cca100","cca101","cca110","cca111","ccu000","ccu001",
    "ccu010","ccu011","ccu100","ccu101","ccu110","ccu111","ccc000","ccc001","ccc010","ccc011",
    "ccc100","ccc101","ccc110","ccc111","ccg000","ccg001","ccg010","ccg011","ccg100","ccg101",
    "ccg110","ccg111","cga000","cga001","cga010","cga011","cga100","cga101","cga110","cga111",
    "cgu000","cgu001","cgu010","cgu011","cgu100","cgu101","cgu110","cgu111","cgc000","cgc001",
    "cgc010","cgc011","cgc100","cgc101","cgc110","cgc111","cgg000","cgg001","cgg010","cgg011",
    "cgg100","cgg101","cgg110","cgg111","gaa000","gaa001","gaa010","gaa011","gaa100","gaa101",
    "gaa110","gaa111","gau000","gau001","gau010","gau011","gau100","gau101","gau110","gau111",
    "gac000","gac001","gac010","gac011","gac100","gac101","gac110","gac111","gag000","gag001",
    "gag010","gag011","gag100","gag101","gag110","gag111","gua000","gua001","gua010","gua011",
    "gua100","gua101","gua110","gua111","guu000","guu001","guu010","guu011","guu100","guu101",
    "guu110","guu111","guc000","guc001","guc010","guc011","guc100","guc101","guc110","guc111",
    "gug000","gug001","gug010","gug011","gug100","gug101","gug110","gug111","gca000","gca001",
    "gca010","gca011","gca100","gca101","gca110","gca111","gcu000","gcu001","gcu010","gcu011",
    "gcu100","gcu101","gcu110","gcu111","gcc000","gcc001","gcc010","gcc011","gcc100","gcc101",
    "gcc110","gcc111","gcg000","gcg001","gcg010","gcg011","gcg100","gcg101","gcg110","gcg111",
    "gga000","gga001","gga010","gga011","gga100","gga101","gga110","gga111","ggu000","ggu001",
    "ggu010","ggu011","ggu100","ggu101","ggu110","ggu111","ggc000","ggc001","ggc010","ggc011",
    "ggc100","ggc101","ggc110","ggc111","ggg000","ggg001","ggg010","ggg011","ggg100","ggg101",
    "ggg110","ggg111","stem_prob", "poly_T"};

  // Hash table, initialize sequences to zero frequencies
  for (unsigned i = 0; i < sequences.size(); i++){
    probabilities[sequences[i]] = 0.0;
  }

  std::string window;
  std::regex pattern("((\\(+)\\.*(\\)+)\\.*$)");

  // Iterate through 1000 random samples of secondary structures
  for (unsigned i = 0; i < 1000; i++){

    std::string struct_seq = Rcpp::as<std::string>(Struct_seq[i]);

    for (unsigned j = 1; j < rna_seq.length() - 1 ; j++){

      window = rna_seq[j-1];
      window += rna_seq[j];
      window += rna_seq[j+1];

      char dot = '.';
      char one = '1';
      char zero = '0';

      if (struct_seq[j-1] == dot ){
        window += zero;
      } else {
        window += one;
      }

      if (struct_seq[j] == dot){
        window += zero;
      } else {
        window += one;
      }

      if (struct_seq[j+1] == dot){
        window += zero;
      } else {
        window += one;
      }

      probabilities[window] += 1;

    }

    // Check for stem loop
    std::smatch m;
    if (std::regex_search(struct_seq, m, pattern) && m[2].length() == m[3].length()){
      probabilities["stem_prob"] += 1;
    }

  }

  // rho independent terminator feature
  if (rna_seq.length() > 12){
    int llen = rna_seq.length();
    std::string proximal_part = rna_seq.substr(llen -13, 13);
    std::string distal_part = rna_seq.substr(llen - 7, 7);
    std::regex UVVUU("u(a|c|g)(a|c|g)uu");

    if (std::count(proximal_part.begin(), proximal_part.end(), 'u') > 2 && !std::regex_search(proximal_part, UVVUU) && std::count(distal_part.begin(), distal_part.end(), 'u') > 0){
      probabilities["poly_T"] += 1;
    }
  }

  for (unsigned i = 0; i < sequences.size(); i++){
    features[i] = probabilities[sequences[i]];
  }

  // compute probabilities
  double total = 0.0;
  for (int i = 0; i < 512; i++){
    total += features[i];
  }

  for(int i = 0; i < 512; i++){
    features[i] = features[i]/total;
  }

  features[512] = features[512]/1000.0;

  return(features);
}
