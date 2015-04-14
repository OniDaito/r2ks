/**
 * A remake of the R2KS code, optimised for the apocrita cluster
 */


#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <list>
#include <cstdlib>
#include <getopt.h>

using namespace std;


// Options Struct
struct Options {
  std::string filename;

  // Internal
  unsigned int num_genes;
  unsigned int num_lists;
  unsigned int pivot;

};


/**
 * Calculate the weight for this term
 * For now, we are ignoring weight
 */

float calculateWeight(unsigned int idx, unsigned int pivot) {
  return 1.0;
  int h = (pivot - idx);
  float w = h * (h + 1.0) / 2.0;
  w = w < 0 ? 1 : w;
  return w;
}


/**
 * Perform the r2ks score for two lists
 * Take the two lists, the corresponding weights and the total weight, then perform
 * the scoring.
 */

float scoreLists(Options & options, std::vector<unsigned int> & gene_list0, std::vector<unsigned int> & gene_list1) {

  // First, calculate total weeight
  // TODO - we can so cache this

  // TODO - Given how weights are calculated this formula might work
  // totalweight = (pivot * ( pivot + 1) * (pivot + 2)) / 6.0 + (list length - pivot);

  float total_weight = 0.0;
  for (unsigned int i = 0; i < gene_list0.size(); ++i ){
    total_weight += calculateWeight(i, options.pivot);
  }

  // Create two diagonals that hold our current results.
  // Init with zeros. We add +1 as we want a border of 0s across the top and bottom of this matrix.

  std::vector< vector<float> > rmatrix(2, vector<float>(options.num_genes + 1));

  for (unsigned int i = 0; i < options.num_genes + 1; ++i ){
    rmatrix[0][i] = rmatrix[1][i] = 0.0;
  }

  // Our nasty double loop Compare all the things!

  float rvalue = 0;

  for (unsigned int i = 1; i < options.num_genes + 1; ++i ){

    unsigned int ivalue = gene_list0[i-1];

    for (unsigned int j = 1; j < options.num_genes + 1; ++j ){
      
      unsigned int jvalue = gene_list1[j-1];

      if (ivalue == jvalue){
          float iw = calculateWeight(i, options.pivot);
          float jw = calculateWeight(j, options.pivot);
          float w = iw < jw ? iw : jw;
          rmatrix[1][j] = rmatrix[0][j-1] + w;
          
      } else {

        rmatrix[1][j] = rmatrix[0][j] +  rmatrix[1][j-1] - rmatrix[0][j-1];
        
      }
      // Now check for the largest R score
      float second_term =  static_cast<float>(i * j) / (options.num_genes * options.num_genes);
      float nvalue = (rmatrix[1][j] / total_weight) - second_term; 

      if (rvalue < nvalue) {
        rvalue = nvalue;
      }
    }


    // Flip over the memory in the rmatrix (front row becomes rear row)
    for (unsigned int k = 0; k < options.num_genes + 1; ++k ){
      rmatrix[0][k] = rmatrix[1][k];
    }
  }

  return rvalue * sqrt(options.num_genes);
}


/**
 * Read the header block from our file
 */

void readHeaderBlock(Options & options) {
  std::ifstream is (options.filename.c_str());
  is >> options.num_genes >> options.num_lists;
  is.close();
}


/**
 * Read a line from a file.
 * Assume gene_list is pre-populated with zeros
 */

void readLineIndex(Options & options, int idx, std::vector<unsigned int> & gene_list ) {  

  std::ifstream fin;
  fin.open(options.filename.c_str());
  fin.seekg(std::ios::beg);

  for(int i=0; i < idx; ++i){
    fin.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
  }

  // The list we are reading in is a set of gene expressions.
  // The position/index of this number is the gene itself.
  // We want a list that ranks these indices

  unsigned int gene_expression;
  unsigned int gidx = 0;
  while (fin >> gene_expression && gidx < options.num_genes) {
    gene_list[gene_expression] = gidx;
    gidx++;
  }
}


/**
 * parse our command line options
 */

void parseCommandOptions (int argc, const char * argv[], Options &options) {
  int c;
  int digit_optind = 0;
  static struct option long_options[] = {
  };

  int option_index = 0;

  while ((c = getopt_long(argc, (char **)argv, "f:?", long_options, &option_index)) != -1) {
    int this_option_optind = optind ? optind : 1;
    switch (c) {
      case 0 :
        break;
  
      case 'f' :
        options.filename = std::string(optarg);
        break;

      default:
        std::cout << "?? getopt returned character code" << c << std::endl;
      }
  }
  
  if (optind < argc) {
      std::cout << "non-option ARGV-elements: " << std::endl;
      while (optind < argc)
          std::cout << argv[optind++];
      std::cout << std::endl;
  }

}

/**
 * Main entry point
 */

int main (int argc, const char * argv[]) {
  Options ops;

  // Defaults for Options
  ops.pivot = 0;

  parseCommandOptions(argc,argv,ops);

  readHeaderBlock(ops);

  std::vector<unsigned int> gene_list0(ops.num_genes);
  std::vector<unsigned int> gene_list1(ops.num_genes);

  readLineIndex(ops, 1, gene_list0);
  readLineIndex(ops, 2, gene_list1);

  float rvalue = scoreLists(ops, gene_list0, gene_list1 );

  std::cout << rvalue << std::endl;

}