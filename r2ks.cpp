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
#include <mpi.h>

using namespace std;


// Options Struct
struct Options {
  std::string filename;

  // Internal
  unsigned int num_genes;
  unsigned int num_lists;
  unsigned int pivot;

  // MPI
  int num_procs, mpi_id;

};


// MPI Datatype
typedef struct {
  int i,j;
  double result;
}MPIResult;

MPI_Datatype resultType;


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
 * Calculate the R Value for the given parameters
 */

double calculateR(unsigned int i, unsigned int j, double front, double total_weight, Options & options){
  double second_term = static_cast<double>((j+1) * (i+1)) / (options.num_genes * options.num_genes);
  return (front / total_weight - second_term);
}


/**
 * Perform the r2ks score for two lists
 * Take the two lists, the corresponding weights and the total weight, then perform
 * the scoring.
 */

double scoreLists(Options & options, std::vector<unsigned int> & gene_list0, std::vector<unsigned int> & gene_list1) {

  // First, calculate total weeight
  // TODO - we can so cache this

  // TODO - Given how weights are calculated this formula might work
  // totalweight = (pivot * ( pivot + 1) * (pivot + 2)) / 6.0 + (list length - pivot);

  float total_weight = 0.0;
  for (unsigned int i = 0; i < gene_list0.size(); ++i ){
    total_weight += calculateWeight(i, options.pivot);
  }

  std::vector<unsigned int> buff (gene_list0.size());

  for (unsigned int i = 0; i < gene_list0.size(); ++i ){
    buff[ gene_list1[i] ] = i;
  }

  // The algorithm follows the paper:
  // http://online.liebertpub.com/doi/pdf/10.1089/cmb.2012.0026
  // We optimise the algorithm by storing a history of previous values in the second loop
  // We take advantage of the fact that a gene is unique and always occurs in both lists and only once
  // thus rather than filling in the complete matrix, we only fill in the values we need by keeping
  // a record of where and what the previous values were.

  typedef struct {
    unsigned int pos_y;
    double value;
    double prev_value;
 
  }History;

  std::vector<History> history;

  double rvalue = 0.0;

  // First line has no history, so we do this one separately.
  {

    unsigned int ivalue = gene_list0[0];
    float pivot = buff[ivalue];

    double iw = calculateWeight(0, options.pivot);
    double jw = calculateWeight(pivot, options.pivot);
    double w = iw < jw ? iw : jw;
    double front = w;

    double prev = 0.0;

    History nh;
    nh.value = front;
    nh.pos_y = pivot;
    nh.prev_value = prev;

    history.push_back(nh);

  }

  for (History h : history){
    cout << h.pos_y << ",(" << h.value << "," << h.prev_value << ")" << endl;
  }
  cout << "---" << endl;


  // Now we have a history so we perform the operation normally

  for (unsigned int i = 1; i < options.num_genes; ++i ){

    unsigned int ivalue = gene_list0[i];
    float pivot = buff[ivalue];

    double iw = calculateWeight(i, options.pivot);
    double jw = calculateWeight(pivot, options.pivot);
    double w = iw < jw ? iw : jw;
    
    double front = 0.0;
    double prev = 0.0;

    std::vector<History> new_history;

    unsigned int prev_pos_y = 0;

    // If pivot is less than all the history that will bump things up
    if (pivot < history[0].pos_y){
      front = history[0].prev_value + w;

      // Add pivot as a new history item
      History nh;
      nh.value = front;
      nh.prev_value = prev;
      nh.pos_y = pivot;
      new_history.push_back(nh);
      prev_pos_y = pivot;
    }


    // Loop through the histories
   

    for (History h : history){
      prev = front;

      // Check if out pivot happened between previous and this history
      if (pivot > prev_pos_y && pivot < h.pos_y){
        // It did so insert
        front = h.prev_value + w;
        History nh;
        nh.prev_value = prev;
        nh.value = front;
        nh.pos_y = pivot;
        new_history.push_back(nh);
        prev = front;

      }
     
      front = prev + h.value - h.prev_value;
  
      // Create new history
      History nh;
      nh.prev_value = prev;
      nh.value = front;
      nh.pos_y = h.pos_y;
      new_history.push_back(nh);

      prev_pos_y = nh.pos_y;

    }
    
    prev = front;

    // Finally, we check to see if the pivot is still greater than all the histories
    History lh = history[history.size()-1];
    if (pivot > lh.pos_y){
      History nh;
      nh.prev_value = prev;
      nh.value = front + w;
      nh.pos_y = pivot;
      new_history.push_back(nh);
    }


    // Finally, swap histories
    history.swap(new_history);


    // Find the largest R Value
    // Do the R-Value test
    for (History h : history){
      double nvalue = calculateR(h.pos_y, i, h.value, total_weight, options);
      rvalue = nvalue > rvalue ? nvalue : rvalue; 
    }


    for (History h : history){
      cout << h.pos_y << ",(" << h.value << "," << h.prev_value << ")" << endl;
    }
    cout << "---" << endl;

  }


  // Create two diagonals that hold our current results.
  // Init with zeros. We add +1 as we want a border of 0s across the top and bottom of this matrix.

  /*std::vector<float> back (options.num_genes + 1);
  std::vector<float> front (options.num_genes + 1);


  for (unsigned int i = 0; i < options.num_genes + 1; ++i ){
    front[i] = back[i] = 0.0;
  }

  // Our nasty double loop Compare all the things!

  float rvalue = 0;

  for (unsigned int i = 1; i < options.num_genes + 1; ++i ){

    unsigned int ivalue = gene_list0[i-1];

    // We don't need to loop over all the inner list because genes are unique
    // so we zoom straight to the corresponding gene using our buff lookup
    // table. We compute the new front value and propagate forward while the
    // values

    float pivot = buff[ivalue] + 1;

    for (unsigned int j = 1; j < pivot; ++j){
      front[j] = back[j];
    }

    unsigned int jvalue = gene_list1[pivot-1];

    float iw = calculateWeight(i, options.pivot);
    float jw = calculateWeight(pivot, options.pivot);
    float w = iw < jw ? iw : jw;
    front[pivot] = back[pivot-1] + w;

    for (unsigned int j = pivot + 1; j < options.num_genes + 1; ++j ){
           
      front[j] = back[j] +  front[j-1] - back[j-1];
      
      // Now check for the largest R score
      float second_term =  static_cast<float>(i * j) / (options.num_genes * options.num_genes);
      float nvalue = (front[j] / total_weight) - second_term; 

      if (rvalue < nvalue) {
        rvalue = nvalue;
      } 
    }

    // Flip over the memory in the rmatrix (front row becomes rear row)
    front.swap(back);

  }*/

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

  fin.close();
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
 * Master MPI Process
 */

void masterProcess(Options &options){
  int total_tests =  options.num_lists * ( options.num_lists - 1) / 2;
  int processes_per_node = total_tests / (options.num_procs - 1);
  int extra_processes = total_tests % (options.num_procs - 1);

  vector < int > test_numbers;

  // Create our list of pairs
  for (int i=1; i < options.num_lists + 1; ++i){
    for (int j = i + 1; j < options.num_lists + 1; ++j){
      test_numbers.push_back( i );
      test_numbers.push_back( j );
    }
  }


  // Send the initial values to each client process, its total number of tests
  // and the test numbers themselves

  for (int i = 1; i < options.num_procs; ++i){

    int send_count = processes_per_node;
    if (i == options.num_procs - 1)
      send_count += extra_processes;

    int offset = processes_per_node * 2 * (i-1);
 
    send_count *= 2;
    MPI_Send(&send_count, 1, MPI_INT, i, 999, MPI_COMM_WORLD);
    MPI_Send(&test_numbers[offset], send_count, MPI_INT, i, 999, MPI_COMM_WORLD);

  } 

  // Now wait to receive the results

  vector <MPIResult> results;

  while(results.size() < total_tests){
    MPIResult mp;
    MPI_Status status;
    MPI_Recv(&mp, 1, resultType, MPI_ANY_SOURCE,  MPI_ANY_TAG,MPI_COMM_WORLD, &status );
    results.push_back(mp);
  }

  // Finish with the results
  for (MPIResult mp : results){
    std::cout << mp.i << "_" << mp.j << " " << mp.result << std::endl;
  }


}

/**
 * Client MPI Process
 */

void clientProcess(Options &options){

  MPI_Status status;
  int bsize;

  MPI_Recv(&bsize, 1, MPI_INT, 0, MPI_ANY_TAG,MPI_COMM_WORLD, &status);
  
  int buffer[bsize];
  MPI_Recv(&buffer, bsize, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

  // Now the proper loop begins

  for (int i = 0; i < bsize; i+=2){
    int l0 = buffer[i];
    int l1 = buffer[i+1];

  
    std::vector<unsigned int> gene_list0(options.num_genes);
    std::vector<unsigned int> gene_list1(options.num_genes);

    readLineIndex(options, l0, gene_list0);
    readLineIndex(options, l1, gene_list1);

    float rvalue = scoreLists(options, gene_list0, gene_list1 );

    MPIResult mp;
    mp.i = l0;
    mp.j = l1;
    mp.result = rvalue;

    MPI_Send(&mp, 1, resultType, 0, 999, MPI_COMM_WORLD);

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

  // MPI Init

  MPI_Init(&argc, const_cast<char***>(&argv));
  int length_name;
  char name[200];

  MPI_Comm_size(MPI_COMM_WORLD, &ops.num_procs);
  MPI_Comm_rank(MPI_COMM_WORLD, &ops.mpi_id);
  MPI_Get_processor_name(name, &length_name);
  
  if (ops.mpi_id == 0){
    cout << "MPI NumProcs: " << ops.num_procs << endl;
  }
  
  cout << "MPI ID: " << ops.mpi_id << " Name: " << name << endl;

  /// MPI Datatype for results
  MPI_Aint offsets[2], extent; 
  MPI_Datatype oldtypes[2];
  int blockcounts[2]; 

  offsets[0] = 0; 
  oldtypes[0] = MPI_INT; 
  blockcounts[0] = 2; 
 
  MPI_Type_extent(MPI_INT, &extent); 
  offsets[1] = 2 * extent; 
  oldtypes[1] = MPI_DOUBLE; 
  blockcounts[1] = 1; 

  MPI_Type_struct(2, blockcounts, offsets, oldtypes, &resultType); 
  MPI_Type_commit(&resultType); 


  // Read only, so mpi processes shouldnt clash with that
  readHeaderBlock(ops);

  if (ops.num_procs > 1){

    // The master MPI Process will farm out to the required processes
    if (ops.mpi_id == 0){
      masterProcess(ops);

    } else {
      clientProcess(ops);
    }

  } else {

    for (int i = 0; i < ops.num_lists; ++i ){
      for (int j = i + 1; j < ops.num_lists; ++j ){

        int l0 = i + 1;
        int l1 = j + 1;
        
        std::vector<unsigned int> gene_list0(ops.num_genes);
        std::vector<unsigned int> gene_list1(ops.num_genes);

        readLineIndex(ops, l0, gene_list0);
        readLineIndex(ops, l1, gene_list1);

        float rvalue = scoreLists(ops, gene_list0, gene_list1 );

        cout << l0 << "_" << l1 << " " << rvalue << endl;
      }
    }
  }

  MPI_Finalize();
}