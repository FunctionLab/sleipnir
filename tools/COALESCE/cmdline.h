/** @file cmdline.h
 *  @brief The header file for the command line option parser
 *  generated by GNU Gengetopt version 2.22
 *  http://www.gnu.org/software/gengetopt.
 *  DO NOT modify this file, since it can be overwritten
 *  @author GNU Gengetopt by Lorenzo Bettini */

#ifndef CMDLINE_H
#define CMDLINE_H

/* If we use autoconf.  */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h> /* for FILE */

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#ifndef CMDLINE_PARSER_PACKAGE
/** @brief the program name */
#define CMDLINE_PARSER_PACKAGE "COALESCE"
#endif

#ifndef CMDLINE_PARSER_VERSION
/** @brief the program version */
#define CMDLINE_PARSER_VERSION "1.0"
#endif

/** @brief Where the command line options are stored */
struct gengetopt_args_info
{
  const char *help_help; /**< @brief Print help and exit help description.  */
  const char *version_help; /**< @brief Print version and exit help description.  */
  char * input_arg;	/**< @brief Input PCL file.  */
  char * input_orig;	/**< @brief Input PCL file original value given at command line.  */
  const char *input_help; /**< @brief Input PCL file help description.  */
  char * fasta_arg;	/**< @brief Input FASTA file.  */
  char * fasta_orig;	/**< @brief Input FASTA file original value given at command line.  */
  const char *fasta_help; /**< @brief Input FASTA file help description.  */
  char * datasets_arg;	/**< @brief Condition groupings into dataset blocks.  */
  char * datasets_orig;	/**< @brief Condition groupings into dataset blocks original value given at command line.  */
  const char *datasets_help; /**< @brief Condition groupings into dataset blocks help description.  */
  char * output_arg;	/**< @brief Directory for output files (PCLs/motifs).  */
  char * output_orig;	/**< @brief Directory for output files (PCLs/motifs) original value given at command line.  */
  const char *output_help; /**< @brief Directory for output files (PCLs/motifs) help description.  */
  double prob_gene_arg;	/**< @brief Probability threshhold for gene inclusion (default='0.95').  */
  char * prob_gene_orig;	/**< @brief Probability threshhold for gene inclusion original value given at command line.  */
  const char *prob_gene_help; /**< @brief Probability threshhold for gene inclusion help description.  */
  double pvalue_cond_arg;	/**< @brief P-value threshhold for condition inclusion (default='0.05').  */
  char * pvalue_cond_orig;	/**< @brief P-value threshhold for condition inclusion original value given at command line.  */
  const char *pvalue_cond_help; /**< @brief P-value threshhold for condition inclusion help description.  */
  double pvalue_motif_arg;	/**< @brief P-value threshhold for motif inclusion (default='0.05').  */
  char * pvalue_motif_orig;	/**< @brief P-value threshhold for motif inclusion original value given at command line.  */
  const char *pvalue_motif_help; /**< @brief P-value threshhold for motif inclusion help description.  */
  double zscore_cond_arg;	/**< @brief Z-score threshhold for condition inclusion (default='0.5').  */
  char * zscore_cond_orig;	/**< @brief Z-score threshhold for condition inclusion original value given at command line.  */
  const char *zscore_cond_help; /**< @brief Z-score threshhold for condition inclusion help description.  */
  double zscore_motif_arg;	/**< @brief Z-score threshhold for motif inclusion (default='0.5').  */
  char * zscore_motif_orig;	/**< @brief Z-score threshhold for motif inclusion original value given at command line.  */
  const char *zscore_motif_help; /**< @brief Z-score threshhold for motif inclusion help description.  */
  int k_arg;	/**< @brief Sequence kmer length (default='7').  */
  char * k_orig;	/**< @brief Sequence kmer length original value given at command line.  */
  const char *k_help; /**< @brief Sequence kmer length help description.  */
  double pvalue_merge_arg;	/**< @brief P-value threshhold for motif merging (default='0.05').  */
  char * pvalue_merge_orig;	/**< @brief P-value threshhold for motif merging original value given at command line.  */
  const char *pvalue_merge_help; /**< @brief P-value threshhold for motif merging help description.  */
  double cutoff_merge_arg;	/**< @brief Edit distance cutoff for motif merging (default='2.5').  */
  char * cutoff_merge_orig;	/**< @brief Edit distance cutoff for motif merging original value given at command line.  */
  const char *cutoff_merge_help; /**< @brief Edit distance cutoff for motif merging help description.  */
  double penalty_gap_arg;	/**< @brief Edit distance penalty for gaps (default='1').  */
  char * penalty_gap_orig;	/**< @brief Edit distance penalty for gaps original value given at command line.  */
  const char *penalty_gap_help; /**< @brief Edit distance penalty for gaps help description.  */
  double penalty_mismatch_arg;	/**< @brief Edit distance penalty for mismatches (default='2.1').  */
  char * penalty_mismatch_orig;	/**< @brief Edit distance penalty for mismatches original value given at command line.  */
  const char *penalty_mismatch_help; /**< @brief Edit distance penalty for mismatches help description.  */
  double pvalue_correl_arg;	/**< @brief P-value threshhold for significant correlation (default='0.05').  */
  char * pvalue_correl_orig;	/**< @brief P-value threshhold for significant correlation original value given at command line.  */
  const char *pvalue_correl_help; /**< @brief P-value threshhold for significant correlation help description.  */
  int number_correl_arg;	/**< @brief Maximum number of pairs to sample for significant correlation (default='100000').  */
  char * number_correl_orig;	/**< @brief Maximum number of pairs to sample for significant correlation original value given at command line.  */
  const char *number_correl_help; /**< @brief Maximum number of pairs to sample for significant correlation help description.  */
  char * sequences_arg;	/**< @brief Sequence types to use (comma separated).  */
  char * sequences_orig;	/**< @brief Sequence types to use (comma separated) original value given at command line.  */
  const char *sequences_help; /**< @brief Sequence types to use (comma separated) help description.  */
  int bases_arg;	/**< @brief Resolution of bases per motif match (default='5000').  */
  char * bases_orig;	/**< @brief Resolution of bases per motif match original value given at command line.  */
  const char *bases_help; /**< @brief Resolution of bases per motif match help description.  */
  int size_minimum_arg;	/**< @brief Minimum gene count for clusters of interest (default='5').  */
  char * size_minimum_orig;	/**< @brief Minimum gene count for clusters of interest original value given at command line.  */
  const char *size_minimum_help; /**< @brief Minimum gene count for clusters of interest help description.  */
  int size_merge_arg;	/**< @brief Maximum motif count for realtime merging (default='100').  */
  char * size_merge_orig;	/**< @brief Maximum motif count for realtime merging original value given at command line.  */
  const char *size_merge_help; /**< @brief Maximum motif count for realtime merging help description.  */
  int size_maximum_arg;	/**< @brief Maximum motif count to consider a cluster saturated (default='1000').  */
  char * size_maximum_orig;	/**< @brief Maximum motif count to consider a cluster saturated original value given at command line.  */
  const char *size_maximum_help; /**< @brief Maximum motif count to consider a cluster saturated help description.  */
  char * postprocess_arg;	/**< @brief Input file/directory of clusters to postprocess.  */
  char * postprocess_orig;	/**< @brief Input file/directory of clusters to postprocess original value given at command line.  */
  const char *postprocess_help; /**< @brief Input file/directory of clusters to postprocess help description.  */
  char * known_motifs_arg;	/**< @brief File containing known motifs.  */
  char * known_motifs_orig;	/**< @brief File containing known motifs original value given at command line.  */
  const char *known_motifs_help; /**< @brief File containing known motifs help description.  */
  double known_cutoff_arg;	/**< @brief Score cutoff for known motif labeling (default='0.05').  */
  char * known_cutoff_orig;	/**< @brief Score cutoff for known motif labeling original value given at command line.  */
  const char *known_cutoff_help; /**< @brief Score cutoff for known motif labeling help description.  */
  char * known_type_arg;	/**< @brief Type of known motif matching (default='pvalue').  */
  char * known_type_orig;	/**< @brief Type of known motif matching original value given at command line.  */
  const char *known_type_help; /**< @brief Type of known motif matching help description.  */
  double cutoff_postprocess_arg;	/**< @brief Similarity cutoff for cluster merging (default='1').  */
  char * cutoff_postprocess_orig;	/**< @brief Similarity cutoff for cluster merging original value given at command line.  */
  const char *cutoff_postprocess_help; /**< @brief Similarity cutoff for cluster merging help description.  */
  double fraction_postprocess_arg;	/**< @brief Overlap fraction for postprocessing gene/condition inclusion (default='0.5').  */
  char * fraction_postprocess_orig;	/**< @brief Overlap fraction for postprocessing gene/condition inclusion original value given at command line.  */
  const char *fraction_postprocess_help; /**< @brief Overlap fraction for postprocessing gene/condition inclusion help description.  */
  double cutoff_trim_arg;	/**< @brief Cocluster stdev cutoff for cluster trimming (default='1').  */
  char * cutoff_trim_orig;	/**< @brief Cocluster stdev cutoff for cluster trimming original value given at command line.  */
  const char *cutoff_trim_help; /**< @brief Cocluster stdev cutoff for cluster trimming help description.  */
  int remove_rcs_flag;	/**< @brief Convert RCs and RC-like PSTs to single strand (default=on).  */
  const char *remove_rcs_help; /**< @brief Convert RCs and RC-like PSTs to single strand help description.  */
  double min_info_arg;	/**< @brief Uninformative motif threshhold (bits) (default='0.3').  */
  char * min_info_orig;	/**< @brief Uninformative motif threshhold (bits) original value given at command line.  */
  const char *min_info_help; /**< @brief Uninformative motif threshhold (bits) help description.  */
  int max_motifs_arg;	/**< @brief Maximum motifs to merge exactly (default='2500').  */
  char * max_motifs_orig;	/**< @brief Maximum motifs to merge exactly original value given at command line.  */
  const char *max_motifs_help; /**< @brief Maximum motifs to merge exactly help description.  */
  int normalize_flag;	/**< @brief Automatically detect/normalize single channel data (default=off).  */
  const char *normalize_help; /**< @brief Automatically detect/normalize single channel data help description.  */
  int progressive_flag;	/**< @brief Generate output progressively (default=on).  */
  const char *progressive_help; /**< @brief Generate output progressively help description.  */
  int threads_arg;	/**< @brief Maximum number of concurrent threads (default='1').  */
  char * threads_orig;	/**< @brief Maximum number of concurrent threads original value given at command line.  */
  const char *threads_help; /**< @brief Maximum number of concurrent threads help description.  */
  int skip_arg;	/**< @brief Columns to skip in input PCL (default='2').  */
  char * skip_orig;	/**< @brief Columns to skip in input PCL original value given at command line.  */
  const char *skip_help; /**< @brief Columns to skip in input PCL help description.  */
  int random_arg;	/**< @brief Seed random generator (default='0').  */
  char * random_orig;	/**< @brief Seed random generator original value given at command line.  */
  const char *random_help; /**< @brief Seed random generator help description.  */
  int verbosity_arg;	/**< @brief Message verbosity (default='5').  */
  char * verbosity_orig;	/**< @brief Message verbosity original value given at command line.  */
  const char *verbosity_help; /**< @brief Message verbosity help description.  */
  
  unsigned int help_given ;	/**< @brief Whether help was given.  */
  unsigned int version_given ;	/**< @brief Whether version was given.  */
  unsigned int input_given ;	/**< @brief Whether input was given.  */
  unsigned int fasta_given ;	/**< @brief Whether fasta was given.  */
  unsigned int datasets_given ;	/**< @brief Whether datasets was given.  */
  unsigned int output_given ;	/**< @brief Whether output was given.  */
  unsigned int prob_gene_given ;	/**< @brief Whether prob_gene was given.  */
  unsigned int pvalue_cond_given ;	/**< @brief Whether pvalue_cond was given.  */
  unsigned int pvalue_motif_given ;	/**< @brief Whether pvalue_motif was given.  */
  unsigned int zscore_cond_given ;	/**< @brief Whether zscore_cond was given.  */
  unsigned int zscore_motif_given ;	/**< @brief Whether zscore_motif was given.  */
  unsigned int k_given ;	/**< @brief Whether k was given.  */
  unsigned int pvalue_merge_given ;	/**< @brief Whether pvalue_merge was given.  */
  unsigned int cutoff_merge_given ;	/**< @brief Whether cutoff_merge was given.  */
  unsigned int penalty_gap_given ;	/**< @brief Whether penalty_gap was given.  */
  unsigned int penalty_mismatch_given ;	/**< @brief Whether penalty_mismatch was given.  */
  unsigned int pvalue_correl_given ;	/**< @brief Whether pvalue_correl was given.  */
  unsigned int number_correl_given ;	/**< @brief Whether number_correl was given.  */
  unsigned int sequences_given ;	/**< @brief Whether sequences was given.  */
  unsigned int bases_given ;	/**< @brief Whether bases was given.  */
  unsigned int size_minimum_given ;	/**< @brief Whether size_minimum was given.  */
  unsigned int size_merge_given ;	/**< @brief Whether size_merge was given.  */
  unsigned int size_maximum_given ;	/**< @brief Whether size_maximum was given.  */
  unsigned int postprocess_given ;	/**< @brief Whether postprocess was given.  */
  unsigned int known_motifs_given ;	/**< @brief Whether known_motifs was given.  */
  unsigned int known_cutoff_given ;	/**< @brief Whether known_cutoff was given.  */
  unsigned int known_type_given ;	/**< @brief Whether known_type was given.  */
  unsigned int cutoff_postprocess_given ;	/**< @brief Whether cutoff_postprocess was given.  */
  unsigned int fraction_postprocess_given ;	/**< @brief Whether fraction_postprocess was given.  */
  unsigned int cutoff_trim_given ;	/**< @brief Whether cutoff_trim was given.  */
  unsigned int remove_rcs_given ;	/**< @brief Whether remove_rcs was given.  */
  unsigned int min_info_given ;	/**< @brief Whether min_info was given.  */
  unsigned int max_motifs_given ;	/**< @brief Whether max_motifs was given.  */
  unsigned int normalize_given ;	/**< @brief Whether normalize was given.  */
  unsigned int progressive_given ;	/**< @brief Whether progressive was given.  */
  unsigned int threads_given ;	/**< @brief Whether threads was given.  */
  unsigned int skip_given ;	/**< @brief Whether skip was given.  */
  unsigned int random_given ;	/**< @brief Whether random was given.  */
  unsigned int verbosity_given ;	/**< @brief Whether verbosity was given.  */

  char **inputs ; /**< @brief unamed options (options without names) */
  unsigned inputs_num ; /**< @brief unamed options number */
} ;

/** @brief The additional parameters to pass to parser functions */
struct cmdline_parser_params
{
  int override; /**< @brief whether to override possibly already present options (default 0) */
  int initialize; /**< @brief whether to initialize the option structure gengetopt_args_info (default 1) */
  int check_required; /**< @brief whether to check that all required options were provided (default 1) */
  int check_ambiguity; /**< @brief whether to check for options already specified in the option structure gengetopt_args_info (default 0) */
  int print_errors; /**< @brief whether getopt_long should print an error message for a bad option (default 1) */
} ;

/** @brief the purpose string of the program */
extern const char *gengetopt_args_info_purpose;
/** @brief the usage string of the program */
extern const char *gengetopt_args_info_usage;
/** @brief all the lines making the help output */
extern const char *gengetopt_args_info_help[];

/**
 * The command line parser
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser (int argc, char * const *argv,
  struct gengetopt_args_info *args_info);

/**
 * The command line parser (version with additional parameters - deprecated)
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @param override whether to override possibly already present options
 * @param initialize whether to initialize the option structure my_args_info
 * @param check_required whether to check that all required options were provided
 * @return 0 if everything went fine, NON 0 if an error took place
 * @deprecated use cmdline_parser_ext() instead
 */
int cmdline_parser2 (int argc, char * const *argv,
  struct gengetopt_args_info *args_info,
  int override, int initialize, int check_required);

/**
 * The command line parser (version with additional parameters)
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @param params additional parameters for the parser
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser_ext (int argc, char * const *argv,
  struct gengetopt_args_info *args_info,
  struct cmdline_parser_params *params);

/**
 * Save the contents of the option struct into an already open FILE stream.
 * @param outfile the stream where to dump options
 * @param args_info the option struct to dump
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser_dump(FILE *outfile,
  struct gengetopt_args_info *args_info);

/**
 * Save the contents of the option struct into a (text) file.
 * This file can be read by the config file parser (if generated by gengetopt)
 * @param filename the file where to save
 * @param args_info the option struct to save
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser_file_save(const char *filename,
  struct gengetopt_args_info *args_info);

/**
 * Print the help
 */
void cmdline_parser_print_help(void);
/**
 * Print the version
 */
void cmdline_parser_print_version(void);

/**
 * Initializes all the fields a cmdline_parser_params structure 
 * to their default values
 * @param params the structure to initialize
 */
void cmdline_parser_params_init(struct cmdline_parser_params *params);

/**
 * Allocates dynamically a cmdline_parser_params structure and initializes
 * all its fields to their default values
 * @return the created and initialized cmdline_parser_params structure
 */
struct cmdline_parser_params *cmdline_parser_params_create(void);

/**
 * Initializes the passed gengetopt_args_info structure's fields
 * (also set default values for options that have a default)
 * @param args_info the structure to initialize
 */
void cmdline_parser_init (struct gengetopt_args_info *args_info);
/**
 * Deallocates the string fields of the gengetopt_args_info structure
 * (but does not deallocate the structure itself)
 * @param args_info the structure to deallocate
 */
void cmdline_parser_free (struct gengetopt_args_info *args_info);

/**
 * Checks that all the required options were specified
 * @param args_info the structure to check
 * @param prog_name the name of the program that will be used to print
 *   possible errors
 * @return
 */
int cmdline_parser_required (struct gengetopt_args_info *args_info,
  const char *prog_name);

extern char *cmdline_parser_known_type_values[] ;	/**< @brief Possible values for known_type.  */


#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /* CMDLINE_H */
