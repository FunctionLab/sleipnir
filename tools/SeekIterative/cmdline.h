/** @file cmdline.h
 *  @brief The header file for the command line option parser
 *  generated by GNU Gengetopt version 2.22.5
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
/** @brief the program name (used for printing errors) */
#define CMDLINE_PARSER_PACKAGE "SeekIterative"
#endif

#ifndef CMDLINE_PARSER_PACKAGE_NAME
/** @brief the complete program name (used for help and version) */
#define CMDLINE_PARSER_PACKAGE_NAME "SeekIterative"
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
  int tdab_flag;	/**< @brief Traditional DAB mode (default=off).  */
  const char *tdab_help; /**< @brief Traditional DAB mode help description.  */
  int dab_flag;	/**< @brief Sparse Dab mode (default=off).  */
  const char *dab_help; /**< @brief Sparse Dab mode help description.  */
  int combined_flag;	/**< @brief Combined-dab mode (default=off).  */
  const char *combined_help; /**< @brief Combined-dab mode help description.  */
  int test_flag;	/**< @brief Test mode (default=off).  */
  const char *test_help; /**< @brief Test mode help description.  */
  int testcount_flag;	/**< @brief Test count mode (default=off).  */
  const char *testcount_help; /**< @brief Test count mode help description.  */
  int testcombined_flag;	/**< @brief Test count mode (default=off).  */
  const char *testcombined_help; /**< @brief Test count mode help description.  */
  int visualize_flag;	/**< @brief Visualization mode (default=off).  */
  const char *visualize_help; /**< @brief Visualization mode help description.  */
  char * dab_basename_arg;	/**< @brief Combined-dab basename, also shared with Test Mode.  */
  char * dab_basename_orig;	/**< @brief Combined-dab basename, also shared with Test Mode original value given at command line.  */
  const char *dab_basename_help; /**< @brief Combined-dab basename, also shared with Test Mode help description.  */
  int top_genes_arg;	/**< @brief Top genes to visualize (for combined-dab) (default='100').  */
  char * top_genes_orig;	/**< @brief Top genes to visualize (for combined-dab) original value given at command line.  */
  const char *top_genes_help; /**< @brief Top genes to visualize (for combined-dab) help description.  */
  int generate_dot_flag;	/**< @brief Generate a dot file (for combined-dab) (default=off).  */
  const char *generate_dot_help; /**< @brief Generate a dot file (for combined-dab) help description.  */
  int print_distr_flag;	/**< @brief Print distribution of edge values (default=off).  */
  const char *print_distr_help; /**< @brief Print distribution of edge values help description.  */
  float cutoff_arg;	/**< @brief Cutoff value (default='0.0001').  */
  char * cutoff_orig;	/**< @brief Cutoff value original value given at command line.  */
  const char *cutoff_help; /**< @brief Cutoff value help description.  */
  char * genome_arg;	/**< @brief Genome mapping file.  */
  char * genome_orig;	/**< @brief Genome mapping file original value given at command line.  */
  const char *genome_help; /**< @brief Genome mapping file help description.  */
  char * tdab_list_arg;	/**< @brief DAB list.  */
  char * tdab_list_orig;	/**< @brief DAB list original value given at command line.  */
  const char *tdab_list_help; /**< @brief DAB list help description.  */
  char * gavg_dir_arg;	/**< @brief Gene average directory (.gavg files) (default='NA').  */
  char * gavg_dir_orig;	/**< @brief Gene average directory (.gavg files) original value given at command line.  */
  const char *gavg_dir_help; /**< @brief Gene average directory (.gavg files) help description.  */
  char * dab_list_arg;	/**< @brief DAB list.  */
  char * dab_list_orig;	/**< @brief DAB list original value given at command line.  */
  const char *dab_list_help; /**< @brief DAB list help description.  */
  int num_iter_arg;	/**< @brief Number of iterations (default='0').  */
  char * num_iter_orig;	/**< @brief Number of iterations original value given at command line.  */
  const char *num_iter_help; /**< @brief Number of iterations help description.  */
  int default_type_arg;	/**< @brief Default gene index type (choose unsigned short for genes, or unsigned int (32-bit) for transcripts) (required for DAB mode) (0 - unsigned int, 1 - unsigned short) (default='-1').  */
  char * default_type_orig;	/**< @brief Default gene index type (choose unsigned short for genes, or unsigned int (32-bit) for transcripts) (required for DAB mode) (0 - unsigned int, 1 - unsigned short) original value given at command line.  */
  const char *default_type_help; /**< @brief Default gene index type (choose unsigned short for genes, or unsigned int (32-bit) for transcripts) (required for DAB mode) (0 - unsigned int, 1 - unsigned short) help description.  */
  float rbp_p_arg;	/**< @brief RBP p parameter (must be specified) (p<1.0) (recommended > 0.95) (default='-1').  */
  char * rbp_p_orig;	/**< @brief RBP p parameter (must be specified) (p<1.0) (recommended > 0.95) original value given at command line.  */
  const char *rbp_p_help; /**< @brief RBP p parameter (must be specified) (p<1.0) (recommended > 0.95) help description.  */
  int max_rank_arg;	/**< @brief Maximum rank number in the sparse DAB matrix (must be specified) (default='-1').  */
  char * max_rank_orig;	/**< @brief Maximum rank number in the sparse DAB matrix (must be specified) original value given at command line.  */
  const char *max_rank_help; /**< @brief Maximum rank number in the sparse DAB matrix (must be specified) help description.  */
  char * dset_cutoff_file_arg;	/**< @brief Dataset score cutoff file (default='NA').  */
  char * dset_cutoff_file_orig;	/**< @brief Dataset score cutoff file original value given at command line.  */
  const char *dset_cutoff_file_help; /**< @brief Dataset score cutoff file help description.  */
  char * norm_mode_arg;	/**< @brief Normalization method: rank - rank-normalize matrix, subtract_z - subtract-z-normalize matrix (default='NA').  */
  char * norm_mode_orig;	/**< @brief Normalization method: rank - rank-normalize matrix, subtract_z - subtract-z-normalize matrix original value given at command line.  */
  const char *norm_mode_help; /**< @brief Normalization method: rank - rank-normalize matrix, subtract_z - subtract-z-normalize matrix help description.  */
  float exp_arg;	/**< @brief Raise the z-score to the power of this value (for --norm_mode=subtract_z) (default='-1.0').  */
  char * exp_orig;	/**< @brief Raise the z-score to the power of this value (for --norm_mode=subtract_z) original value given at command line.  */
  const char *exp_help; /**< @brief Raise the z-score to the power of this value (for --norm_mode=subtract_z) help description.  */
  char * input_arg;	/**< @brief Gene mapping file.  */
  char * input_orig;	/**< @brief Gene mapping file original value given at command line.  */
  const char *input_help; /**< @brief Gene mapping file help description.  */
  char * query_arg;	/**< @brief Query file.  */
  char * query_orig;	/**< @brief Query file original value given at command line.  */
  const char *query_help; /**< @brief Query file help description.  */
  char * dab_dir_arg;	/**< @brief DAB directory.  */
  char * dab_dir_orig;	/**< @brief DAB directory original value given at command line.  */
  const char *dab_dir_help; /**< @brief DAB directory help description.  */
  char * not_query_arg;	/**< @brief NOT Query file (optional, for combined-DAB) (default='NA').  */
  char * not_query_orig;	/**< @brief NOT Query file (optional, for combined-DAB) original value given at command line.  */
  const char *not_query_help; /**< @brief NOT Query file (optional, for combined-DAB) help description.  */
  float threshold_q_arg;	/**< @brief Fraction of query genes need to be present in a dataset (default='0').  */
  char * threshold_q_orig;	/**< @brief Fraction of query genes need to be present in a dataset original value given at command line.  */
  const char *threshold_q_help; /**< @brief Fraction of query genes need to be present in a dataset help description.  */
  float threshold_g_arg;	/**< @brief Fraction of datasets that must contain a gene to put it in ranking (important if individual datasets have very different gene coverage, and for datasets with small gene-size) (default='0.50').  */
  char * threshold_g_orig;	/**< @brief Fraction of datasets that must contain a gene to put it in ranking (important if individual datasets have very different gene coverage, and for datasets with small gene-size) original value given at command line.  */
  const char *threshold_g_help; /**< @brief Fraction of datasets that must contain a gene to put it in ranking (important if individual datasets have very different gene coverage, and for datasets with small gene-size) help description.  */
  char * tsearch_mode_arg;	/**< @brief Search mode: equal weighted (eq) or CV LOI (cv_loi) or SPELL (spell) (Applicable if DAB list contains more than 1 dataset). (Required for --tdab and --dab modes) (default='NA').  */
  char * tsearch_mode_orig;	/**< @brief Search mode: equal weighted (eq) or CV LOI (cv_loi) or SPELL (spell) (Applicable if DAB list contains more than 1 dataset). (Required for --tdab and --dab modes) original value given at command line.  */
  const char *tsearch_mode_help; /**< @brief Search mode: equal weighted (eq) or CV LOI (cv_loi) or SPELL (spell) (Applicable if DAB list contains more than 1 dataset). (Required for --tdab and --dab modes) help description.  */
  int debug_flag;	/**< @brief Debug mode (default=off).  */
  const char *debug_help; /**< @brief Debug mode help description.  */
  char * dir_out_arg;	/**< @brief Output directory.  */
  char * dir_out_orig;	/**< @brief Output directory original value given at command line.  */
  const char *dir_out_help; /**< @brief Output directory help description.  */
  
  unsigned int help_given ;	/**< @brief Whether help was given.  */
  unsigned int version_given ;	/**< @brief Whether version was given.  */
  unsigned int tdab_given ;	/**< @brief Whether tdab was given.  */
  unsigned int dab_given ;	/**< @brief Whether dab was given.  */
  unsigned int combined_given ;	/**< @brief Whether combined was given.  */
  unsigned int test_given ;	/**< @brief Whether test was given.  */
  unsigned int testcount_given ;	/**< @brief Whether testcount was given.  */
  unsigned int testcombined_given ;	/**< @brief Whether testcombined was given.  */
  unsigned int visualize_given ;	/**< @brief Whether visualize was given.  */
  unsigned int dab_basename_given ;	/**< @brief Whether dab_basename was given.  */
  unsigned int top_genes_given ;	/**< @brief Whether top_genes was given.  */
  unsigned int generate_dot_given ;	/**< @brief Whether generate_dot was given.  */
  unsigned int print_distr_given ;	/**< @brief Whether print_distr was given.  */
  unsigned int cutoff_given ;	/**< @brief Whether cutoff was given.  */
  unsigned int genome_given ;	/**< @brief Whether genome was given.  */
  unsigned int tdab_list_given ;	/**< @brief Whether tdab_list was given.  */
  unsigned int gavg_dir_given ;	/**< @brief Whether gavg_dir was given.  */
  unsigned int dab_list_given ;	/**< @brief Whether dab_list was given.  */
  unsigned int num_iter_given ;	/**< @brief Whether num_iter was given.  */
  unsigned int default_type_given ;	/**< @brief Whether default_type was given.  */
  unsigned int rbp_p_given ;	/**< @brief Whether rbp_p was given.  */
  unsigned int max_rank_given ;	/**< @brief Whether max_rank was given.  */
  unsigned int dset_cutoff_file_given ;	/**< @brief Whether dset_cutoff_file was given.  */
  unsigned int norm_mode_given ;	/**< @brief Whether norm_mode was given.  */
  unsigned int exp_given ;	/**< @brief Whether exp was given.  */
  unsigned int input_given ;	/**< @brief Whether input was given.  */
  unsigned int query_given ;	/**< @brief Whether query was given.  */
  unsigned int dab_dir_given ;	/**< @brief Whether dab_dir was given.  */
  unsigned int not_query_given ;	/**< @brief Whether not_query was given.  */
  unsigned int threshold_q_given ;	/**< @brief Whether threshold_q was given.  */
  unsigned int threshold_g_given ;	/**< @brief Whether threshold_g was given.  */
  unsigned int tsearch_mode_given ;	/**< @brief Whether tsearch_mode was given.  */
  unsigned int debug_given ;	/**< @brief Whether debug was given.  */
  unsigned int dir_out_given ;	/**< @brief Whether dir_out was given.  */

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
int cmdline_parser (int argc, char **argv,
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
int cmdline_parser2 (int argc, char **argv,
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
int cmdline_parser_ext (int argc, char **argv,
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

extern const char *cmdline_parser_norm_mode_values[];  /**< @brief Possible values for norm_mode. */
extern const char *cmdline_parser_tsearch_mode_values[];  /**< @brief Possible values for tsearch_mode. */


#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /* CMDLINE_H */
