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
#define CMDLINE_PARSER_PACKAGE "SeekAggregatedDataset"
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
  int pcl_flag;	/**< @brief PCL mode, suitable for dataset gene variance calculation (default=off).  */
  const char *pcl_help; /**< @brief PCL mode, suitable for dataset gene variance calculation help description.  */
  char * pcl_list_arg;	/**< @brief PCL list.  */
  char * pcl_list_orig;	/**< @brief PCL list original value given at command line.  */
  const char *pcl_list_help; /**< @brief PCL list help description.  */
  char * pcl_dir_arg;	/**< @brief PCL directory.  */
  char * pcl_dir_orig;	/**< @brief PCL directory original value given at command line.  */
  const char *pcl_dir_help; /**< @brief PCL directory help description.  */
  int step_num_arg;	/**< @brief Step Number (4 steps) (1: separate pairs to batches, 2: calculate Pearson for pairs in each batch (need a batch number), 3: merge Pearson from all batches and output a DAB) (default='0').  */
  char * step_num_orig;	/**< @brief Step Number (4 steps) (1: separate pairs to batches, 2: calculate Pearson for pairs in each batch (need a batch number), 3: merge Pearson from all batches and output a DAB) original value given at command line.  */
  const char *step_num_help; /**< @brief Step Number (4 steps) (1: separate pairs to batches, 2: calculate Pearson for pairs in each batch (need a batch number), 3: merge Pearson from all batches and output a DAB) help description.  */
  char * input_arg;	/**< @brief Gene mapping file.  */
  char * input_orig;	/**< @brief Gene mapping file original value given at command line.  */
  const char *input_help; /**< @brief Gene mapping file help description.  */
  char * query_arg;	/**< @brief Query file.  */
  char * query_orig;	/**< @brief Query file original value given at command line.  */
  const char *query_help; /**< @brief Query file help description.  */
  int num_batch_arg;	/**< @brief Number of batches to split pairs to (for step 1) (default='10').  */
  char * num_batch_orig;	/**< @brief Number of batches to split pairs to (for step 1) original value given at command line.  */
  const char *num_batch_help; /**< @brief Number of batches to split pairs to (for step 1) help description.  */
  char * pairs_dir_arg;	/**< @brief Pairs directory (for steps 1, 2). Pearson for the pairs will also be stored here..  */
  char * pairs_dir_orig;	/**< @brief Pairs directory (for steps 1, 2). Pearson for the pairs will also be stored here. original value given at command line.  */
  const char *pairs_dir_help; /**< @brief Pairs directory (for steps 1, 2). Pearson for the pairs will also be stored here. help description.  */
  int batch_num_arg;	/**< @brief Batch number (for step 2) (default='0').  */
  char * batch_num_orig;	/**< @brief Batch number (for step 2) original value given at command line.  */
  const char *batch_num_help; /**< @brief Batch number (for step 2) help description.  */
  char * dir_out_arg;	/**< @brief DAB output directory (for step 3).  */
  char * dir_out_orig;	/**< @brief DAB output directory (for step 3) original value given at command line.  */
  const char *dir_out_help; /**< @brief DAB output directory (for step 3) help description.  */
  
  unsigned int help_given ;	/**< @brief Whether help was given.  */
  unsigned int version_given ;	/**< @brief Whether version was given.  */
  unsigned int pcl_given ;	/**< @brief Whether pcl was given.  */
  unsigned int pcl_list_given ;	/**< @brief Whether pcl_list was given.  */
  unsigned int pcl_dir_given ;	/**< @brief Whether pcl_dir was given.  */
  unsigned int step_num_given ;	/**< @brief Whether step_num was given.  */
  unsigned int input_given ;	/**< @brief Whether input was given.  */
  unsigned int query_given ;	/**< @brief Whether query was given.  */
  unsigned int num_batch_given ;	/**< @brief Whether num_batch was given.  */
  unsigned int pairs_dir_given ;	/**< @brief Whether pairs_dir was given.  */
  unsigned int batch_num_given ;	/**< @brief Whether batch_num was given.  */
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


#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /* CMDLINE_H */
