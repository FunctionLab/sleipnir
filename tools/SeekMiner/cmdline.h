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
#define CMDLINE_PARSER_PACKAGE "SeekMiner"
#endif

#ifndef CMDLINE_PARSER_PACKAGE_NAME
/** @brief the complete program name (used for help and version) */
#define CMDLINE_PARSER_PACKAGE_NAME "SeekMiner"
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
  char * dset_arg;	/**< @brief Input a set of datasets.  */
  char * dset_orig;	/**< @brief Input a set of datasets original value given at command line.  */
  const char *dset_help; /**< @brief Input a set of datasets help description.  */
  char * search_dset_arg;	/**< @brief A set of datasets to search.  */
  char * search_dset_orig;	/**< @brief A set of datasets to search original value given at command line.  */
  const char *search_dset_help; /**< @brief A set of datasets to search help description.  */
  char * input_arg;	/**< @brief Input gene mapping.  */
  char * input_orig;	/**< @brief Input gene mapping original value given at command line.  */
  const char *input_help; /**< @brief Input gene mapping help description.  */
  char * query_arg;	/**< @brief Query gene list.  */
  char * query_orig;	/**< @brief Query gene list original value given at command line.  */
  const char *query_help; /**< @brief Query gene list help description.  */
  char * dir_in_arg;	/**< @brief Database directory.  */
  char * dir_in_orig;	/**< @brief Database directory original value given at command line.  */
  const char *dir_in_help; /**< @brief Database directory help description.  */
  char * dir_prep_in_arg;	/**< @brief Prep directory (containing .gavg, .gpres files).  */
  char * dir_prep_in_orig;	/**< @brief Prep directory (containing .gavg, .gpres files) original value given at command line.  */
  const char *dir_prep_in_help; /**< @brief Prep directory (containing .gavg, .gpres files) help description.  */
  char * dir_platform_arg;	/**< @brief Platform directory (containing .gplatavg, .gplatstdev, .gplatorder files).  */
  char * dir_platform_orig;	/**< @brief Platform directory (containing .gplatavg, .gplatstdev, .gplatorder files) original value given at command line.  */
  const char *dir_platform_help; /**< @brief Platform directory (containing .gplatavg, .gplatstdev, .gplatorder files) help description.  */
  char * quant_arg;	/**< @brief quant file (assuming all datasets use the same quantization).  */
  char * quant_orig;	/**< @brief quant file (assuming all datasets use the same quantization) original value given at command line.  */
  const char *quant_help; /**< @brief quant file (assuming all datasets use the same quantization) help description.  */
  int num_db_arg;	/**< @brief Number of databaselets in database (default='1000').  */
  char * num_db_orig;	/**< @brief Number of databaselets in database original value given at command line.  */
  const char *num_db_help; /**< @brief Number of databaselets in database help description.  */
  char * func_db_arg;	/**< @brief Functional network db path.  */
  char * func_db_orig;	/**< @brief Functional network db path original value given at command line.  */
  const char *func_db_help; /**< @brief Functional network db path help description.  */
  int func_n_arg;	/**< @brief Functional network number of databaselets (default='1000').  */
  char * func_n_orig;	/**< @brief Functional network number of databaselets original value given at command line.  */
  const char *func_n_help; /**< @brief Functional network number of databaselets help description.  */
  char * func_prep_arg;	/**< @brief Functional network prep & platform directory.  */
  char * func_prep_orig;	/**< @brief Functional network prep & platform directory original value given at command line.  */
  const char *func_prep_help; /**< @brief Functional network prep & platform directory help description.  */
  char * func_quant_arg;	/**< @brief Functional network quant file.  */
  char * func_quant_orig;	/**< @brief Functional network quant file original value given at command line.  */
  const char *func_quant_help; /**< @brief Functional network quant file help description.  */
  char * func_dset_arg;	/**< @brief Functional network dset-list file (1 dataset).  */
  char * func_dset_orig;	/**< @brief Functional network dset-list file (1 dataset) original value given at command line.  */
  const char *func_dset_help; /**< @brief Functional network dset-list file (1 dataset) help description.  */
  int func_logit_flag;	/**< @brief Functional network, integrate using logit values (default=on).  */
  const char *func_logit_help; /**< @brief Functional network, integrate using logit values help description.  */
  int norm_subavg_flag;	/**< @brief Per dataset, normalize z-scores by subtracting average of result gene (default=on).  */
  const char *norm_subavg_help; /**< @brief Per dataset, normalize z-scores by subtracting average of result gene help description.  */
  int norm_platsubavg_flag;	/**< @brief Per platform, normalize z-scores by subtracting average of query gene across platform (default=on).  */
  const char *norm_platsubavg_help; /**< @brief Per platform, normalize z-scores by subtracting average of query gene across platform help description.  */
  int norm_platstdev_flag;	/**< @brief Per platform, normalize z-scores by dividing stdev of query gene across platform (default=on).  */
  const char *norm_platstdev_help; /**< @brief Per platform, normalize z-scores by dividing stdev of query gene across platform help description.  */
  int is_nibble_flag;	/**< @brief Whether the input DB is nibble type (default=off).  */
  const char *is_nibble_help; /**< @brief Whether the input DB is nibble type help description.  */
  int buffer_arg;	/**< @brief Number of Databaselets to store in memory (default='20').  */
  char * buffer_orig;	/**< @brief Number of Databaselets to store in memory original value given at command line.  */
  const char *buffer_help; /**< @brief Number of Databaselets to store in memory help description.  */
  
  unsigned int help_given ;	/**< @brief Whether help was given.  */
  unsigned int version_given ;	/**< @brief Whether version was given.  */
  unsigned int dset_given ;	/**< @brief Whether dset was given.  */
  unsigned int search_dset_given ;	/**< @brief Whether search_dset was given.  */
  unsigned int input_given ;	/**< @brief Whether input was given.  */
  unsigned int query_given ;	/**< @brief Whether query was given.  */
  unsigned int dir_in_given ;	/**< @brief Whether dir_in was given.  */
  unsigned int dir_prep_in_given ;	/**< @brief Whether dir_prep_in was given.  */
  unsigned int dir_platform_given ;	/**< @brief Whether dir_platform was given.  */
  unsigned int quant_given ;	/**< @brief Whether quant was given.  */
  unsigned int num_db_given ;	/**< @brief Whether num_db was given.  */
  unsigned int func_db_given ;	/**< @brief Whether func_db was given.  */
  unsigned int func_n_given ;	/**< @brief Whether func_n was given.  */
  unsigned int func_prep_given ;	/**< @brief Whether func_prep was given.  */
  unsigned int func_quant_given ;	/**< @brief Whether func_quant was given.  */
  unsigned int func_dset_given ;	/**< @brief Whether func_dset was given.  */
  unsigned int func_logit_given ;	/**< @brief Whether func_logit was given.  */
  unsigned int norm_subavg_given ;	/**< @brief Whether norm_subavg was given.  */
  unsigned int norm_platsubavg_given ;	/**< @brief Whether norm_platsubavg was given.  */
  unsigned int norm_platstdev_given ;	/**< @brief Whether norm_platstdev was given.  */
  unsigned int is_nibble_given ;	/**< @brief Whether is_nibble was given.  */
  unsigned int buffer_given ;	/**< @brief Whether buffer was given.  */

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


#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /* CMDLINE_H */
