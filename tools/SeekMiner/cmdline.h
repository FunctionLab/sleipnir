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
#define CMDLINE_PARSER_PACKAGE "SeekMiner"
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
  char * dir_sinfo_arg;	/**< @brief Sinfo Directory (containing .sinfo files) (default='NA').  */
  char * dir_sinfo_orig;	/**< @brief Sinfo Directory (containing .sinfo files) original value given at command line.  */
  const char *dir_sinfo_help; /**< @brief Sinfo Directory (containing .sinfo files) help description.  */
  char * dir_gvar_arg;	/**< @brief Gene variance directory (containing .gexpvar files) (default='NA').  */
  char * dir_gvar_orig;	/**< @brief Gene variance directory (containing .gexpvar files) original value given at command line.  */
  const char *dir_gvar_help; /**< @brief Gene variance directory (containing .gexpvar files) help description.  */
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
  int func_logit_flag;	/**< @brief Functional network, integrate using logit values (default=off).  */
  const char *func_logit_help; /**< @brief Functional network, integrate using logit values help description.  */
  int correlation_flag;	/**< @brief Use Pearson correlation values, instead of z-score. -m, -M, -r do not apply (default=off).  */
  const char *correlation_help; /**< @brief Use Pearson correlation values, instead of z-score. -m, -M, -r do not apply help description.  */
  int norm_subavg_flag;	/**< @brief Per dataset, normalize z-scores by subtracting average of result gene (default=off).  */
  const char *norm_subavg_help; /**< @brief Per dataset, normalize z-scores by subtracting average of result gene help description.  */
  int norm_platsubavg_flag;	/**< @brief Per platform, normalize z-scores by subtracting average of query gene across platform (default=off).  */
  const char *norm_platsubavg_help; /**< @brief Per platform, normalize z-scores by subtracting average of query gene across platform help description.  */
  int norm_platstdev_flag;	/**< @brief Per platform, normalize z-scores by dividing stdev of query gene across platform (default=off).  */
  const char *norm_platstdev_help; /**< @brief Per platform, normalize z-scores by dividing stdev of query gene across platform help description.  */
  float score_cutoff_arg;	/**< @brief Cutoff on the gene-gene score before adding, default: no cutoff (default='-9999').  */
  char * score_cutoff_orig;	/**< @brief Cutoff on the gene-gene score before adding, default: no cutoff original value given at command line.  */
  const char *score_cutoff_help; /**< @brief Cutoff on the gene-gene score before adding, default: no cutoff help description.  */
  float per_q_required_arg;	/**< @brief Fraction (max 1.0) of query required to correlate with a gene, in order to count the gene's query score. A gene may not correlate with a query gene if it is absent, or its correlation with query does not pass cut-off (specified by --score_cutoff). Use this with caution. Be careful if using with --score_cutoff. (default='0.0').  */
  char * per_q_required_orig;	/**< @brief Fraction (max 1.0) of query required to correlate with a gene, in order to count the gene's query score. A gene may not correlate with a query gene if it is absent, or its correlation with query does not pass cut-off (specified by --score_cutoff). Use this with caution. Be careful if using with --score_cutoff. original value given at command line.  */
  const char *per_q_required_help; /**< @brief Fraction (max 1.0) of query required to correlate with a gene, in order to count the gene's query score. A gene may not correlate with a query gene if it is absent, or its correlation with query does not pass cut-off (specified by --score_cutoff). Use this with caution. Be careful if using with --score_cutoff. help description.  */
  int square_z_flag;	/**< @brief If using z-score, square-transform z-scores. Usually used in conjunction with --score-cutoff (default=off).  */
  const char *square_z_help; /**< @brief If using z-score, square-transform z-scores. Usually used in conjunction with --score-cutoff help description.  */
  float rank_biased_precision_p_arg;	/**< @brief Rank biased precision p parameter (default='0.99').  */
  char * rank_biased_precision_p_orig;	/**< @brief Rank biased precision p parameter original value given at command line.  */
  const char *rank_biased_precision_p_help; /**< @brief Rank biased precision p parameter help description.  */
  int random_flag;	/**< @brief Generate random ranking score (default=off).  */
  const char *random_help; /**< @brief Generate random ranking score help description.  */
  int num_random_arg;	/**< @brief Number of repetitions of generating random rankings (default='10').  */
  char * num_random_orig;	/**< @brief Number of repetitions of generating random rankings original value given at command line.  */
  const char *num_random_help; /**< @brief Number of repetitions of generating random rankings help description.  */
  int is_nibble_flag;	/**< @brief Whether the input DB is nibble type (default=off).  */
  const char *is_nibble_help; /**< @brief Whether the input DB is nibble type help description.  */
  int buffer_arg;	/**< @brief Number of Databaselets to store in memory (default='20').  */
  char * buffer_orig;	/**< @brief Number of Databaselets to store in memory original value given at command line.  */
  const char *buffer_help; /**< @brief Number of Databaselets to store in memory help description.  */
  int output_text_flag;	/**< @brief Output results (gene list and dataset weights) as text (default=off).  */
  const char *output_text_help; /**< @brief Output results (gene list and dataset weights) as text help description.  */
  char * output_dir_arg;	/**< @brief Output directory.  */
  char * output_dir_orig;	/**< @brief Output directory original value given at command line.  */
  const char *output_dir_help; /**< @brief Output directory help description.  */
  
  unsigned int help_given ;	/**< @brief Whether help was given.  */
  unsigned int version_given ;	/**< @brief Whether version was given.  */
  unsigned int dset_given ;	/**< @brief Whether dset was given.  */
  unsigned int search_dset_given ;	/**< @brief Whether search_dset was given.  */
  unsigned int input_given ;	/**< @brief Whether input was given.  */
  unsigned int query_given ;	/**< @brief Whether query was given.  */
  unsigned int dir_in_given ;	/**< @brief Whether dir_in was given.  */
  unsigned int dir_prep_in_given ;	/**< @brief Whether dir_prep_in was given.  */
  unsigned int dir_platform_given ;	/**< @brief Whether dir_platform was given.  */
  unsigned int dir_sinfo_given ;	/**< @brief Whether dir_sinfo was given.  */
  unsigned int dir_gvar_given ;	/**< @brief Whether dir_gvar was given.  */
  unsigned int quant_given ;	/**< @brief Whether quant was given.  */
  unsigned int num_db_given ;	/**< @brief Whether num_db was given.  */
  unsigned int func_db_given ;	/**< @brief Whether func_db was given.  */
  unsigned int func_n_given ;	/**< @brief Whether func_n was given.  */
  unsigned int func_prep_given ;	/**< @brief Whether func_prep was given.  */
  unsigned int func_quant_given ;	/**< @brief Whether func_quant was given.  */
  unsigned int func_dset_given ;	/**< @brief Whether func_dset was given.  */
  unsigned int func_logit_given ;	/**< @brief Whether func_logit was given.  */
  unsigned int correlation_given ;	/**< @brief Whether correlation was given.  */
  unsigned int norm_subavg_given ;	/**< @brief Whether norm_subavg was given.  */
  unsigned int norm_platsubavg_given ;	/**< @brief Whether norm_platsubavg was given.  */
  unsigned int norm_platstdev_given ;	/**< @brief Whether norm_platstdev was given.  */
  unsigned int score_cutoff_given ;	/**< @brief Whether score_cutoff was given.  */
  unsigned int per_q_required_given ;	/**< @brief Whether per_q_required was given.  */
  unsigned int square_z_given ;	/**< @brief Whether square_z was given.  */
  unsigned int rank_biased_precision_p_given ;	/**< @brief Whether rank_biased_precision_p was given.  */
  unsigned int random_given ;	/**< @brief Whether random was given.  */
  unsigned int num_random_given ;	/**< @brief Whether num_random was given.  */
  unsigned int is_nibble_given ;	/**< @brief Whether is_nibble was given.  */
  unsigned int buffer_given ;	/**< @brief Whether buffer was given.  */
  unsigned int output_text_given ;	/**< @brief Whether output_text was given.  */
  unsigned int output_dir_given ;	/**< @brief Whether output_dir was given.  */

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