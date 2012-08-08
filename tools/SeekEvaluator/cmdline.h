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
#define CMDLINE_PARSER_PACKAGE "SeekEvaluator"
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
  int single_flag;	/**< @brief Evaluate one query's rank list result (default=off).  */
  const char *single_help; /**< @brief Evaluate one query's rank list result help description.  */
  int aggregate_flag;	/**< @brief Evaluate multiple queries and aggregates results (default=off).  */
  const char *aggregate_help; /**< @brief Evaluate multiple queries and aggregates results help description.  */
  int rbp_flag;	/**< @brief Rank biased precision (requires parameter p to be set) (default=off).  */
  const char *rbp_help; /**< @brief Rank biased precision (requires parameter p to be set) help description.  */
  int avgp_flag;	/**< @brief Average precision for the top X positives, where X = integer, or % of total positives (default=off).  */
  const char *avgp_help; /**< @brief Average precision for the top X positives, where X = integer, or % of total positives help description.  */
  int pr_flag;	/**< @brief Precision at X-th positive, where X = integer, or % of total positives (default=off).  */
  const char *pr_help; /**< @brief Precision at X-th positive, where X = integer, or % of total positives help description.  */
  int pr_all_flag;	/**< @brief Precision at all positive depths (useful for drawing precision-recall curve) (default=off).  */
  const char *pr_all_help; /**< @brief Precision at all positive depths (useful for drawing precision-recall curve) help description.  */
  int auc_flag;	/**< @brief AUC (default=off).  */
  const char *auc_help; /**< @brief AUC help description.  */
  int x_int_arg;	/**< @brief Parameter X = integer, for --avgp, --pr (default='-1').  */
  char * x_int_orig;	/**< @brief Parameter X = integer, for --avgp, --pr original value given at command line.  */
  const char *x_int_help; /**< @brief Parameter X = integer, for --avgp, --pr help description.  */
  float x_per_arg;	/**< @brief Parameter X = percentage, for --avgp, --pr (default='0').  */
  char * x_per_orig;	/**< @brief Parameter X = percentage, for --avgp, --pr original value given at command line.  */
  const char *x_per_help; /**< @brief Parameter X = percentage, for --avgp, --pr help description.  */
  float rbp_p_arg;	/**< @brief Parameter p, for --rbp (default='0.95').  */
  char * rbp_p_orig;	/**< @brief Parameter p, for --rbp original value given at command line.  */
  const char *rbp_p_help; /**< @brief Parameter p, for --rbp help description.  */
  int dislay_only_flag;	/**< @brief Display the genes sorted by score (top 500 is shown, for single mode only) (default=off).  */
  const char *dislay_only_help; /**< @brief Display the genes sorted by score (top 500 is shown, for single mode only) help description.  */
  int display_weight_flag;	/**< @brief Display dataset weights (top 100) (default=off).  */
  const char *display_weight_help; /**< @brief Display dataset weights (top 100) help description.  */
  int agg_avg_flag;	/**< @brief Show the average, standard deviation of the metric for all queries (default=off).  */
  const char *agg_avg_help; /**< @brief Show the average, standard deviation of the metric for all queries help description.  */
  int agg_quartile_flag;	/**< @brief Show the min, max, as well as the 1st, 2nd, 3rd quartile of the metric for all queries (default=off).  */
  const char *agg_quartile_help; /**< @brief Show the min, max, as well as the 1st, 2nd, 3rd quartile of the metric for all queries help description.  */
  int agg_ranksum_flag;	/**< @brief Sum up the ranks of genes in all query rankings to produce a master list sorted by summed rank, and perform metric on this list (default=off).  */
  const char *agg_ranksum_help; /**< @brief Sum up the ranks of genes in all query rankings to produce a master list sorted by summed rank, and perform metric on this list help description.  */
  int agg_scoresum_flag;	/**< @brief Sum up the scores of genes in all query rankings to produce a master list sorted by summed score, and perform metric on this list (default=off).  */
  const char *agg_scoresum_help; /**< @brief Sum up the scores of genes in all query rankings to produce a master list sorted by summed score, and perform metric on this list help description.  */
  int display_all_flag;	/**< @brief Display the metric for all queries (default=off).  */
  const char *display_all_help; /**< @brief Display the metric for all queries help description.  */
  char * input_arg;	/**< @brief Gene mapping file.  */
  char * input_orig;	/**< @brief Gene mapping file original value given at command line.  */
  const char *input_help; /**< @brief Gene mapping file help description.  */
  char * dataset_map_arg;	/**< @brief Dataset mapping file, only required for displaying dataset weights.  */
  char * dataset_map_orig;	/**< @brief Dataset mapping file, only required for displaying dataset weights original value given at command line.  */
  const char *dataset_map_help; /**< @brief Dataset mapping file, only required for displaying dataset weights help description.  */
  char * weight_arg;	/**< @brief Dataset weight file, (*.dweight).  */
  char * weight_orig;	/**< @brief Dataset weight file, (*.dweight) original value given at command line.  */
  const char *weight_help; /**< @brief Dataset weight file, (*.dweight) help description.  */
  char * goldstd_arg;	/**< @brief Gold standard gene set file (one line, space delimited).  */
  char * goldstd_orig;	/**< @brief Gold standard gene set file (one line, space delimited) original value given at command line.  */
  const char *goldstd_help; /**< @brief Gold standard gene set file (one line, space delimited) help description.  */
  char * gscore_arg;	/**< @brief Gene score file (.gscore).  */
  char * gscore_orig;	/**< @brief Gene score file (.gscore) original value given at command line.  */
  const char *gscore_help; /**< @brief Gene score file (.gscore) help description.  */
  char * query_arg;	/**< @brief Query gene set file (to be excluded from evaluation) (.query).  */
  char * query_orig;	/**< @brief Query gene set file (to be excluded from evaluation) (.query) original value given at command line.  */
  const char *query_help; /**< @brief Query gene set file (to be excluded from evaluation) (.query) help description.  */
  float nan_arg;	/**< @brief Define NaN score (default='-320').  */
  char * nan_orig;	/**< @brief Define NaN score original value given at command line.  */
  const char *nan_help; /**< @brief Define NaN score help description.  */
  char * goldstd_list_arg;	/**< @brief List of gold standard gene set files.  */
  char * goldstd_list_orig;	/**< @brief List of gold standard gene set files original value given at command line.  */
  const char *goldstd_list_help; /**< @brief List of gold standard gene set files help description.  */
  char * gscore_list_arg;	/**< @brief List of gene score files.  */
  char * gscore_list_orig;	/**< @brief List of gene score files original value given at command line.  */
  const char *gscore_list_help; /**< @brief List of gene score files help description.  */
  char * query_list_arg;	/**< @brief List of query gene set files.  */
  char * query_list_orig;	/**< @brief List of query gene set files original value given at command line.  */
  const char *query_list_help; /**< @brief List of query gene set files help description.  */
  char * dir_out_arg;	/**< @brief Output directory.  */
  char * dir_out_orig;	/**< @brief Output directory original value given at command line.  */
  const char *dir_out_help; /**< @brief Output directory help description.  */
  
  unsigned int help_given ;	/**< @brief Whether help was given.  */
  unsigned int version_given ;	/**< @brief Whether version was given.  */
  unsigned int single_given ;	/**< @brief Whether single was given.  */
  unsigned int aggregate_given ;	/**< @brief Whether aggregate was given.  */
  unsigned int rbp_given ;	/**< @brief Whether rbp was given.  */
  unsigned int avgp_given ;	/**< @brief Whether avgp was given.  */
  unsigned int pr_given ;	/**< @brief Whether pr was given.  */
  unsigned int pr_all_given ;	/**< @brief Whether pr_all was given.  */
  unsigned int auc_given ;	/**< @brief Whether auc was given.  */
  unsigned int x_int_given ;	/**< @brief Whether x_int was given.  */
  unsigned int x_per_given ;	/**< @brief Whether x_per was given.  */
  unsigned int rbp_p_given ;	/**< @brief Whether rbp_p was given.  */
  unsigned int dislay_only_given ;	/**< @brief Whether dislay_only was given.  */
  unsigned int display_weight_given ;	/**< @brief Whether display_weight was given.  */
  unsigned int agg_avg_given ;	/**< @brief Whether agg_avg was given.  */
  unsigned int agg_quartile_given ;	/**< @brief Whether agg_quartile was given.  */
  unsigned int agg_ranksum_given ;	/**< @brief Whether agg_ranksum was given.  */
  unsigned int agg_scoresum_given ;	/**< @brief Whether agg_scoresum was given.  */
  unsigned int display_all_given ;	/**< @brief Whether display_all was given.  */
  unsigned int input_given ;	/**< @brief Whether input was given.  */
  unsigned int dataset_map_given ;	/**< @brief Whether dataset_map was given.  */
  unsigned int weight_given ;	/**< @brief Whether weight was given.  */
  unsigned int goldstd_given ;	/**< @brief Whether goldstd was given.  */
  unsigned int gscore_given ;	/**< @brief Whether gscore was given.  */
  unsigned int query_given ;	/**< @brief Whether query was given.  */
  unsigned int nan_given ;	/**< @brief Whether nan was given.  */
  unsigned int goldstd_list_given ;	/**< @brief Whether goldstd_list was given.  */
  unsigned int gscore_list_given ;	/**< @brief Whether gscore_list was given.  */
  unsigned int query_list_given ;	/**< @brief Whether query_list was given.  */
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
