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
  char * search_dset_arg;	/**< @brief A set of datasets to search. If not specified, search all datasets. (default='NA').  */
  char * search_dset_orig;	/**< @brief A set of datasets to search. If not specified, search all datasets. original value given at command line.  */
  const char *search_dset_help; /**< @brief A set of datasets to search. If not specified, search all datasets. help description.  */
  char * input_arg;	/**< @brief Input gene mapping.  */
  char * input_orig;	/**< @brief Input gene mapping original value given at command line.  */
  const char *input_help; /**< @brief Input gene mapping help description.  */
  char * query_arg;	/**< @brief A list of queries (each line is a query, and contains the names of the query genes delimited by spaces).  */
  char * query_orig;	/**< @brief A list of queries (each line is a query, and contains the names of the query genes delimited by spaces) original value given at command line.  */
  const char *query_help; /**< @brief A list of queries (each line is a query, and contains the names of the query genes delimited by spaces) help description.  */
  char * dir_in_arg;	/**< @brief Database directory (containing .db files).  */
  char * dir_in_orig;	/**< @brief Database directory (containing .db files) original value given at command line.  */
  const char *dir_in_help; /**< @brief Database directory (containing .db files) help description.  */
  char * dir_prep_in_arg;	/**< @brief Prep directory (containing .gavg, .gpres files).  */
  char * dir_prep_in_orig;	/**< @brief Prep directory (containing .gavg, .gpres files) original value given at command line.  */
  const char *dir_prep_in_help; /**< @brief Prep directory (containing .gavg, .gpres files) help description.  */
  char * dir_platform_arg;	/**< @brief Platform directory (containing .gplatavg, .gplatstdev, .gplatorder files).  */
  char * dir_platform_orig;	/**< @brief Platform directory (containing .gplatavg, .gplatstdev, .gplatorder files) original value given at command line.  */
  const char *dir_platform_help; /**< @brief Platform directory (containing .gplatavg, .gplatstdev, .gplatorder files) help description.  */
  char * dir_sinfo_arg;	/**< @brief Sinfo directory (containing .sinfo files) (default='NA').  */
  char * dir_sinfo_orig;	/**< @brief Sinfo directory (containing .sinfo files) original value given at command line.  */
  const char *dir_sinfo_help; /**< @brief Sinfo directory (containing .sinfo files) help description.  */
  char * dir_gvar_arg;	/**< @brief Gene variance directory (containing .gexpvar files) (default='NA').  */
  char * dir_gvar_orig;	/**< @brief Gene variance directory (containing .gexpvar files) original value given at command line.  */
  const char *dir_gvar_help; /**< @brief Gene variance directory (containing .gexpvar files) help description.  */
  char * quant_arg;	/**< @brief Quant file (assuming all datasets use the same quantization).  */
  char * quant_orig;	/**< @brief Quant file (assuming all datasets use the same quantization) original value given at command line.  */
  const char *quant_help; /**< @brief Quant file (assuming all datasets use the same quantization) help description.  */
  int num_db_arg;	/**< @brief Number of databaselets in database (default='1000').  */
  char * num_db_orig;	/**< @brief Number of databaselets in database original value given at command line.  */
  const char *num_db_help; /**< @brief Number of databaselets in database help description.  */
  int num_threads_arg;	/**< @brief Number of threads (default='8').  */
  char * num_threads_orig;	/**< @brief Number of threads original value given at command line.  */
  const char *num_threads_help; /**< @brief Number of threads help description.  */
  float per_g_required_arg;	/**< @brief Fraction (max 1.0) of genome required to be present in a dataset. Datasets not meeting this requirement are skipped. (default='0').  */
  char * per_g_required_orig;	/**< @brief Fraction (max 1.0) of genome required to be present in a dataset. Datasets not meeting this requirement are skipped. original value given at command line.  */
  const char *per_g_required_help; /**< @brief Fraction (max 1.0) of genome required to be present in a dataset. Datasets not meeting this requirement are skipped. help description.  */
  int neg_cor_flag;	/**< @brief Rank genes and datasets by negative correlations (default=off).  */
  const char *neg_cor_help; /**< @brief Rank genes and datasets by negative correlations help description.  */
  int check_dset_size_flag;	/**< @brief Check the number of samples per dataset (if on, needs --dset_size_file) (default=off).  */
  const char *check_dset_size_help; /**< @brief Check the number of samples per dataset (if on, needs --dset_size_file) help description.  */
  char * weighting_method_arg;	/**< @brief Weighting method: query cross-validated weighting (CV), equal weighting (EQUAL), order statistics weighting (ORDER_STAT), variance weighting (VAR), user-given weighting (USER), SPELL weighting (AVERAGE_Z), user cross-validated weighting (CV_CUSTOM) (default='CV').  */
  char * weighting_method_orig;	/**< @brief Weighting method: query cross-validated weighting (CV), equal weighting (EQUAL), order statistics weighting (ORDER_STAT), variance weighting (VAR), user-given weighting (USER), SPELL weighting (AVERAGE_Z), user cross-validated weighting (CV_CUSTOM) original value given at command line.  */
  const char *weighting_method_help; /**< @brief Weighting method: query cross-validated weighting (CV), equal weighting (EQUAL), order statistics weighting (ORDER_STAT), variance weighting (VAR), user-given weighting (USER), SPELL weighting (AVERAGE_Z), user cross-validated weighting (CV_CUSTOM) help description.  */
  char * user_gene_list_arg;	/**< @brief List of user's gene-sets to which queries will be judged upon (text, 1 line per query).  */
  char * user_gene_list_orig;	/**< @brief List of user's gene-sets to which queries will be judged upon (text, 1 line per query) original value given at command line.  */
  const char *user_gene_list_help; /**< @brief List of user's gene-sets to which queries will be judged upon (text, 1 line per query) help description.  */
  char * user_weight_list_arg;	/**< @brief List of pre-computed dataset weight files (.dweight).  */
  char * user_weight_list_orig;	/**< @brief List of pre-computed dataset weight files (.dweight) original value given at command line.  */
  const char *user_weight_list_help; /**< @brief List of pre-computed dataset weight files (.dweight) help description.  */
  int random_flag;	/**< @brief Generate random ranking score (default=off).  */
  const char *random_help; /**< @brief Generate random ranking score help description.  */
  int num_random_arg;	/**< @brief Number of repetitions of generating random rankings (default='10').  */
  char * num_random_orig;	/**< @brief Number of repetitions of generating random rankings original value given at command line.  */
  const char *num_random_help; /**< @brief Number of repetitions of generating random rankings help description.  */
  char * dist_measure_arg;	/**< @brief Distance measure (default='z_score').  */
  char * dist_measure_orig;	/**< @brief Distance measure original value given at command line.  */
  const char *dist_measure_help; /**< @brief Distance measure help description.  */
  int norm_subavg_flag;	/**< @brief If z_score is selected, subtract each gene's average z-score in the dataset. (default=off).  */
  const char *norm_subavg_help; /**< @brief If z_score is selected, subtract each gene's average z-score in the dataset. help description.  */
  int norm_subavg_plat_flag;	/**< @brief If z_score is selected, subtract each query gene's average score across platforms and divide by its stdev. Performed after --norm_subavg. (default=off).  */
  const char *norm_subavg_plat_help; /**< @brief If z_score is selected, subtract each query gene's average score across platforms and divide by its stdev. Performed after --norm_subavg. help description.  */
  float score_cutoff_arg;	/**< @brief Cutoff on the gene-gene distance score before adding, default: no cutoff (default='-9999').  */
  char * score_cutoff_orig;	/**< @brief Cutoff on the gene-gene distance score before adding, default: no cutoff original value given at command line.  */
  const char *score_cutoff_help; /**< @brief Cutoff on the gene-gene distance score before adding, default: no cutoff help description.  */
  int square_z_flag;	/**< @brief If z_score is selected, square the z-scores. Usually used in conjunction with --score-cutoff=1.0. (default=off).  */
  const char *square_z_help; /**< @brief If z_score is selected, square the z-scores. Usually used in conjunction with --score-cutoff=1.0. help description.  */
  float per_q_required_arg;	/**< @brief Fraction (max 1.0) of query genes required to be present in a dataset, if --score_cutoff is NOT specified. Datasets not meeting this requirement are skipped. However, if --score_cutoff is specified, this is the fraction of query genes required to pass the correlation cutoff, indicated by --score_cutoff (shortened SC) (note: Advanced usage!). In this case, in scoring a gene g, we first calculate the number of query genes (QS) to which g's correlation exceeds SC. Then if the fraction (QS) / (number of query genes) < per_q_required, do not score this gene in dataset's coexpressed gene-list. Use this with caution if --score_cutoff is specified. (default='0').  */
  char * per_q_required_orig;	/**< @brief Fraction (max 1.0) of query genes required to be present in a dataset, if --score_cutoff is NOT specified. Datasets not meeting this requirement are skipped. However, if --score_cutoff is specified, this is the fraction of query genes required to pass the correlation cutoff, indicated by --score_cutoff (shortened SC) (note: Advanced usage!). In this case, in scoring a gene g, we first calculate the number of query genes (QS) to which g's correlation exceeds SC. Then if the fraction (QS) / (number of query genes) < per_q_required, do not score this gene in dataset's coexpressed gene-list. Use this with caution if --score_cutoff is specified. original value given at command line.  */
  const char *per_q_required_help; /**< @brief Fraction (max 1.0) of query genes required to be present in a dataset, if --score_cutoff is NOT specified. Datasets not meeting this requirement are skipped. However, if --score_cutoff is specified, this is the fraction of query genes required to pass the correlation cutoff, indicated by --score_cutoff (shortened SC) (note: Advanced usage!). In this case, in scoring a gene g, we first calculate the number of query genes (QS) to which g's correlation exceeds SC. Then if the fraction (QS) / (number of query genes) < per_q_required, do not score this gene in dataset's coexpressed gene-list. Use this with caution if --score_cutoff is specified. help description.  */
  char * CV_partition_arg;	/**< @brief The query partitioning method (for CV weighting): Leave-One-In, Leave-One-Out, X-Fold. (default='LOI').  */
  char * CV_partition_orig;	/**< @brief The query partitioning method (for CV weighting): Leave-One-In, Leave-One-Out, X-Fold. original value given at command line.  */
  const char *CV_partition_help; /**< @brief The query partitioning method (for CV weighting): Leave-One-In, Leave-One-Out, X-Fold. help description.  */
  int CV_fold_arg;	/**< @brief The number of folds (for X-fold partitioning). (default='5').  */
  char * CV_fold_orig;	/**< @brief The number of folds (for X-fold partitioning). original value given at command line.  */
  const char *CV_fold_help; /**< @brief The number of folds (for X-fold partitioning). help description.  */
  float CV_rbp_p_arg;	/**< @brief The parameter p for RBP scoring of each partition for its query gene retrieval (for CV weighting). (default='0.99').  */
  char * CV_rbp_p_orig;	/**< @brief The parameter p for RBP scoring of each partition for its query gene retrieval (for CV weighting). original value given at command line.  */
  const char *CV_rbp_p_help; /**< @brief The parameter p for RBP scoring of each partition for its query gene retrieval (for CV weighting). help description.  */
  char * dset_size_file_arg;	/**< @brief Dataset size file (required if --check_dset_size) (default='NA').  */
  char * dset_size_file_orig;	/**< @brief Dataset size file (required if --check_dset_size) original value given at command line.  */
  const char *dset_size_file_help; /**< @brief Dataset size file (required if --check_dset_size) help description.  */
  int is_nibble_flag;	/**< @brief The input CDatabase collection is nibble type (default=off).  */
  const char *is_nibble_help; /**< @brief The input CDatabase collection is nibble type help description.  */
  int buffer_arg;	/**< @brief Number of query genes to load in memory (recommended: 50-100) (default='50').  */
  char * buffer_orig;	/**< @brief Number of query genes to load in memory (recommended: 50-100) original value given at command line.  */
  const char *buffer_help; /**< @brief Number of query genes to load in memory (recommended: 50-100) help description.  */
  int output_text_flag;	/**< @brief Output results (gene scores and dataset weights) as text (default=off).  */
  const char *output_text_help; /**< @brief Output results (gene scores and dataset weights) as text help description.  */
  char * output_dir_arg;	/**< @brief Output directory.  */
  char * output_dir_orig;	/**< @brief Output directory original value given at command line.  */
  const char *output_dir_help; /**< @brief Output directory help description.  */
  int output_w_comp_flag;	/**< @brief Output dataset weight components (generates .dweight_comp file) (default=off).  */
  const char *output_w_comp_help; /**< @brief Output dataset weight components (generates .dweight_comp file) help description.  */
  int simulate_w_flag;	/**< @brief If equal weighting or order-statistics weighting is selected, output simulated dataset weights (default=off).  */
  const char *simulate_w_help; /**< @brief If equal weighting or order-statistics weighting is selected, output simulated dataset weights help description.  */
  char * additional_db_arg;	/**< @brief Utilize a second CDatabase collection. Path to the second CDatabase's setting file. (default='NA').  */
  char * additional_db_orig;	/**< @brief Utilize a second CDatabase collection. Path to the second CDatabase's setting file. original value given at command line.  */
  const char *additional_db_help; /**< @brief Utilize a second CDatabase collection. Path to the second CDatabase's setting file. help description.  */
  
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
  unsigned int num_threads_given ;	/**< @brief Whether num_threads was given.  */
  unsigned int per_g_required_given ;	/**< @brief Whether per_g_required was given.  */
  unsigned int neg_cor_given ;	/**< @brief Whether neg_cor was given.  */
  unsigned int check_dset_size_given ;	/**< @brief Whether check_dset_size was given.  */
  unsigned int weighting_method_given ;	/**< @brief Whether weighting_method was given.  */
  unsigned int user_gene_list_given ;	/**< @brief Whether user_gene_list was given.  */
  unsigned int user_weight_list_given ;	/**< @brief Whether user_weight_list was given.  */
  unsigned int random_given ;	/**< @brief Whether random was given.  */
  unsigned int num_random_given ;	/**< @brief Whether num_random was given.  */
  unsigned int dist_measure_given ;	/**< @brief Whether dist_measure was given.  */
  unsigned int norm_subavg_given ;	/**< @brief Whether norm_subavg was given.  */
  unsigned int norm_subavg_plat_given ;	/**< @brief Whether norm_subavg_plat was given.  */
  unsigned int score_cutoff_given ;	/**< @brief Whether score_cutoff was given.  */
  unsigned int square_z_given ;	/**< @brief Whether square_z was given.  */
  unsigned int per_q_required_given ;	/**< @brief Whether per_q_required was given.  */
  unsigned int CV_partition_given ;	/**< @brief Whether CV_partition was given.  */
  unsigned int CV_fold_given ;	/**< @brief Whether CV_fold was given.  */
  unsigned int CV_rbp_p_given ;	/**< @brief Whether CV_rbp_p was given.  */
  unsigned int dset_size_file_given ;	/**< @brief Whether dset_size_file was given.  */
  unsigned int is_nibble_given ;	/**< @brief Whether is_nibble was given.  */
  unsigned int buffer_given ;	/**< @brief Whether buffer was given.  */
  unsigned int output_text_given ;	/**< @brief Whether output_text was given.  */
  unsigned int output_dir_given ;	/**< @brief Whether output_dir was given.  */
  unsigned int output_w_comp_given ;	/**< @brief Whether output_w_comp was given.  */
  unsigned int simulate_w_given ;	/**< @brief Whether simulate_w was given.  */
  unsigned int additional_db_given ;	/**< @brief Whether additional_db was given.  */

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

extern const char *cmdline_parser_weighting_method_values[];  /**< @brief Possible values for weighting_method. */
extern const char *cmdline_parser_dist_measure_values[];  /**< @brief Possible values for dist_measure. */
extern const char *cmdline_parser_CV_partition_values[];  /**< @brief Possible values for CV_partition. */


#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /* CMDLINE_H */
