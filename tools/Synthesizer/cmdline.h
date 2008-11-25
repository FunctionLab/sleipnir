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
#define CMDLINE_PARSER_PACKAGE "Synthesizer"
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
  char * output_pcl_arg;	/**< @brief PCL expression output file.  */
  char * output_pcl_orig;	/**< @brief PCL expression output file original value given at command line.  */
  const char *output_pcl_help; /**< @brief PCL expression output file help description.  */
  char * output_fasta_arg;	/**< @brief FASTA sequence output file.  */
  char * output_fasta_orig;	/**< @brief FASTA sequence output file original value given at command line.  */
  const char *output_fasta_help; /**< @brief FASTA sequence output file help description.  */
  int genes_arg;	/**< @brief Number of synthesized genes (default='5000').  */
  char * genes_orig;	/**< @brief Number of synthesized genes original value given at command line.  */
  const char *genes_help; /**< @brief Number of synthesized genes help description.  */
  int conditions_arg;	/**< @brief Number of synthesized conditions (default='100').  */
  char * conditions_orig;	/**< @brief Number of synthesized conditions original value given at command line.  */
  const char *conditions_help; /**< @brief Number of synthesized conditions help description.  */
  int tfs_arg;	/**< @brief Number of transcription factors (default='10').  */
  char * tfs_orig;	/**< @brief Number of transcription factors original value given at command line.  */
  const char *tfs_help; /**< @brief Number of transcription factors help description.  */
  double tf_gene_arg;	/**< @brief Probability of TF activity in a gene (default='0.01').  */
  char * tf_gene_orig;	/**< @brief Probability of TF activity in a gene original value given at command line.  */
  const char *tf_gene_help; /**< @brief Probability of TF activity in a gene help description.  */
  double tf_condition_arg;	/**< @brief Probability of TF activity in a condition (default='0.1').  */
  char * tf_condition_orig;	/**< @brief Probability of TF activity in a condition original value given at command line.  */
  const char *tf_condition_help; /**< @brief Probability of TF activity in a condition help description.  */
  int tf_min_arg;	/**< @brief Minimum TFBS length (default='5').  */
  char * tf_min_orig;	/**< @brief Minimum TFBS length original value given at command line.  */
  const char *tf_min_help; /**< @brief Minimum TFBS length help description.  */
  int tf_max_arg;	/**< @brief Maximum TFBS length (default='12').  */
  char * tf_max_orig;	/**< @brief Maximum TFBS length original value given at command line.  */
  const char *tf_max_help; /**< @brief Maximum TFBS length help description.  */
  double mean_arg;	/**< @brief Expression mean (default='0').  */
  char * mean_orig;	/**< @brief Expression mean original value given at command line.  */
  const char *mean_help; /**< @brief Expression mean help description.  */
  double stdev_arg;	/**< @brief Expression standard deviation (default='1').  */
  char * stdev_orig;	/**< @brief Expression standard deviation original value given at command line.  */
  const char *stdev_help; /**< @brief Expression standard deviation help description.  */
  double tf_mean_arg;	/**< @brief Up/downregulation mean (default='2').  */
  char * tf_mean_orig;	/**< @brief Up/downregulation mean original value given at command line.  */
  const char *tf_mean_help; /**< @brief Up/downregulation mean help description.  */
  double tf_stdev_arg;	/**< @brief Up/downregulation standard deviation (default='1').  */
  char * tf_stdev_orig;	/**< @brief Up/downregulation standard deviation original value given at command line.  */
  const char *tf_stdev_help; /**< @brief Up/downregulation standard deviation help description.  */
  char * fasta_arg;	/**< @brief Input FASTA file.  */
  char * fasta_orig;	/**< @brief Input FASTA file original value given at command line.  */
  const char *fasta_help; /**< @brief Input FASTA file help description.  */
  int degree_arg;	/**< @brief Degree of sequence model HMM (default='3').  */
  char * degree_orig;	/**< @brief Degree of sequence model HMM original value given at command line.  */
  const char *degree_help; /**< @brief Degree of sequence model HMM help description.  */
  int seq_min_arg;	/**< @brief Minimum sequence length (default='1000').  */
  char * seq_min_orig;	/**< @brief Minimum sequence length original value given at command line.  */
  const char *seq_min_help; /**< @brief Minimum sequence length help description.  */
  int seq_max_arg;	/**< @brief Maximum sequence length (default='3000').  */
  char * seq_max_orig;	/**< @brief Maximum sequence length original value given at command line.  */
  const char *seq_max_help; /**< @brief Maximum sequence length help description.  */
  int tf_copm_arg;	/**< @brief Minimum TFBS copies (default='1').  */
  char * tf_copm_orig;	/**< @brief Minimum TFBS copies original value given at command line.  */
  const char *tf_copm_help; /**< @brief Minimum TFBS copies help description.  */
  int tf_copx_arg;	/**< @brief Maximum TFBS copies (default='5').  */
  char * tf_copx_orig;	/**< @brief Maximum TFBS copies original value given at command line.  */
  const char *tf_copx_help; /**< @brief Maximum TFBS copies help description.  */
  char * tf_types_arg;	/**< @brief Sequence types containing TFBSs, comma separated.  */
  char * tf_types_orig;	/**< @brief Sequence types containing TFBSs, comma separated original value given at command line.  */
  const char *tf_types_help; /**< @brief Sequence types containing TFBSs, comma separated help description.  */
  int wrap_arg;	/**< @brief Wrap width for FASTA output (default='60').  */
  char * wrap_orig;	/**< @brief Wrap width for FASTA output original value given at command line.  */
  const char *wrap_help; /**< @brief Wrap width for FASTA output help description.  */
  int random_arg;	/**< @brief Seed random generator (default='0').  */
  char * random_orig;	/**< @brief Seed random generator original value given at command line.  */
  const char *random_help; /**< @brief Seed random generator help description.  */
  int verbosity_arg;	/**< @brief Message verbosity (default='5').  */
  char * verbosity_orig;	/**< @brief Message verbosity original value given at command line.  */
  const char *verbosity_help; /**< @brief Message verbosity help description.  */
  
  unsigned int help_given ;	/**< @brief Whether help was given.  */
  unsigned int version_given ;	/**< @brief Whether version was given.  */
  unsigned int output_pcl_given ;	/**< @brief Whether output_pcl was given.  */
  unsigned int output_fasta_given ;	/**< @brief Whether output_fasta was given.  */
  unsigned int genes_given ;	/**< @brief Whether genes was given.  */
  unsigned int conditions_given ;	/**< @brief Whether conditions was given.  */
  unsigned int tfs_given ;	/**< @brief Whether tfs was given.  */
  unsigned int tf_gene_given ;	/**< @brief Whether tf_gene was given.  */
  unsigned int tf_condition_given ;	/**< @brief Whether tf_condition was given.  */
  unsigned int tf_min_given ;	/**< @brief Whether tf_min was given.  */
  unsigned int tf_max_given ;	/**< @brief Whether tf_max was given.  */
  unsigned int mean_given ;	/**< @brief Whether mean was given.  */
  unsigned int stdev_given ;	/**< @brief Whether stdev was given.  */
  unsigned int tf_mean_given ;	/**< @brief Whether tf_mean was given.  */
  unsigned int tf_stdev_given ;	/**< @brief Whether tf_stdev was given.  */
  unsigned int fasta_given ;	/**< @brief Whether fasta was given.  */
  unsigned int degree_given ;	/**< @brief Whether degree was given.  */
  unsigned int seq_min_given ;	/**< @brief Whether seq_min was given.  */
  unsigned int seq_max_given ;	/**< @brief Whether seq_max was given.  */
  unsigned int tf_copm_given ;	/**< @brief Whether tf_copm was given.  */
  unsigned int tf_copx_given ;	/**< @brief Whether tf_copx was given.  */
  unsigned int tf_types_given ;	/**< @brief Whether tf_types was given.  */
  unsigned int wrap_given ;	/**< @brief Whether wrap was given.  */
  unsigned int random_given ;	/**< @brief Whether random was given.  */
  unsigned int verbosity_given ;	/**< @brief Whether verbosity was given.  */

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
