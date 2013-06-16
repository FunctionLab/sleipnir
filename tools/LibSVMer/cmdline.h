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
#define CMDLINE_PARSER_PACKAGE "LibSVMer"
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
  char * labels_arg;	/**< @brief Labels file.  */
  char * labels_orig;	/**< @brief Labels file original value given at command line.  */
  const char *labels_help; /**< @brief Labels file help description.  */
  char * output_arg;	/**< @brief Output file .  */
  char * output_orig;	/**< @brief Output file  original value given at command line.  */
  const char *output_help; /**< @brief Output file  help description.  */
  char * input_arg;	/**< @brief Input PCL file .  */
  char * input_orig;	/**< @brief Input PCL file  original value given at command line.  */
  const char *input_help; /**< @brief Input PCL file  help description.  */
  char * model_arg;	/**< @brief Model file.  */
  char * model_orig;	/**< @brief Model file original value given at command line.  */
  const char *model_help; /**< @brief Model file help description.  */
  int all_flag;	/**< @brief Always classify all genes in PCLs (default=off).  */
  const char *all_help; /**< @brief Always classify all genes in PCLs help description.  */
  int skip_arg;	/**< @brief Number of columns to skip in input pcls (default='2').  */
  char * skip_orig;	/**< @brief Number of columns to skip in input pcls original value given at command line.  */
  const char *skip_help; /**< @brief Number of columns to skip in input pcls help description.  */
  int normalize_flag;	/**< @brief Normalize PCLS to 0 mean 1 variance (default=off).  */
  const char *normalize_help; /**< @brief Normalize PCLS to 0 mean 1 variance help description.  */
  int cross_validation_arg;	/**< @brief Number of cross-validation sets ( arg of 1 will turn off cross-validation ) (default='5').  */
  char * cross_validation_orig;	/**< @brief Number of cross-validation sets ( arg of 1 will turn off cross-validation ) original value given at command line.  */
  const char *cross_validation_help; /**< @brief Number of cross-validation sets ( arg of 1 will turn off cross-validation ) help description.  */
  int num_cv_runs_arg;	/**< @brief Number of cross-validation runs (default='1').  */
  char * num_cv_runs_orig;	/**< @brief Number of cross-validation runs original value given at command line.  */
  const char *num_cv_runs_help; /**< @brief Number of cross-validation runs help description.  */
  int svm_type_arg;	/**< @brief Sets type of SVM (default 0)
  0\tC-SVC
  1\tnu-SVC
  2\tone-class SVM\n (default='0').  */
  char * svm_type_orig;	/**< @brief Sets type of SVM (default 0)
  0\tC-SVC
  1\tnu-SVC
  2\tone-class SVM\n original value given at command line.  */
  const char *svm_type_help; /**< @brief Sets type of SVM (default 0)
  0\tC-SVC
  1\tnu-SVC
  2\tone-class SVM\n help description.  */
  int balance_flag;	/**< @brief weight classes such that C_P * n_P = C_N * n_N (default=off).  */
  const char *balance_help; /**< @brief weight classes such that C_P * n_P = C_N * n_N help description.  */
  float tradeoff_arg;	/**< @brief SVM tradeoff constant C of C-SVC (default='1').  */
  char * tradeoff_orig;	/**< @brief SVM tradeoff constant C of C-SVC original value given at command line.  */
  const char *tradeoff_help; /**< @brief SVM tradeoff constant C of C-SVC help description.  */
  float nu_arg;	/**< @brief nu parameter of nu-SVC, one-class SVM (default='0.5').  */
  char * nu_orig;	/**< @brief nu parameter of nu-SVC, one-class SVM original value given at command line.  */
  const char *nu_help; /**< @brief nu parameter of nu-SVC, one-class SVM help description.  */
  int mmap_flag;	/**< @brief Memory map binary input (default=off).  */
  const char *mmap_help; /**< @brief Memory map binary input help description.  */
  
  unsigned int help_given ;	/**< @brief Whether help was given.  */
  unsigned int version_given ;	/**< @brief Whether version was given.  */
  unsigned int labels_given ;	/**< @brief Whether labels was given.  */
  unsigned int output_given ;	/**< @brief Whether output was given.  */
  unsigned int input_given ;	/**< @brief Whether input was given.  */
  unsigned int model_given ;	/**< @brief Whether model was given.  */
  unsigned int all_given ;	/**< @brief Whether all was given.  */
  unsigned int skip_given ;	/**< @brief Whether skip was given.  */
  unsigned int normalize_given ;	/**< @brief Whether normalize was given.  */
  unsigned int cross_validation_given ;	/**< @brief Whether cross_validation was given.  */
  unsigned int num_cv_runs_given ;	/**< @brief Whether num_cv_runs was given.  */
  unsigned int svm_type_given ;	/**< @brief Whether svm_type was given.  */
  unsigned int balance_given ;	/**< @brief Whether balance was given.  */
  unsigned int tradeoff_given ;	/**< @brief Whether tradeoff was given.  */
  unsigned int nu_given ;	/**< @brief Whether nu was given.  */
  unsigned int mmap_given ;	/**< @brief Whether mmap was given.  */

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
