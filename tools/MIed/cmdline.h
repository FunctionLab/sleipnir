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
#define CMDLINE_PARSER_PACKAGE "MIed"
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
  char * output_arg;	/**< @brief Output file.  */
  char * output_orig;	/**< @brief Output file original value given at command line.  */
  const char *output_help; /**< @brief Output file help description.  */
  int verbosity_arg;	/**< @brief Message verbosity (default='5').  */
  char * verbosity_orig;	/**< @brief Message verbosity original value given at command line.  */
  const char *verbosity_help; /**< @brief Message verbosity help description.  */
  int random_arg;	/**< @brief Seed random generator (default='0').  */
  char * random_orig;	/**< @brief Seed random generator original value given at command line.  */
  const char *random_help; /**< @brief Seed random generator help description.  */
  int start_arg;	/**< @brief Process only the starting at dataset index (default='-1').  */
  char * start_orig;	/**< @brief Process only the starting at dataset index original value given at command line.  */
  const char *start_help; /**< @brief Process only the starting at dataset index help description.  */
  int end_arg;	/**< @brief Process only up to this ending dataset index (default='-1').  */
  char * end_orig;	/**< @brief Process only up to this ending dataset index original value given at command line.  */
  const char *end_help; /**< @brief Process only up to this ending dataset index help description.  */
  int union_flag;	/**< @brief Use the union genes of two datasets while assign missing values randomly (default=off).  */
  const char *union_help; /**< @brief Use the union genes of two datasets while assign missing values randomly help description.  */
  char * directory_arg;	/**< @brief input directory.  */
  char * directory_orig;	/**< @brief input directory original value given at command line.  */
  const char *directory_help; /**< @brief input directory help description.  */
  char * datasets_arg;	/**< @brief Calculate MI for datasets in given file against all datasets.  */
  char * datasets_orig;	/**< @brief Calculate MI for datasets in given file against all datasets original value given at command line.  */
  const char *datasets_help; /**< @brief Calculate MI for datasets in given file against all datasets help description.  */
  char * zeros_arg;	/**< @brief Read zeroed node IDs/outputs from the given file.  */
  char * zeros_orig;	/**< @brief Read zeroed node IDs/outputs from the given file original value given at command line.  */
  const char *zeros_help; /**< @brief Read zeroed node IDs/outputs from the given file help description.  */
  char * edges_arg;	/**< @brief Process only edges from the given DAT/DAB.  */
  char * edges_orig;	/**< @brief Process only edges from the given DAT/DAB original value given at command line.  */
  const char *edges_help; /**< @brief Process only edges from the given DAT/DAB help description.  */
  
  unsigned int help_given ;	/**< @brief Whether help was given.  */
  unsigned int version_given ;	/**< @brief Whether version was given.  */
  unsigned int output_given ;	/**< @brief Whether output was given.  */
  unsigned int verbosity_given ;	/**< @brief Whether verbosity was given.  */
  unsigned int random_given ;	/**< @brief Whether random was given.  */
  unsigned int start_given ;	/**< @brief Whether start was given.  */
  unsigned int end_given ;	/**< @brief Whether end was given.  */
  unsigned int union_given ;	/**< @brief Whether union was given.  */
  unsigned int directory_given ;	/**< @brief Whether directory was given.  */
  unsigned int datasets_given ;	/**< @brief Whether datasets was given.  */
  unsigned int zeros_given ;	/**< @brief Whether zeros was given.  */
  unsigned int edges_given ;	/**< @brief Whether edges was given.  */

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
