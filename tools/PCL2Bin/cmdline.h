/** @file cmdline.h
 *  @brief The header file for the command line option parser
 *  generated by GNU Gengetopt version 2.22.4
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
#define CMDLINE_PARSER_PACKAGE "PCL2Bin"
#endif

#ifndef CMDLINE_PARSER_PACKAGE_NAME
/** @brief the complete program name (used for help and version) */
#define CMDLINE_PARSER_PACKAGE_NAME "PCL2Bin"
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
  char * input_arg;	/**< @brief Input PCL/BIN file.  */
  char * input_orig;	/**< @brief Input PCL/BIN file original value given at command line.  */
  const char *input_help; /**< @brief Input PCL/BIN file help description.  */
  char * output_arg;	/**< @brief Output PCL/BIN file.  */
  char * output_orig;	/**< @brief Output PCL/BIN file original value given at command line.  */
  const char *output_help; /**< @brief Output PCL/BIN file help description.  */
  char * genes_arg;	/**< @brief Process only genes from the given set.  */
  char * genes_orig;	/**< @brief Process only genes from the given set original value given at command line.  */
  const char *genes_help; /**< @brief Process only genes from the given set help description.  */
  char * genex_arg;	/**< @brief Exclude all genes from the given set.  */
  char * genex_orig;	/**< @brief Exclude all genes from the given set original value given at command line.  */
  const char *genex_help; /**< @brief Exclude all genes from the given set help description.  */
  int Genes_flag;	/**< @brief Only print genes and their indecies (default=off).  */
  const char *Genes_help; /**< @brief Only print genes and their indecies help description.  */
  int zrow_flag;	/**< @brief Normalize Z score row (default=off).  */
  const char *zrow_help; /**< @brief Normalize Z score row help description.  */
  int zcol_flag;	/**< @brief Normalize Z score column (default=off).  */
  const char *zcol_help; /**< @brief Normalize Z score column help description.  */
  int scol_flag;	/**< @brief Subtract mean per column (default=off).  */
  const char *scol_help; /**< @brief Subtract mean per column help description.  */
  int normalize_flag;	/**< @brief Normalize 0 1 range (default=off).  */
  const char *normalize_help; /**< @brief Normalize 0 1 range help description.  */
  int transpose_flag;	/**< @brief transpose dataset (prints text to stdout) (default=off).  */
  const char *transpose_help; /**< @brief transpose dataset (prints text to stdout) help description.  */
  int skip_arg;	/**< @brief Columns to skip in input PCL (default='0').  */
  char * skip_orig;	/**< @brief Columns to skip in input PCL original value given at command line.  */
  const char *skip_help; /**< @brief Columns to skip in input PCL help description.  */
  int verbosity_arg;	/**< @brief Message verbosity (default='5').  */
  char * verbosity_orig;	/**< @brief Message verbosity original value given at command line.  */
  const char *verbosity_help; /**< @brief Message verbosity help description.  */
  int mmap_flag;	/**< @brief Memmap binary input (default=off).  */
  const char *mmap_help; /**< @brief Memmap binary input help description.  */
  
  unsigned int help_given ;	/**< @brief Whether help was given.  */
  unsigned int version_given ;	/**< @brief Whether version was given.  */
  unsigned int input_given ;	/**< @brief Whether input was given.  */
  unsigned int output_given ;	/**< @brief Whether output was given.  */
  unsigned int genes_given ;	/**< @brief Whether genes was given.  */
  unsigned int genex_given ;	/**< @brief Whether genex was given.  */
  unsigned int Genes_given ;	/**< @brief Whether Genes was given.  */
  unsigned int zrow_given ;	/**< @brief Whether zrow was given.  */
  unsigned int zcol_given ;	/**< @brief Whether zcol was given.  */
  unsigned int scol_given ;	/**< @brief Whether scol was given.  */
  unsigned int normalize_given ;	/**< @brief Whether normalize was given.  */
  unsigned int transpose_given ;	/**< @brief Whether transpose was given.  */
  unsigned int skip_given ;	/**< @brief Whether skip was given.  */
  unsigned int verbosity_given ;	/**< @brief Whether verbosity was given.  */
  unsigned int mmap_given ;	/**< @brief Whether mmap was given.  */

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
