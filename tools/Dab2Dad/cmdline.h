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
#define CMDLINE_PARSER_PACKAGE "Dab2Dad"
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
  char * input_arg;	/**< @brief Input DAD file.  */
  char * input_orig;	/**< @brief Input DAD file original value given at command line.  */
  const char *input_help; /**< @brief Input DAD file help description.  */
  char * load_arg;	/**< @brief Persistent load DAD file.  */
  char * load_orig;	/**< @brief Persistent load DAD file original value given at command line.  */
  const char *load_help; /**< @brief Persistent load DAD file help description.  */
  char * network_arg;	/**< @brief Input Bayesian network (X)DSL.  */
  char * network_orig;	/**< @brief Input Bayesian network (X)DSL original value given at command line.  */
  const char *network_help; /**< @brief Input Bayesian network (X)DSL help description.  */
  char * output_arg;	/**< @brief Output DAD file.  */
  char * output_orig;	/**< @brief Output DAD file original value given at command line.  */
  const char *output_help; /**< @brief Output DAD file help description.  */
  char * answers_arg;	/**< @brief Answer DAT/DAB file.  */
  char * answers_orig;	/**< @brief Answer DAT/DAB file original value given at command line.  */
  const char *answers_help; /**< @brief Answer DAT/DAB file help description.  */
  char * directory_arg;	/**< @brief Directory with DAB files (default='.').  */
  char * directory_orig;	/**< @brief Directory with DAB files original value given at command line.  */
  const char *directory_help; /**< @brief Directory with DAB files help description.  */
  int everything_flag;	/**< @brief Include pairs without answers (default=off).  */
  const char *everything_help; /**< @brief Include pairs without answers help description.  */
  char * genes_arg;	/**< @brief Gene inclusion file.  */
  char * genes_orig;	/**< @brief Gene inclusion file original value given at command line.  */
  const char *genes_help; /**< @brief Gene inclusion file help description.  */
  char * genex_arg;	/**< @brief Gene exclusion file.  */
  char * genex_orig;	/**< @brief Gene exclusion file original value given at command line.  */
  const char *genex_help; /**< @brief Gene exclusion file help description.  */
  char * lookup1_arg;	/**< @brief First lookup gene.  */
  char * lookup1_orig;	/**< @brief First lookup gene original value given at command line.  */
  const char *lookup1_help; /**< @brief First lookup gene help description.  */
  char * lookup2_arg;	/**< @brief Second lookup gene.  */
  char * lookup2_orig;	/**< @brief Second lookup gene original value given at command line.  */
  const char *lookup2_help; /**< @brief Second lookup gene help description.  */
  char * lookups_arg;	/**< @brief Lookup gene set.  */
  char * lookups_orig;	/**< @brief Lookup gene set original value given at command line.  */
  const char *lookups_help; /**< @brief Lookup gene set help description.  */
  char * lookupp_arg;	/**< @brief Lookup pair set.  */
  char * lookupp_orig;	/**< @brief Lookup pair set original value given at command line.  */
  const char *lookupp_help; /**< @brief Lookup pair set help description.  */
  int quantize_flag;	/**< @brief Discretize lookups (default=off).  */
  const char *quantize_help; /**< @brief Discretize lookups help description.  */
  char * mask_arg;	/**< @brief Mask DAT/DAB file.  */
  char * mask_orig;	/**< @brief Mask DAT/DAB file original value given at command line.  */
  const char *mask_help; /**< @brief Mask DAT/DAB file help description.  */
  int memmap_flag;	/**< @brief Memory map input/output (default=off).  */
  const char *memmap_help; /**< @brief Memory map input/output help description.  */
  int verbosity_arg;	/**< @brief Message verbosity (default='5').  */
  char * verbosity_orig;	/**< @brief Message verbosity original value given at command line.  */
  const char *verbosity_help; /**< @brief Message verbosity help description.  */
  
  unsigned int help_given ;	/**< @brief Whether help was given.  */
  unsigned int version_given ;	/**< @brief Whether version was given.  */
  unsigned int input_given ;	/**< @brief Whether input was given.  */
  unsigned int load_given ;	/**< @brief Whether load was given.  */
  unsigned int network_given ;	/**< @brief Whether network was given.  */
  unsigned int output_given ;	/**< @brief Whether output was given.  */
  unsigned int answers_given ;	/**< @brief Whether answers was given.  */
  unsigned int directory_given ;	/**< @brief Whether directory was given.  */
  unsigned int everything_given ;	/**< @brief Whether everything was given.  */
  unsigned int genes_given ;	/**< @brief Whether genes was given.  */
  unsigned int genex_given ;	/**< @brief Whether genex was given.  */
  unsigned int lookup1_given ;	/**< @brief Whether lookup1 was given.  */
  unsigned int lookup2_given ;	/**< @brief Whether lookup2 was given.  */
  unsigned int lookups_given ;	/**< @brief Whether lookups was given.  */
  unsigned int lookupp_given ;	/**< @brief Whether lookupp was given.  */
  unsigned int quantize_given ;	/**< @brief Whether quantize was given.  */
  unsigned int mask_given ;	/**< @brief Whether mask was given.  */
  unsigned int memmap_given ;	/**< @brief Whether memmap was given.  */
  unsigned int verbosity_given ;	/**< @brief Whether verbosity was given.  */

  char **inputs ; /**< @brief unamed options (options without names) */
  unsigned inputs_num ; /**< @brief unamed options number */
  int Input_Output_group_counter; /**< @brief Counter for group Input_Output */
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
