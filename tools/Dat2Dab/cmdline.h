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
#define CMDLINE_PARSER_PACKAGE "Dat2Dab"
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
  char * input_arg;	/**< @brief Input DAT/DAB file.  */
  char * input_orig;	/**< @brief Input DAT/DAB file original value given at command line.  */
  const char *input_help; /**< @brief Input DAT/DAB file help description.  */
  char * output_arg;	/**< @brief Output DAT/DAB file.  */
  char * output_orig;	/**< @brief Output DAT/DAB file original value given at command line.  */
  const char *output_help; /**< @brief Output DAT/DAB file help description.  */
  int flip_flag;	/**< @brief Calculate one minus values (default=off).  */
  const char *flip_help; /**< @brief Calculate one minus values help description.  */
  int normalize_flag;	/**< @brief Normalize to the range [0,1] (default=off).  */
  const char *normalize_help; /**< @brief Normalize to the range [0,1] help description.  */
  int zscore_flag;	/**< @brief Convert values to z-scores (default=off).  */
  const char *zscore_help; /**< @brief Convert values to z-scores help description.  */
  int rank_flag;	/**< @brief Rank transform data (default=off).  */
  const char *rank_help; /**< @brief Rank transform data help description.  */
  int randomize_flag;	/**< @brief Randomize data (default=off).  */
  const char *randomize_help; /**< @brief Randomize data help description.  */
  char * genes_arg;	/**< @brief Process only genes from the given set.  */
  char * genes_orig;	/**< @brief Process only genes from the given set original value given at command line.  */
  const char *genes_help; /**< @brief Process only genes from the given set help description.  */
  char * genex_arg;	/**< @brief Exclude all genes from the given set.  */
  char * genex_orig;	/**< @brief Exclude all genes from the given set original value given at command line.  */
  const char *genex_help; /**< @brief Exclude all genes from the given set help description.  */
  char * edges_arg;	/**< @brief Process only edges from the given DAT/DAB.  */
  char * edges_orig;	/**< @brief Process only edges from the given DAT/DAB original value given at command line.  */
  const char *edges_help; /**< @brief Process only edges from the given DAT/DAB help description.  */
  double cutoff_arg;	/**< @brief Exclude edges below cutoff.  */
  char * cutoff_orig;	/**< @brief Exclude edges below cutoff original value given at command line.  */
  const char *cutoff_help; /**< @brief Exclude edges below cutoff help description.  */
  int zero_flag;	/**< @brief Zero missing values (default=off).  */
  const char *zero_help; /**< @brief Zero missing values help description.  */
  int duplicates_flag;	/**< @brief Allow dissimilar duplicate values (default=off).  */
  const char *duplicates_help; /**< @brief Allow dissimilar duplicate values help description.  */
  float subsample_arg;	/**< @brief Fraction of output to randomly subsample (default='1').  */
  char * subsample_orig;	/**< @brief Fraction of output to randomly subsample original value given at command line.  */
  const char *subsample_help; /**< @brief Fraction of output to randomly subsample help description.  */
  char * lookup1_arg;	/**< @brief First lookup gene.  */
  char * lookup1_orig;	/**< @brief First lookup gene original value given at command line.  */
  const char *lookup1_help; /**< @brief First lookup gene help description.  */
  char * lookup2_arg;	/**< @brief Second lookup gene.  */
  char * lookup2_orig;	/**< @brief Second lookup gene original value given at command line.  */
  const char *lookup2_help; /**< @brief Second lookup gene help description.  */
  char * lookups1_arg;	/**< @brief First lookup gene set.  */
  char * lookups1_orig;	/**< @brief First lookup gene set original value given at command line.  */
  const char *lookups1_help; /**< @brief First lookup gene set help description.  */
  char * lookups2_arg;	/**< @brief First lookup gene set.  */
  char * lookups2_orig;	/**< @brief First lookup gene set original value given at command line.  */
  const char *lookups2_help; /**< @brief First lookup gene set help description.  */
  int genelist_flag;	/**< @brief Only list genes (default=off).  */
  const char *genelist_help; /**< @brief Only list genes help description.  */
  int paircount_flag;	/**< @brief Only count pairs above cutoff (default=off).  */
  const char *paircount_help; /**< @brief Only count pairs above cutoff help description.  */
  char * remap_arg;	/**< @brief Gene name remapping file.  */
  char * remap_orig;	/**< @brief Gene name remapping file original value given at command line.  */
  const char *remap_help; /**< @brief Gene name remapping file help description.  */
  int table_flag;	/**< @brief Produce table formatted output (default=off).  */
  const char *table_help; /**< @brief Produce table formatted output help description.  */
  int skip_arg;	/**< @brief Columns to skip in input PCL (default='2').  */
  char * skip_orig;	/**< @brief Columns to skip in input PCL original value given at command line.  */
  const char *skip_help; /**< @brief Columns to skip in input PCL help description.  */
  int memmap_flag;	/**< @brief Memory map input/output (default=off).  */
  const char *memmap_help; /**< @brief Memory map input/output help description.  */
  int verbosity_arg;	/**< @brief Message verbosity (default='5').  */
  char * verbosity_orig;	/**< @brief Message verbosity original value given at command line.  */
  const char *verbosity_help; /**< @brief Message verbosity help description.  */
  
  unsigned int help_given ;	/**< @brief Whether help was given.  */
  unsigned int version_given ;	/**< @brief Whether version was given.  */
  unsigned int input_given ;	/**< @brief Whether input was given.  */
  unsigned int output_given ;	/**< @brief Whether output was given.  */
  unsigned int flip_given ;	/**< @brief Whether flip was given.  */
  unsigned int normalize_given ;	/**< @brief Whether normalize was given.  */
  unsigned int zscore_given ;	/**< @brief Whether zscore was given.  */
  unsigned int rank_given ;	/**< @brief Whether rank was given.  */
  unsigned int randomize_given ;	/**< @brief Whether randomize was given.  */
  unsigned int genes_given ;	/**< @brief Whether genes was given.  */
  unsigned int genex_given ;	/**< @brief Whether genex was given.  */
  unsigned int edges_given ;	/**< @brief Whether edges was given.  */
  unsigned int cutoff_given ;	/**< @brief Whether cutoff was given.  */
  unsigned int zero_given ;	/**< @brief Whether zero was given.  */
  unsigned int duplicates_given ;	/**< @brief Whether duplicates was given.  */
  unsigned int subsample_given ;	/**< @brief Whether subsample was given.  */
  unsigned int lookup1_given ;	/**< @brief Whether lookup1 was given.  */
  unsigned int lookup2_given ;	/**< @brief Whether lookup2 was given.  */
  unsigned int lookups1_given ;	/**< @brief Whether lookups1 was given.  */
  unsigned int lookups2_given ;	/**< @brief Whether lookups2 was given.  */
  unsigned int genelist_given ;	/**< @brief Whether genelist was given.  */
  unsigned int paircount_given ;	/**< @brief Whether paircount was given.  */
  unsigned int remap_given ;	/**< @brief Whether remap was given.  */
  unsigned int table_given ;	/**< @brief Whether table was given.  */
  unsigned int skip_given ;	/**< @brief Whether skip was given.  */
  unsigned int memmap_given ;	/**< @brief Whether memmap was given.  */
  unsigned int verbosity_given ;	/**< @brief Whether verbosity was given.  */

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
