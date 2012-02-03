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
#define CMDLINE_PARSER_PACKAGE "DChecker"
#endif

#ifndef CMDLINE_PARSER_PACKAGE_NAME
/** @brief the complete program name (used for help and version) */
#define CMDLINE_PARSER_PACKAGE_NAME "DChecker"
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
  char * input_arg;	/**< @brief Similarity DAT/DAB file.  */
  char * input_orig;	/**< @brief Similarity DAT/DAB file original value given at command line.  */
  const char *input_help; /**< @brief Similarity DAT/DAB file help description.  */
  char * answers_arg;	/**< @brief Answer DAT/DAB file.  */
  char * answers_orig;	/**< @brief Answer DAT/DAB file original value given at command line.  */
  const char *answers_help; /**< @brief Answer DAT/DAB file help description.  */
  char * directory_arg;	/**< @brief Output directory (default='.').  */
  char * directory_orig;	/**< @brief Output directory original value given at command line.  */
  const char *directory_help; /**< @brief Output directory help description.  */
  float auc_arg;	/**< @brief Use alternative AUCn calculation (default='0').  */
  char * auc_orig;	/**< @brief Use alternative AUCn calculation original value given at command line.  */
  const char *auc_help; /**< @brief Use alternative AUCn calculation help description.  */
  int randomize_arg;	/**< @brief Calculate specified number of randomized scores (default='0').  */
  char * randomize_orig;	/**< @brief Calculate specified number of randomized scores original value given at command line.  */
  const char *randomize_help; /**< @brief Calculate specified number of randomized scores help description.  */
  int bins_arg;	/**< @brief Bins for quantile sorting (default='1000').  */
  char * bins_orig;	/**< @brief Bins for quantile sorting original value given at command line.  */
  const char *bins_help; /**< @brief Bins for quantile sorting help description.  */
  int finite_flag;	/**< @brief Count finitely many bins (default=off).  */
  const char *finite_help; /**< @brief Count finitely many bins help description.  */
  float min_arg;	/**< @brief Minimum correlation to process (default='0').  */
  char * min_orig;	/**< @brief Minimum correlation to process original value given at command line.  */
  const char *min_help; /**< @brief Minimum correlation to process help description.  */
  float max_arg;	/**< @brief Maximum correlation to process (default='1').  */
  char * max_orig;	/**< @brief Maximum correlation to process original value given at command line.  */
  const char *max_help; /**< @brief Maximum correlation to process help description.  */
  double delta_arg;	/**< @brief Size of correlation bins (default='0.01').  */
  char * delta_orig;	/**< @brief Size of correlation bins original value given at command line.  */
  const char *delta_help; /**< @brief Size of correlation bins help description.  */
  char * genes_arg;	/**< @brief Gene inclusion file.  */
  char * genes_orig;	/**< @brief Gene inclusion file original value given at command line.  */
  const char *genes_help; /**< @brief Gene inclusion file help description.  */
  char * genex_arg;	/**< @brief Gene exclusion file.  */
  char * genex_orig;	/**< @brief Gene exclusion file original value given at command line.  */
  const char *genex_help; /**< @brief Gene exclusion file help description.  */
  char * ubiqg_arg;	/**< @brief Ubiquitous gene file (-j and -J refer to connections to ubiq instead of all bridging pairs).  */
  char * ubiqg_orig;	/**< @brief Ubiquitous gene file (-j and -J refer to connections to ubiq instead of all bridging pairs) original value given at command line.  */
  const char *ubiqg_help; /**< @brief Ubiquitous gene file (-j and -J refer to connections to ubiq instead of all bridging pairs) help description.  */
  char * genet_arg;	/**< @brief Term inclusion file.  */
  char * genet_orig;	/**< @brief Term inclusion file original value given at command line.  */
  const char *genet_help; /**< @brief Term inclusion file help description.  */
  char * genee_arg;	/**< @brief Edge inclusion file.  */
  char * genee_orig;	/**< @brief Edge inclusion file original value given at command line.  */
  const char *genee_help; /**< @brief Edge inclusion file help description.  */
  char * genep_arg;	/**< @brief Gene inclusion file for positives.  */
  char * genep_orig;	/**< @brief Gene inclusion file for positives original value given at command line.  */
  const char *genep_help; /**< @brief Gene inclusion file for positives help description.  */
  int ctxtpos_flag;	/**< @brief Use positive edges between context genes (default=on).  */
  const char *ctxtpos_help; /**< @brief Use positive edges between context genes help description.  */
  int ctxtneg_flag;	/**< @brief Use negative edges between context genes (default=on).  */
  const char *ctxtneg_help; /**< @brief Use negative edges between context genes help description.  */
  int bridgepos_flag;	/**< @brief Use bridging positives between context and non-context genes (default=off).  */
  const char *bridgepos_help; /**< @brief Use bridging positives between context and non-context genes help description.  */
  int bridgeneg_flag;	/**< @brief Use bridging negatives between context and non-context genes (default=on).  */
  const char *bridgeneg_help; /**< @brief Use bridging negatives between context and non-context genes help description.  */
  int outpos_flag;	/**< @brief Use positive edges outside the context (default=off).  */
  const char *outpos_help; /**< @brief Use positive edges outside the context help description.  */
  int outneg_flag;	/**< @brief Use negative edges outside the context (default=off).  */
  const char *outneg_help; /**< @brief Use negative edges outside the context help description.  */
  int normalize_flag;	/**< @brief Normalize scores before processing (default=off).  */
  const char *normalize_help; /**< @brief Normalize scores before processing help description.  */
  int invert_flag;	/**< @brief Invert correlations to distances (default=off).  */
  const char *invert_help; /**< @brief Invert correlations to distances help description.  */
  int sse_flag;	/**< @brief Calculate sum of squared errors (default=off).  */
  const char *sse_help; /**< @brief Calculate sum of squared errors help description.  */
  int memmap_flag;	/**< @brief Memory map input DABs (default=off).  */
  const char *memmap_help; /**< @brief Memory map input DABs help description.  */
  int verbosity_arg;	/**< @brief Message verbosity (default='5').  */
  char * verbosity_orig;	/**< @brief Message verbosity original value given at command line.  */
  const char *verbosity_help; /**< @brief Message verbosity help description.  */
  
  unsigned int help_given ;	/**< @brief Whether help was given.  */
  unsigned int version_given ;	/**< @brief Whether version was given.  */
  unsigned int input_given ;	/**< @brief Whether input was given.  */
  unsigned int answers_given ;	/**< @brief Whether answers was given.  */
  unsigned int directory_given ;	/**< @brief Whether directory was given.  */
  unsigned int auc_given ;	/**< @brief Whether auc was given.  */
  unsigned int randomize_given ;	/**< @brief Whether randomize was given.  */
  unsigned int bins_given ;	/**< @brief Whether bins was given.  */
  unsigned int finite_given ;	/**< @brief Whether finite was given.  */
  unsigned int min_given ;	/**< @brief Whether min was given.  */
  unsigned int max_given ;	/**< @brief Whether max was given.  */
  unsigned int delta_given ;	/**< @brief Whether delta was given.  */
  unsigned int genes_given ;	/**< @brief Whether genes was given.  */
  unsigned int genex_given ;	/**< @brief Whether genex was given.  */
  unsigned int ubiqg_given ;	/**< @brief Whether ubiqg was given.  */
  unsigned int genet_given ;	/**< @brief Whether genet was given.  */
  unsigned int genee_given ;	/**< @brief Whether genee was given.  */
  unsigned int genep_given ;	/**< @brief Whether genep was given.  */
  unsigned int ctxtpos_given ;	/**< @brief Whether ctxtpos was given.  */
  unsigned int ctxtneg_given ;	/**< @brief Whether ctxtneg was given.  */
  unsigned int bridgepos_given ;	/**< @brief Whether bridgepos was given.  */
  unsigned int bridgeneg_given ;	/**< @brief Whether bridgeneg was given.  */
  unsigned int outpos_given ;	/**< @brief Whether outpos was given.  */
  unsigned int outneg_given ;	/**< @brief Whether outneg was given.  */
  unsigned int normalize_given ;	/**< @brief Whether normalize was given.  */
  unsigned int invert_given ;	/**< @brief Whether invert was given.  */
  unsigned int sse_given ;	/**< @brief Whether sse was given.  */
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
