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
#define CMDLINE_PARSER_PACKAGE "Counter"
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
  char * answers_arg;	/**< @brief Answer file (-w triggers counts mode).  */
  char * answers_orig;	/**< @brief Answer file (-w triggers counts mode) original value given at command line.  */
  const char *answers_help; /**< @brief Answer file (-w triggers counts mode) help description.  */
  char * counts_arg;	/**< @brief Directory containing count files (-k triggers learning mode).  */
  char * counts_orig;	/**< @brief Directory containing count files (-k triggers learning mode) original value given at command line.  */
  const char *counts_help; /**< @brief Directory containing count files (-k triggers learning mode) help description.  */
  char * networks_arg;	/**< @brief Bayes nets (-n triggers inference mode).  */
  char * networks_orig;	/**< @brief Bayes nets (-n triggers inference mode) original value given at command line.  */
  const char *networks_help; /**< @brief Bayes nets (-n triggers inference mode) help description.  */
  char * output_arg;	/**< @brief Output count directory, Bayes nets, or inferences.  */
  char * output_orig;	/**< @brief Output count directory, Bayes nets, or inferences original value given at command line.  */
  const char *output_help; /**< @brief Output count directory, Bayes nets, or inferences help description.  */
  char * countname_arg;	/**< @brief For learning stage, what should the countname be called if no contexts are used (default: global). (default='global').  */
  char * countname_orig;	/**< @brief For learning stage, what should the countname be called if no contexts are used (default: global). original value given at command line.  */
  const char *countname_help; /**< @brief For learning stage, what should the countname be called if no contexts are used (default: global). help description.  */
  char * directory_arg;	/**< @brief Data directory (default='.').  */
  char * directory_orig;	/**< @brief Data directory original value given at command line.  */
  const char *directory_help; /**< @brief Data directory help description.  */
  char * datasets_arg;	/**< @brief Dataset ID text file.  */
  char * datasets_orig;	/**< @brief Dataset ID text file original value given at command line.  */
  const char *datasets_help; /**< @brief Dataset ID text file help description.  */
  char * genome_arg;	/**< @brief Gene ID text file.  */
  char * genome_orig;	/**< @brief Gene ID text file original value given at command line.  */
  const char *genome_help; /**< @brief Gene ID text file help description.  */
  char * contexts_arg;	/**< @brief Context ID text file.  */
  char * contexts_orig;	/**< @brief Context ID text file original value given at command line.  */
  const char *contexts_help; /**< @brief Context ID text file help description.  */
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
  int weights_flag;	/**< @brief Use weighted context file (default=off).  */
  const char *weights_help; /**< @brief Use weighted context file help description.  */
  int dweightpos_flag;	/**< @brief Use weights to divide positive counts into correponding proportion of positives and negatives (default=off).  */
  const char *dweightpos_help; /**< @brief Use weights to divide positive counts into correponding proportion of positives and negatives help description.  */
  int flipneg_flag;	/**< @brief Flip weights(one minus original) for negative standards (default=on).  */
  const char *flipneg_help; /**< @brief Flip weights(one minus original) for negative standards help description.  */
  int noweightneg_flag;	/**< @brief Use weight one for all negative standards (default=off).  */
  const char *noweightneg_help; /**< @brief Use weight one for all negative standards help description.  */
  char * default_arg;	/**< @brief Count file containing defaults for cases with missing data.  */
  char * default_orig;	/**< @brief Count file containing defaults for cases with missing data original value given at command line.  */
  const char *default_help; /**< @brief Count file containing defaults for cases with missing data help description.  */
  char * zeros_arg;	/**< @brief Read zeroed node IDs/outputs from the given file.  */
  char * zeros_orig;	/**< @brief Read zeroed node IDs/outputs from the given file original value given at command line.  */
  const char *zeros_help; /**< @brief Read zeroed node IDs/outputs from the given file help description.  */
  int genewise_flag;	/**< @brief Evaluate networks assuming genewise contexts (default=off).  */
  const char *genewise_help; /**< @brief Evaluate networks assuming genewise contexts help description.  */
  float pseudocounts_arg;	/**< @brief Effective number of pseudocounts to use (default='-1').  */
  char * pseudocounts_orig;	/**< @brief Effective number of pseudocounts to use original value given at command line.  */
  const char *pseudocounts_help; /**< @brief Effective number of pseudocounts to use help description.  */
  char * alphas_arg;	/**< @brief File containing equivalent sample sizes (alphas) for each node.  */
  char * alphas_orig;	/**< @brief File containing equivalent sample sizes (alphas) for each node original value given at command line.  */
  const char *alphas_help; /**< @brief File containing equivalent sample sizes (alphas) for each node help description.  */
  int regularize_flag;	/**< @brief Automatically regularize based on similarity (default=off).  */
  const char *regularize_help; /**< @brief Automatically regularize based on similarity help description.  */
  char * reggroups_arg;	/**< @brief Automatically regularize based on given groups.  */
  char * reggroups_orig;	/**< @brief Automatically regularize based on given groups original value given at command line.  */
  const char *reggroups_help; /**< @brief Automatically regularize based on given groups help description.  */
  char * temporary_arg;	/**< @brief Directory for temporary files (default='.').  */
  char * temporary_orig;	/**< @brief Directory for temporary files original value given at command line.  */
  const char *temporary_help; /**< @brief Directory for temporary files help description.  */
  int smile_flag;	/**< @brief Output SMILE (X)DSL files rather than minimal networks (default=off).  */
  const char *smile_help; /**< @brief Output SMILE (X)DSL files rather than minimal networks help description.  */
  int xdsl_flag;	/**< @brief Generate XDSL output rather than DSL (default=on).  */
  const char *xdsl_help; /**< @brief Generate XDSL output rather than DSL help description.  */
  int memmap_flag;	/**< @brief Memory map input files (default=off).  */
  const char *memmap_help; /**< @brief Memory map input files help description.  */
  int memmapout_flag;	/**< @brief Memory map output files (only for inference mode) (default=off).  */
  const char *memmapout_help; /**< @brief Memory map output files (only for inference mode) help description.  */
  int threads_arg;	/**< @brief Maximum number of threads to spawn (default='-1').  */
  char * threads_orig;	/**< @brief Maximum number of threads to spawn original value given at command line.  */
  const char *threads_help; /**< @brief Maximum number of threads to spawn help description.  */
  int verbosity_arg;	/**< @brief Message verbosity (default='5').  */
  char * verbosity_orig;	/**< @brief Message verbosity original value given at command line.  */
  const char *verbosity_help; /**< @brief Message verbosity help description.  */
  int logratio_flag;	/**< @brief Output log ratios (instead of posteriors) (default=off).  */
  const char *logratio_help; /**< @brief Output log ratios (instead of posteriors) help description.  */
  
  unsigned int help_given ;	/**< @brief Whether help was given.  */
  unsigned int version_given ;	/**< @brief Whether version was given.  */
  unsigned int answers_given ;	/**< @brief Whether answers was given.  */
  unsigned int counts_given ;	/**< @brief Whether counts was given.  */
  unsigned int networks_given ;	/**< @brief Whether networks was given.  */
  unsigned int output_given ;	/**< @brief Whether output was given.  */
  unsigned int countname_given ;	/**< @brief Whether countname was given.  */
  unsigned int directory_given ;	/**< @brief Whether directory was given.  */
  unsigned int datasets_given ;	/**< @brief Whether datasets was given.  */
  unsigned int genome_given ;	/**< @brief Whether genome was given.  */
  unsigned int contexts_given ;	/**< @brief Whether contexts was given.  */
  unsigned int genes_given ;	/**< @brief Whether genes was given.  */
  unsigned int genex_given ;	/**< @brief Whether genex was given.  */
  unsigned int ubiqg_given ;	/**< @brief Whether ubiqg was given.  */
  unsigned int genet_given ;	/**< @brief Whether genet was given.  */
  unsigned int genee_given ;	/**< @brief Whether genee was given.  */
  unsigned int ctxtpos_given ;	/**< @brief Whether ctxtpos was given.  */
  unsigned int ctxtneg_given ;	/**< @brief Whether ctxtneg was given.  */
  unsigned int bridgepos_given ;	/**< @brief Whether bridgepos was given.  */
  unsigned int bridgeneg_given ;	/**< @brief Whether bridgeneg was given.  */
  unsigned int outpos_given ;	/**< @brief Whether outpos was given.  */
  unsigned int outneg_given ;	/**< @brief Whether outneg was given.  */
  unsigned int weights_given ;	/**< @brief Whether weights was given.  */
  unsigned int dweightpos_given ;	/**< @brief Whether dweightpos was given.  */
  unsigned int flipneg_given ;	/**< @brief Whether flipneg was given.  */
  unsigned int noweightneg_given ;	/**< @brief Whether noweightneg was given.  */
  unsigned int default_given ;	/**< @brief Whether default was given.  */
  unsigned int zeros_given ;	/**< @brief Whether zeros was given.  */
  unsigned int genewise_given ;	/**< @brief Whether genewise was given.  */
  unsigned int pseudocounts_given ;	/**< @brief Whether pseudocounts was given.  */
  unsigned int alphas_given ;	/**< @brief Whether alphas was given.  */
  unsigned int regularize_given ;	/**< @brief Whether regularize was given.  */
  unsigned int reggroups_given ;	/**< @brief Whether reggroups was given.  */
  unsigned int temporary_given ;	/**< @brief Whether temporary was given.  */
  unsigned int smile_given ;	/**< @brief Whether smile was given.  */
  unsigned int xdsl_given ;	/**< @brief Whether xdsl was given.  */
  unsigned int memmap_given ;	/**< @brief Whether memmap was given.  */
  unsigned int memmapout_given ;	/**< @brief Whether memmapout was given.  */
  unsigned int threads_given ;	/**< @brief Whether threads was given.  */
  unsigned int verbosity_given ;	/**< @brief Whether verbosity was given.  */
  unsigned int logratio_given ;	/**< @brief Whether logratio was given.  */

  int Mode_group_counter; /**< @brief Counter for group Mode */
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
