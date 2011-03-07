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
#define CMDLINE_PARSER_PACKAGE "BNFunc"
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
  char * input_arg;	/**< @brief Ontology slim file.  */
  char * input_orig;	/**< @brief Ontology slim file original value given at command line.  */
  const char *input_help; /**< @brief Ontology slim file help description.  */
  char * directory_arg;	/**< @brief Output directory (default='.').  */
  char * directory_orig;	/**< @brief Output directory original value given at command line.  */
  const char *directory_help; /**< @brief Output directory help description.  */
  char * output_arg;	/**< @brief Answer file.  */
  char * output_orig;	/**< @brief Answer file original value given at command line.  */
  const char *output_help; /**< @brief Answer file help description.  */
  char * negatives_arg;	/**< @brief Negative slim file.  */
  char * negatives_orig;	/**< @brief Negative slim file original value given at command line.  */
  const char *negatives_help; /**< @brief Negative slim file help description.  */
  char * onto_arg;	/**< @brief ontology (obo file).  */
  char * onto_orig;	/**< @brief ontology (obo file) original value given at command line.  */
  const char *onto_help; /**< @brief ontology (obo file) help description.  */
  char * obo_anno_arg;	/**< @brief Gene annotations that correspond to the OBO ontology for the organism of interest..  */
  char * obo_anno_orig;	/**< @brief Gene annotations that correspond to the OBO ontology for the organism of interest. original value given at command line.  */
  const char *obo_anno_help; /**< @brief Gene annotations that correspond to the OBO ontology for the organism of interest. help description.  */
  char * namespace_arg;	/**< @brief Namespace (the gene ontology namespaces can be abbreviated bp, mf, and cc) (default='').  */
  char * namespace_orig;	/**< @brief Namespace (the gene ontology namespaces can be abbreviated bp, mf, and cc) original value given at command line.  */
  const char *namespace_help; /**< @brief Namespace (the gene ontology namespaces can be abbreviated bp, mf, and cc) help description.  */
  char * kegg_arg;	/**< @brief KEGG ontology.  */
  char * kegg_orig;	/**< @brief KEGG ontology original value given at command line.  */
  const char *kegg_help; /**< @brief KEGG ontology help description.  */
  char * kegg_org_arg;	/**< @brief KEGG organism (default='SCE').  */
  char * kegg_org_orig;	/**< @brief KEGG organism original value given at command line.  */
  const char *kegg_org_help; /**< @brief KEGG organism help description.  */
  char * mips_onto_arg;	/**< @brief MIPS ontology.  */
  char * mips_onto_orig;	/**< @brief MIPS ontology original value given at command line.  */
  const char *mips_onto_help; /**< @brief MIPS ontology help description.  */
  char * mips_anno_arg;	/**< @brief MIPS annotations.  */
  char * mips_anno_orig;	/**< @brief MIPS annotations original value given at command line.  */
  const char *mips_anno_help; /**< @brief MIPS annotations help description.  */
  int synonyms_flag;	/**< @brief Prefer synonym names (default=off).  */
  const char *synonyms_help; /**< @brief Prefer synonym names help description.  */
  int dbids_flag;	/**< @brief Include GO database IDs (default=off).  */
  const char *dbids_help; /**< @brief Include GO database IDs help description.  */
  int allids_flag;	/**< @brief Output all available IDs (default=off).  */
  const char *allids_help; /**< @brief Output all available IDs help description.  */
  double test_arg;	/**< @brief Test fraction (default='0').  */
  char * test_orig;	/**< @brief Test fraction original value given at command line.  */
  const char *test_help; /**< @brief Test fraction help description.  */
  char * sql_arg;	/**< @brief File in which to save SQL tables.  */
  char * sql_orig;	/**< @brief File in which to save SQL tables original value given at command line.  */
  const char *sql_help; /**< @brief File in which to save SQL tables help description.  */
  int nsets_flag;	/**< @brief Generate negative sets for input slim (default=off).  */
  const char *nsets_help; /**< @brief Generate negative sets for input slim help description.  */
  double nsetlap_arg;	/**< @brief P-value of overlap for negative rejection (default='0.05').  */
  char * nsetlap_orig;	/**< @brief P-value of overlap for negative rejection original value given at command line.  */
  const char *nsetlap_help; /**< @brief P-value of overlap for negative rejection help description.  */
  char * config_arg;	/**< @brief Command line config file (default='BNFunc.ini').  */
  char * config_orig;	/**< @brief Command line config file original value given at command line.  */
  const char *config_help; /**< @brief Command line config file help description.  */
  int random_arg;	/**< @brief Seed random generator (default='0').  */
  char * random_orig;	/**< @brief Seed random generator original value given at command line.  */
  const char *random_help; /**< @brief Seed random generator help description.  */
  int verbosity_arg;	/**< @brief Message verbosity (default='5').  */
  char * verbosity_orig;	/**< @brief Message verbosity original value given at command line.  */
  const char *verbosity_help; /**< @brief Message verbosity help description.  */
  char * annotations_arg;	/**< @brief File for propogated annotations.  */
  char * annotations_orig;	/**< @brief File for propogated annotations original value given at command line.  */
  const char *annotations_help; /**< @brief File for propogated annotations help description.  */
  
  unsigned int help_given ;	/**< @brief Whether help was given.  */
  unsigned int version_given ;	/**< @brief Whether version was given.  */
  unsigned int input_given ;	/**< @brief Whether input was given.  */
  unsigned int directory_given ;	/**< @brief Whether directory was given.  */
  unsigned int output_given ;	/**< @brief Whether output was given.  */
  unsigned int negatives_given ;	/**< @brief Whether negatives was given.  */
  unsigned int onto_given ;	/**< @brief Whether onto was given.  */
  unsigned int obo_anno_given ;	/**< @brief Whether obo_anno was given.  */
  unsigned int namespace_given ;	/**< @brief Whether namespace was given.  */
  unsigned int kegg_given ;	/**< @brief Whether kegg was given.  */
  unsigned int kegg_org_given ;	/**< @brief Whether kegg_org was given.  */
  unsigned int mips_onto_given ;	/**< @brief Whether mips_onto was given.  */
  unsigned int mips_anno_given ;	/**< @brief Whether mips_anno was given.  */
  unsigned int synonyms_given ;	/**< @brief Whether synonyms was given.  */
  unsigned int dbids_given ;	/**< @brief Whether dbids was given.  */
  unsigned int allids_given ;	/**< @brief Whether allids was given.  */
  unsigned int test_given ;	/**< @brief Whether test was given.  */
  unsigned int sql_given ;	/**< @brief Whether sql was given.  */
  unsigned int nsets_given ;	/**< @brief Whether nsets was given.  */
  unsigned int nsetlap_given ;	/**< @brief Whether nsetlap was given.  */
  unsigned int config_given ;	/**< @brief Whether config was given.  */
  unsigned int random_given ;	/**< @brief Whether random was given.  */
  unsigned int verbosity_given ;	/**< @brief Whether verbosity was given.  */
  unsigned int annotations_given ;	/**< @brief Whether annotations was given.  */

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
 * The config file parser (deprecated version)
 * @param filename the name of the config file
 * @param args_info the structure where option information will be stored
 * @param override whether to override possibly already present options
 * @param initialize whether to initialize the option structure my_args_info
 * @param check_required whether to check that all required options were provided
 * @return 0 if everything went fine, NON 0 if an error took place
 * @deprecated use cmdline_parser_config_file() instead
 */
int cmdline_parser_configfile (char * const filename,
  struct gengetopt_args_info *args_info,
  int override, int initialize, int check_required);

/**
 * The config file parser
 * @param filename the name of the config file
 * @param args_info the structure where option information will be stored
 * @param params additional parameters for the parser
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser_config_file (char * const filename,
  struct gengetopt_args_info *args_info,
  struct cmdline_parser_params *params);

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
