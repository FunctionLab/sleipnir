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
#define CMDLINE_PARSER_PACKAGE "Dat2Graph"
#endif

#ifndef CMDLINE_PARSER_VERSION
/** @brief the program version */
#define CMDLINE_PARSER_VERSION "1.0"
#endif

/** @brief Where the command line options are stored */
struct gengetopt_args_info {
    const char *help_help; /**< @brief Print help and exit help description.  */
    const char *version_help; /**< @brief Print version and exit help description.  */
    char *input_arg;    /**< @brief Input DAT/DAB file.  */
    char *input_orig;    /**< @brief Input DAT/DAB file original value given at command line.  */
    const char *input_help; /**< @brief Input DAT/DAB file help description.  */
    char *format_arg;    /**< @brief Output graph format (default='dot').  */
    char *format_orig;    /**< @brief Output graph format original value given at command line.  */
    const char *format_help; /**< @brief Output graph format help description.  */
    char *geneq_arg;    /**< @brief Query inclusion file.  */
    char *geneq_orig;    /**< @brief Query inclusion file original value given at command line.  */
    const char *geneq_help; /**< @brief Query inclusion file help description.  */
    char *genew_arg;    /**< @brief Query weights file.  */
    char *genew_orig;    /**< @brief Query weights file original value given at command line.  */
    const char *genew_help; /**< @brief Query weights file help description.  */
    int neighbors_arg;    /**< @brief Size of query neighborhood (default='-1').  */
    char *neighbors_orig;    /**< @brief Size of query neighborhood original value given at command line.  */
    const char *neighbors_help; /**< @brief Size of query neighborhood help description.  */
    int hefalmp_flag;    /**< @brief Perform HEFalMp query instead of bioPIXIE query (default=on).  */
    const char *hefalmp_help; /**< @brief Perform HEFalMp query instead of bioPIXIE query help description.  */
    double edges_arg;    /**< @brief Aggressiveness of edge trimming after query (default='1').  */
    char *edges_orig;    /**< @brief Aggressiveness of edge trimming after query original value given at command line.  */
    const char *edges_help; /**< @brief Aggressiveness of edge trimming after query help description.  */
    int hubs_arg;    /**< @brief Number of neighbors to query hubs (default='-1').  */
    char *hubs_orig;    /**< @brief Number of neighbors to query hubs original value given at command line.  */
    const char *hubs_help; /**< @brief Number of neighbors to query hubs help description.  */
    double cutoff_arg;    /**< @brief Minimum edge weight for output.  */
    char *cutoff_orig;    /**< @brief Minimum edge weight for output original value given at command line.  */
    const char *cutoff_help; /**< @brief Minimum edge weight for output help description.  */
    char *genes_arg;    /**< @brief Gene inclusion file.  */
    char *genes_orig;    /**< @brief Gene inclusion file original value given at command line.  */
    const char *genes_help; /**< @brief Gene inclusion file help description.  */
    char *genex_arg;    /**< @brief Gene exclusion file.  */
    char *genex_orig;    /**< @brief Gene exclusion file original value given at command line.  */
    const char *genex_help; /**< @brief Gene exclusion file help description.  */
    char *knowns_arg;    /**< @brief Known interactions (DAT/DAB) to ignore.  */
    char *knowns_orig;    /**< @brief Known interactions (DAT/DAB) to ignore original value given at command line.  */
    const char *knowns_help; /**< @brief Known interactions (DAT/DAB) to ignore help description.  */
    char *features_arg;    /**< @brief SGD gene features.  */
    char *features_orig;    /**< @brief SGD gene features original value given at command line.  */
    const char *features_help; /**< @brief SGD gene features help description.  */
    char *colors_arg;    /**< @brief Colors for graph nodes.  */
    char *colors_orig;    /**< @brief Colors for graph nodes original value given at command line.  */
    const char *colors_help; /**< @brief Colors for graph nodes help description.  */
    char *borders_arg;    /**< @brief Borders for graph nodes.  */
    char *borders_orig;    /**< @brief Borders for graph nodes original value given at command line.  */
    const char *borders_help; /**< @brief Borders for graph nodes help description.  */
    int normalize_flag;    /**< @brief Normalize edge weights before processing (default=off).  */
    const char *normalize_help; /**< @brief Normalize edge weights before processing help description.  */
    int absolute_flag;    /**< @brief Use absolute value of edge weights (default=off).  */
    const char *absolute_help; /**< @brief Use absolute value of edge weights help description.  */
    int memmap_flag;    /**< @brief Memory map input file (default=off).  */
    const char *memmap_help; /**< @brief Memory map input file help description.  */
    char *config_arg;    /**< @brief Command line config file (default='Dat2Graph.ini').  */
    char *config_orig;    /**< @brief Command line config file original value given at command line.  */
    const char *config_help; /**< @brief Command line config file help description.  */
    int verbosity_arg;    /**< @brief Message verbosity (default='5').  */
    char *verbosity_orig;    /**< @brief Message verbosity original value given at command line.  */
    const char *verbosity_help; /**< @brief Message verbosity help description.  */

    unsigned int help_given;    /**< @brief Whether help was given.  */
    unsigned int version_given;    /**< @brief Whether version was given.  */
    unsigned int input_given;    /**< @brief Whether input was given.  */
    unsigned int format_given;    /**< @brief Whether format was given.  */
    unsigned int geneq_given;    /**< @brief Whether geneq was given.  */
    unsigned int genew_given;    /**< @brief Whether genew was given.  */
    unsigned int neighbors_given;    /**< @brief Whether neighbors was given.  */
    unsigned int hefalmp_given;    /**< @brief Whether hefalmp was given.  */
    unsigned int edges_given;    /**< @brief Whether edges was given.  */
    unsigned int hubs_given;    /**< @brief Whether hubs was given.  */
    unsigned int cutoff_given;    /**< @brief Whether cutoff was given.  */
    unsigned int genes_given;    /**< @brief Whether genes was given.  */
    unsigned int genex_given;    /**< @brief Whether genex was given.  */
    unsigned int knowns_given;    /**< @brief Whether knowns was given.  */
    unsigned int features_given;    /**< @brief Whether features was given.  */
    unsigned int colors_given;    /**< @brief Whether colors was given.  */
    unsigned int borders_given;    /**< @brief Whether borders was given.  */
    unsigned int normalize_given;    /**< @brief Whether normalize was given.  */
    unsigned int absolute_given;    /**< @brief Whether absolute was given.  */
    unsigned int memmap_given;    /**< @brief Whether memmap was given.  */
    unsigned int config_given;    /**< @brief Whether config was given.  */
    unsigned int verbosity_given;    /**< @brief Whether verbosity was given.  */

};

/** @brief The additional parameters to pass to parser functions */
struct cmdline_parser_params {
    int override; /**< @brief whether to override possibly already present options (default 0) */
    int initialize; /**< @brief whether to initialize the option structure gengetopt_args_info (default 1) */
    int check_required; /**< @brief whether to check that all required options were provided (default 1) */
    int check_ambiguity; /**< @brief whether to check for options already specified in the option structure gengetopt_args_info (default 0) */
    int print_errors; /**< @brief whether getopt_long should print an error message for a bad option (default 1) */
};

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
int cmdline_parser(int argc, char *const *argv,
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
int cmdline_parser2(int argc, char *const *argv,
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
int cmdline_parser_ext(int argc, char *const *argv,
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
void cmdline_parser_init(struct gengetopt_args_info *args_info);

/**
 * Deallocates the string fields of the gengetopt_args_info structure
 * (but does not deallocate the structure itself)
 * @param args_info the structure to deallocate
 */
void cmdline_parser_free(struct gengetopt_args_info *args_info);

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
int cmdline_parser_configfile(char *const filename,
                              struct gengetopt_args_info *args_info,
                              int override, int initialize, int check_required);

/**
 * The config file parser
 * @param filename the name of the config file
 * @param args_info the structure where option information will be stored
 * @param params additional parameters for the parser
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser_config_file(char *const filename,
                               struct gengetopt_args_info *args_info,
                               struct cmdline_parser_params *params);

/**
 * Checks that all the required options were specified
 * @param args_info the structure to check
 * @param prog_name the name of the program that will be used to print
 *   possible errors
 * @return
 */
int cmdline_parser_required(struct gengetopt_args_info *args_info,
                            const char *prog_name);

extern char *cmdline_parser_format_values[];    /**< @brief Possible values for format.  */


#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /* CMDLINE_H */
