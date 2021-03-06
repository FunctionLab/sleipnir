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
#define CMDLINE_PARSER_PACKAGE "Clusterer"
#endif

#ifndef CMDLINE_PARSER_PACKAGE_NAME
/** @brief the complete program name (used for help and version) */
#define CMDLINE_PARSER_PACKAGE_NAME "Clusterer"
#endif

#ifndef CMDLINE_PARSER_VERSION
/** @brief the program version */
#define CMDLINE_PARSER_VERSION "1.0"
#endif

/** @brief Where the command line options are stored */
struct gengetopt_args_info {
    const char *help_help; /**< @brief Print help and exit help description.  */
    const char *version_help; /**< @brief Print version and exit help description.  */
    char *input_arg;    /**< @brief Input PCL/DAB file.  */
    char *input_orig;    /**< @brief Input PCL/DAB file original value given at command line.  */
    const char *input_help; /**< @brief Input PCL/DAB file help description.  */
    char *algorithm_arg;    /**< @brief Clustering algorithm (default='kmeans').  */
    char *algorithm_orig;    /**< @brief Clustering algorithm original value given at command line.  */
    const char *algorithm_help; /**< @brief Clustering algorithm help description.  */
    char *weights_arg;    /**< @brief Input weights file.  */
    char *weights_orig;    /**< @brief Input weights file original value given at command line.  */
    const char *weights_help; /**< @brief Input weights file help description.  */
    char *distance_arg;    /**< @brief Similarity measure (default='pearson').  */
    char *distance_orig;    /**< @brief Similarity measure original value given at command line.  */
    const char *distance_help; /**< @brief Similarity measure help description.  */
    int size_arg;    /**< @brief Number of clusters/minimum cluster size (default='10').  */
    char *size_orig;    /**< @brief Number of clusters/minimum cluster size original value given at command line.  */
    const char *size_help; /**< @brief Number of clusters/minimum cluster size help description.  */
    double diameter_arg;    /**< @brief Maximum cluster diameter (default='0.5').  */
    char *diameter_orig;    /**< @brief Maximum cluster diameter original value given at command line.  */
    const char *diameter_help; /**< @brief Maximum cluster diameter help description.  */
    char *output_arg;    /**< @brief Output DAB file.  */
    char *output_orig;    /**< @brief Output DAB file original value given at command line.  */
    const char *output_help; /**< @brief Output DAB file help description.  */
    double diamineter_arg;    /**< @brief Minimum cluster diameter (default='0').  */
    char *diamineter_orig;    /**< @brief Minimum cluster diameter original value given at command line.  */
    const char *diamineter_help; /**< @brief Minimum cluster diameter help description.  */
    double delta_arg;    /**< @brief Cluster diameter step size (default='0').  */
    char *delta_orig;    /**< @brief Cluster diameter step size original value given at command line.  */
    const char *delta_help; /**< @brief Cluster diameter step size help description.  */
    char *output_info_arg;    /**< @brief Output file for clustering info (membership or summary).  */
    char *output_info_orig;    /**< @brief Output file for clustering info (membership or summary) original value given at command line.  */
    const char *output_info_help; /**< @brief Output file for clustering info (membership or summary) help description.  */
    char *pcl_arg;    /**< @brief PCL input if precalculated DAB provided.  */
    char *pcl_orig;    /**< @brief PCL input if precalculated DAB provided original value given at command line.  */
    const char *pcl_help; /**< @brief PCL input if precalculated DAB provided help description.  */
    int skip_arg;    /**< @brief Columns to skip in input PCL (default='2').  */
    char *skip_orig;    /**< @brief Columns to skip in input PCL original value given at command line.  */
    const char *skip_help; /**< @brief Columns to skip in input PCL help description.  */
    int normalize_flag;    /**< @brief Normalize distances before clustering (default=on).  */
    const char *normalize_help; /**< @brief Normalize distances before clustering help description.  */
    int autocorrelate_flag;    /**< @brief Autocorrelate similarity measures (default=off).  */
    const char *autocorrelate_help; /**< @brief Autocorrelate similarity measures help description.  */
    int summary_flag;    /**< @brief Summarize cluster info (default=off).  */
    const char *summary_help; /**< @brief Summarize cluster info help description.  */
    int pcl_out_flag;    /**< @brief Output PCL and clusters a single PCL (default=off).  */
    const char *pcl_out_help; /**< @brief Output PCL and clusters a single PCL help description.  */
    int random_arg;    /**< @brief Seed random generator (default='0').  */
    char *random_orig;    /**< @brief Seed random generator original value given at command line.  */
    const char *random_help; /**< @brief Seed random generator help description.  */
    int verbosity_arg;    /**< @brief Message verbosity (default='5').  */
    char *verbosity_orig;    /**< @brief Message verbosity original value given at command line.  */
    const char *verbosity_help; /**< @brief Message verbosity help description.  */

    unsigned int help_given;    /**< @brief Whether help was given.  */
    unsigned int version_given;    /**< @brief Whether version was given.  */
    unsigned int input_given;    /**< @brief Whether input was given.  */
    unsigned int algorithm_given;    /**< @brief Whether algorithm was given.  */
    unsigned int weights_given;    /**< @brief Whether weights was given.  */
    unsigned int distance_given;    /**< @brief Whether distance was given.  */
    unsigned int size_given;    /**< @brief Whether size was given.  */
    unsigned int diameter_given;    /**< @brief Whether diameter was given.  */
    unsigned int output_given;    /**< @brief Whether output was given.  */
    unsigned int diamineter_given;    /**< @brief Whether diamineter was given.  */
    unsigned int delta_given;    /**< @brief Whether delta was given.  */
    unsigned int output_info_given;    /**< @brief Whether output_info was given.  */
    unsigned int pcl_given;    /**< @brief Whether pcl was given.  */
    unsigned int skip_given;    /**< @brief Whether skip was given.  */
    unsigned int normalize_given;    /**< @brief Whether normalize was given.  */
    unsigned int autocorrelate_given;    /**< @brief Whether autocorrelate was given.  */
    unsigned int summary_given;    /**< @brief Whether summary was given.  */
    unsigned int pcl_out_given;    /**< @brief Whether pcl_out was given.  */
    unsigned int random_given;    /**< @brief Whether random was given.  */
    unsigned int verbosity_given;    /**< @brief Whether verbosity was given.  */

    char **inputs; /**< @brief unamed options (options without names) */
    unsigned inputs_num; /**< @brief unamed options number */
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
int cmdline_parser(int argc, char **argv,
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
int cmdline_parser2(int argc, char **argv,
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
int cmdline_parser_ext(int argc, char **argv,
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
 * Checks that all the required options were specified
 * @param args_info the structure to check
 * @param prog_name the name of the program that will be used to print
 *   possible errors
 * @return
 */
int cmdline_parser_required(struct gengetopt_args_info *args_info,
                            const char *prog_name);

extern const char *cmdline_parser_algorithm_values[];  /**< @brief Possible values for algorithm. */
extern const char *cmdline_parser_distance_values[];  /**< @brief Possible values for distance. */


#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /* CMDLINE_H */
