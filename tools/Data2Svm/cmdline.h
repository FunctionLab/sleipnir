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
#define CMDLINE_PARSER_PACKAGE "Data2Svm"
#endif

#ifndef CMDLINE_PARSER_VERSION
/** @brief the program version */
#define CMDLINE_PARSER_VERSION "1.0"
#endif

/** @brief Where the command line options are stored */
struct gengetopt_args_info {
    const char *help_help; /**< @brief Print help and exit help description.  */
    const char *version_help; /**< @brief Print version and exit help description.  */
    char *input_arg;    /**< @brief Data set to analyze (PCL).  */
    char *input_orig;    /**< @brief Data set to analyze (PCL) original value given at command line.  */
    const char *input_help; /**< @brief Data set to analyze (PCL) help description.  */
    char *model_arg;    /**< @brief SVM model file.  */
    char *model_orig;    /**< @brief SVM model file original value given at command line.  */
    const char *model_help; /**< @brief SVM model file help description.  */
    char *genes_arg;    /**< @brief List of positive genes.  */
    char *genes_orig;    /**< @brief List of positive genes original value given at command line.  */
    const char *genes_help; /**< @brief List of positive genes help description.  */
    char *genex_arg;    /**< @brief List of test genes.  */
    char *genex_orig;    /**< @brief List of test genes original value given at command line.  */
    const char *genex_help; /**< @brief List of test genes help description.  */
    int heldout_flag;    /**< @brief Evaluate only test genes (default=off).  */
    const char *heldout_help; /**< @brief Evaluate only test genes help description.  */
    int random_features_flag;    /**< @brief Randomize input features (default=off).  */
    const char *random_features_help; /**< @brief Randomize input features help description.  */
    int random_output_flag;    /**< @brief Randomize output values (default=off).  */
    const char *random_output_help; /**< @brief Randomize output values help description.  */
    int cache_arg;    /**< @brief SVM cache size (default='40').  */
    char *cache_orig;    /**< @brief SVM cache size original value given at command line.  */
    const char *cache_help; /**< @brief SVM cache size help description.  */
    char *kernel_arg;    /**< @brief SVM kernel function (default='linear').  */
    char *kernel_orig;    /**< @brief SVM kernel function original value given at command line.  */
    const char *kernel_help; /**< @brief SVM kernel function help description.  */
    float tradeoff_arg;    /**< @brief Classification tradeoff.  */
    char *tradeoff_orig;    /**< @brief Classification tradeoff original value given at command line.  */
    const char *tradeoff_help; /**< @brief Classification tradeoff help description.  */
    float gamma_arg;    /**< @brief RBF gamma (default='1').  */
    char *gamma_orig;    /**< @brief RBF gamma original value given at command line.  */
    const char *gamma_help; /**< @brief RBF gamma help description.  */
    int degree_arg;    /**< @brief Polynomial degree (default='3').  */
    char *degree_orig;    /**< @brief Polynomial degree original value given at command line.  */
    const char *degree_help; /**< @brief Polynomial degree help description.  */
    char *alphas_arg;    /**< @brief SVM alphas file.  */
    char *alphas_orig;    /**< @brief SVM alphas file original value given at command line.  */
    const char *alphas_help; /**< @brief SVM alphas file help description.  */
    int iterations_arg;    /**< @brief SVM iterations (default='100000').  */
    char *iterations_orig;    /**< @brief SVM iterations original value given at command line.  */
    const char *iterations_help; /**< @brief SVM iterations help description.  */
    int normalize_flag;    /**< @brief Z-score normalize feature values (default=off).  */
    const char *normalize_help; /**< @brief Z-score normalize feature values help description.  */
    int skip_arg;    /**< @brief Columns to skip in input PCL (default='2').  */
    char *skip_orig;    /**< @brief Columns to skip in input PCL original value given at command line.  */
    const char *skip_help; /**< @brief Columns to skip in input PCL help description.  */
    int random_arg;    /**< @brief Seed random generator (default='0').  */
    char *random_orig;    /**< @brief Seed random generator original value given at command line.  */
    const char *random_help; /**< @brief Seed random generator help description.  */
    int verbosity_arg;    /**< @brief Message verbosity (default='5').  */
    char *verbosity_orig;    /**< @brief Message verbosity original value given at command line.  */
    const char *verbosity_help; /**< @brief Message verbosity help description.  */

    unsigned int help_given;    /**< @brief Whether help was given.  */
    unsigned int version_given;    /**< @brief Whether version was given.  */
    unsigned int input_given;    /**< @brief Whether input was given.  */
    unsigned int model_given;    /**< @brief Whether model was given.  */
    unsigned int genes_given;    /**< @brief Whether genes was given.  */
    unsigned int genex_given;    /**< @brief Whether genex was given.  */
    unsigned int heldout_given;    /**< @brief Whether heldout was given.  */
    unsigned int random_features_given;    /**< @brief Whether random_features was given.  */
    unsigned int random_output_given;    /**< @brief Whether random_output was given.  */
    unsigned int cache_given;    /**< @brief Whether cache was given.  */
    unsigned int kernel_given;    /**< @brief Whether kernel was given.  */
    unsigned int tradeoff_given;    /**< @brief Whether tradeoff was given.  */
    unsigned int gamma_given;    /**< @brief Whether gamma was given.  */
    unsigned int degree_given;    /**< @brief Whether degree was given.  */
    unsigned int alphas_given;    /**< @brief Whether alphas was given.  */
    unsigned int iterations_given;    /**< @brief Whether iterations was given.  */
    unsigned int normalize_given;    /**< @brief Whether normalize was given.  */
    unsigned int skip_given;    /**< @brief Whether skip was given.  */
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
 * Checks that all the required options were specified
 * @param args_info the structure to check
 * @param prog_name the name of the program that will be used to print
 *   possible errors
 * @return
 */
int cmdline_parser_required(struct gengetopt_args_info *args_info,
                            const char *prog_name);

extern char *cmdline_parser_kernel_values[];    /**< @brief Possible values for kernel.  */


#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /* CMDLINE_H */
