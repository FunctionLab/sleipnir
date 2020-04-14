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
#define CMDLINE_PARSER_PACKAGE "SeekReader"
#endif

#ifndef CMDLINE_PARSER_VERSION
/** @brief the program version */
#define CMDLINE_PARSER_VERSION "1.0"
#endif

/** @brief Where the command line options are stored */
struct gengetopt_args_info {
    const char *help_help; /**< @brief Print help and exit help description.  */
    const char *version_help; /**< @brief Print version and exit help description.  */
    int databaselet_flag;    /**< @brief Display values from databaselet(s) (default=off).  */
    const char *databaselet_help; /**< @brief Display values from databaselet(s) help description.  */
    int dataset_flag;    /**< @brief Check which datasets contain query of interest, based on .gpres file (default=off).  */
    const char *dataset_help; /**< @brief Check which datasets contain query of interest, based on .gpres file help description.  */
    int dataset2_flag;    /**< @brief Read datasets' sinfo files, requires -s, -x (default=off).  */
    const char *dataset2_help; /**< @brief Read datasets' sinfo files, requires -s, -x help description.  */
    int weight_flag;    /**< @brief Test dataset weights (default=off).  */
    const char *weight_help; /**< @brief Test dataset weights help description.  */
    int weight2_flag;    /**< @brief Test dataset weights 2 (default=off).  */
    const char *weight2_help; /**< @brief Test dataset weights 2 help description.  */
    int comp_ranking_flag;    /**< @brief Compare two rankings (*.gscore files) (default=off).  */
    const char *comp_ranking_help; /**< @brief Compare two rankings (*.gscore files) help description.  */
    int convert_aracne_flag;    /**< @brief Convert Aracne output (.txt) to DAB file (default=off).  */
    const char *convert_aracne_help; /**< @brief Convert Aracne output (.txt) to DAB file help description.  */
    int convert_dab_flag;    /**< @brief Convert DAB to matrix (default=off).  */
    const char *convert_dab_help; /**< @brief Convert DAB to matrix help description.  */
    int limit_hub_flag;    /**< @brief Limit genes in the DAB to those that are hubby (default=off).  */
    const char *limit_hub_help; /**< @brief Limit genes in the DAB to those that are hubby help description.  */
    int combine_pcl_flag;    /**< @brief Combine PCL bin files (default=off).  */
    const char *combine_pcl_help; /**< @brief Combine PCL bin files help description.  */
    int increase_gscore_flag;    /**< @brief Increase the gene scores (default=off).  */
    const char *increase_gscore_help; /**< @brief Increase the gene scores help description.  */
    int add_gscore_flag;    /**< @brief Add the gene scores (default=off).  */
    const char *add_gscore_help; /**< @brief Add the gene scores help description.  */
    int read_bin_vec_flag;    /**< @brief Display float vector (default=off).  */
    const char *read_bin_vec_help; /**< @brief Display float vector help description.  */
    char *input_vec_arg;    /**< @brief Input file (default='NA').  */
    char *input_vec_orig;    /**< @brief Input file original value given at command line.  */
    const char *input_vec_help; /**< @brief Input file help description.  */
    char *gscore_file_arg;    /**< @brief Gene score file (input) (default='NA').  */
    char *gscore_file_orig;    /**< @brief Gene score file (input) original value given at command line.  */
    const char *gscore_file_help; /**< @brief Gene score file (input) help description.  */
    char *gscore_file_2_arg;    /**< @brief Gene score file (input 2) (default='NA').  */
    char *gscore_file_2_orig;    /**< @brief Gene score file (input 2) original value given at command line.  */
    const char *gscore_file_2_help; /**< @brief Gene score file (input 2) help description.  */
    char *gscore_output_arg;    /**< @brief Gene score output file (default='NA').  */
    char *gscore_output_orig;    /**< @brief Gene score output file original value given at command line.  */
    const char *gscore_output_help; /**< @brief Gene score output file help description.  */
    char *gscore_list_arg;    /**< @brief Gene score list (default='NA').  */
    char *gscore_list_orig;    /**< @brief Gene score list original value given at command line.  */
    const char *gscore_list_help; /**< @brief Gene score list help description.  */
    char *gscore_dir_arg;    /**< @brief Gene score directory (default='NA').  */
    char *gscore_dir_orig;    /**< @brief Gene score directory original value given at command line.  */
    const char *gscore_dir_help; /**< @brief Gene score directory help description.  */
    char *gscore_output2_arg;    /**< @brief Gene score output file (default='NA').  */
    char *gscore_output2_orig;    /**< @brief Gene score output file original value given at command line.  */
    const char *gscore_output2_help; /**< @brief Gene score output file help description.  */
    char *pcl_list_arg;    /**< @brief File containing a list of pcl bin files (including path) (default='NA').  */
    char *pcl_list_orig;    /**< @brief File containing a list of pcl bin files (including path) original value given at command line.  */
    const char *pcl_list_help; /**< @brief File containing a list of pcl bin files (including path) help description.  */
    int binarize_flag;    /**< @brief Binarize the output matrix (default=off).  */
    const char *binarize_help; /**< @brief Binarize the output matrix help description.  */
    char *output_pcl_arg;    /**< @brief Output file (default='NA').  */
    char *output_pcl_orig;    /**< @brief Output file original value given at command line.  */
    const char *output_pcl_help; /**< @brief Output file help description.  */
    char *dabinput_arg;    /**< @brief DAB input file (default='NA').  */
    char *dabinput_orig;    /**< @brief DAB input file original value given at command line.  */
    const char *dabinput_help; /**< @brief DAB input file help description.  */
    char *hub_dab_output_arg;    /**< @brief DAB output file (default='NA').  */
    char *hub_dab_output_orig;    /**< @brief DAB output file original value given at command line.  */
    const char *hub_dab_output_help; /**< @brief DAB output file help description.  */
    char *aracne_file_arg;    /**< @brief Aracne .txt output file (default='NA').  */
    char *aracne_file_orig;    /**< @brief Aracne .txt output file original value given at command line.  */
    const char *aracne_file_help; /**< @brief Aracne .txt output file help description.  */
    char *output_dab_file_arg;    /**< @brief DAB file (default='NA').  */
    char *output_dab_file_orig;    /**< @brief DAB file original value given at command line.  */
    const char *output_dab_file_help; /**< @brief DAB file help description.  */
    char *dab_file_arg;    /**< @brief DAB file (default='NA').  */
    char *dab_file_orig;    /**< @brief DAB file original value given at command line.  */
    const char *dab_file_help; /**< @brief DAB file help description.  */
    char *output_matrix_arg;    /**< @brief Output matrix filename (default='NA').  */
    char *output_matrix_orig;    /**< @brief Output matrix filename original value given at command line.  */
    const char *output_matrix_help; /**< @brief Output matrix filename help description.  */
    char *dweight_dir_arg;    /**< @brief Dataset weight directory (default='NA').  */
    char *dweight_dir_orig;    /**< @brief Dataset weight directory original value given at command line.  */
    const char *dweight_dir_help; /**< @brief Dataset weight directory help description.  */
    int dweight_num_arg;    /**< @brief Number of .dweight files (default='1000').  */
    char *dweight_num_orig;    /**< @brief Number of .dweight files original value given at command line.  */
    const char *dweight_num_help; /**< @brief Number of .dweight files help description.  */
    char *dweight_map_arg;    /**< @brief Dataset mapping file (default='NA').  */
    char *dweight_map_orig;    /**< @brief Dataset mapping file original value given at command line.  */
    const char *dweight_map_help; /**< @brief Dataset mapping file help description.  */
    char *dweight_test_dir_arg;    /**< @brief Test dataset weight directory (default='NA').  */
    char *dweight_test_dir_orig;    /**< @brief Test dataset weight directory original value given at command line.  */
    const char *dweight_test_dir_help; /**< @brief Test dataset weight directory help description.  */
    int dweight_test_num_arg;    /**< @brief Test number of .dweight files (default='1000').  */
    char *dweight_test_num_orig;    /**< @brief Test number of .dweight files original value given at command line.  */
    const char *dweight_test_num_help; /**< @brief Test number of .dweight files help description.  */
    char *gscore_dir1_arg;    /**< @brief Gene score directory 1 (default='NA').  */
    char *gscore_dir1_orig;    /**< @brief Gene score directory 1 original value given at command line.  */
    const char *gscore_dir1_help; /**< @brief Gene score directory 1 help description.  */
    char *gscore_dir2_arg;    /**< @brief Gene score directory 2 (default='NA').  */
    char *gscore_dir2_orig;    /**< @brief Gene score directory 2 original value given at command line.  */
    const char *gscore_dir2_help; /**< @brief Gene score directory 2 help description.  */
    int gscore_num1_arg;    /**< @brief Number of .gscore files (default='1000').  */
    char *gscore_num1_orig;    /**< @brief Number of .gscore files original value given at command line.  */
    const char *gscore_num1_help; /**< @brief Number of .gscore files help description.  */
    int order_stat_single_gene_query_flag;    /**< @brief Order statistics mode (single-gene query) (default=off).  */
    const char *order_stat_single_gene_query_help; /**< @brief Order statistics mode (single-gene query) help description.  */
    char *db_arg;    /**< @brief Input dataset-platform definition.  */
    char *db_orig;    /**< @brief Input dataset-platform definition original value given at command line.  */
    const char *db_help; /**< @brief Input dataset-platform definition help description.  */
    char *dset_list_arg;    /**< @brief Input a set of datasets.  */
    char *dset_list_orig;    /**< @brief Input a set of datasets original value given at command line.  */
    const char *dset_list_help; /**< @brief Input a set of datasets help description.  */
    char *input_arg;    /**< @brief Input gene mapping.  */
    char *input_orig;    /**< @brief Input gene mapping original value given at command line.  */
    const char *input_help; /**< @brief Input gene mapping help description.  */
    char *single_query_arg;    /**< @brief Query gene list.  */
    char *single_query_orig;    /**< @brief Query gene list original value given at command line.  */
    const char *single_query_help; /**< @brief Query gene list help description.  */
    char *dir_in_arg;    /**< @brief Database directory.  */
    char *dir_in_orig;    /**< @brief Database directory original value given at command line.  */
    const char *dir_in_help; /**< @brief Database directory help description.  */
    char *dir_prep_in_arg;    /**< @brief Prep directory (containing .gavg, .gpres files).  */
    char *dir_prep_in_orig;    /**< @brief Prep directory (containing .gavg, .gpres files) original value given at command line.  */
    const char *dir_prep_in_help; /**< @brief Prep directory (containing .gavg, .gpres files) help description.  */
    char *dir_gvar_in_arg;    /**< @brief Prep directory (containing .gexpvar files) (default='NA').  */
    char *dir_gvar_in_orig;    /**< @brief Prep directory (containing .gexpvar files) original value given at command line.  */
    const char *dir_gvar_in_help; /**< @brief Prep directory (containing .gexpvar files) help description.  */
    char *dir_sinfo_in_arg;    /**< @brief Sinfo directory (containing .sinfo files) (default='NA').  */
    char *dir_sinfo_in_orig;    /**< @brief Sinfo directory (containing .sinfo files) original value given at command line.  */
    const char *dir_sinfo_in_help; /**< @brief Sinfo directory (containing .sinfo files) help description.  */
    int is_nibble_flag;    /**< @brief Whether the input DB is nibble type (default=off).  */
    const char *is_nibble_help; /**< @brief Whether the input DB is nibble type help description.  */
    char *platform_dir_arg;    /**< @brief Platform directory.  */
    char *platform_dir_orig;    /**< @brief Platform directory original value given at command line.  */
    const char *platform_dir_help; /**< @brief Platform directory help description.  */
    float gvar_cutoff_arg;    /**< @brief Query gene's variance in the dataset cutoff (default='-1').  */
    char *gvar_cutoff_orig;    /**< @brief Query gene's variance in the dataset cutoff original value given at command line.  */
    const char *gvar_cutoff_help; /**< @brief Query gene's variance in the dataset cutoff help description.  */
    char *multi_query_arg;    /**< @brief File containing multiple queries (default='NA').  */
    char *multi_query_orig;    /**< @brief File containing multiple queries original value given at command line.  */
    const char *multi_query_help; /**< @brief File containing multiple queries help description.  */
    char *output_file_arg;    /**< @brief Output file (default='NA').  */
    char *output_file_orig;    /**< @brief Output file original value given at command line.  */
    const char *output_file_help; /**< @brief Output file help description.  */

    unsigned int help_given;    /**< @brief Whether help was given.  */
    unsigned int version_given;    /**< @brief Whether version was given.  */
    unsigned int databaselet_given;    /**< @brief Whether databaselet was given.  */
    unsigned int dataset_given;    /**< @brief Whether dataset was given.  */
    unsigned int dataset2_given;    /**< @brief Whether dataset2 was given.  */
    unsigned int weight_given;    /**< @brief Whether weight was given.  */
    unsigned int weight2_given;    /**< @brief Whether weight2 was given.  */
    unsigned int comp_ranking_given;    /**< @brief Whether comp_ranking was given.  */
    unsigned int convert_aracne_given;    /**< @brief Whether convert_aracne was given.  */
    unsigned int convert_dab_given;    /**< @brief Whether convert_dab was given.  */
    unsigned int limit_hub_given;    /**< @brief Whether limit_hub was given.  */
    unsigned int combine_pcl_given;    /**< @brief Whether combine_pcl was given.  */
    unsigned int increase_gscore_given;    /**< @brief Whether increase_gscore was given.  */
    unsigned int add_gscore_given;    /**< @brief Whether add_gscore was given.  */
    unsigned int read_bin_vec_given;    /**< @brief Whether read_bin_vec was given.  */
    unsigned int input_vec_given;    /**< @brief Whether input_vec was given.  */
    unsigned int gscore_file_given;    /**< @brief Whether gscore_file was given.  */
    unsigned int gscore_file_2_given;    /**< @brief Whether gscore_file_2 was given.  */
    unsigned int gscore_output_given;    /**< @brief Whether gscore_output was given.  */
    unsigned int gscore_list_given;    /**< @brief Whether gscore_list was given.  */
    unsigned int gscore_dir_given;    /**< @brief Whether gscore_dir was given.  */
    unsigned int gscore_output2_given;    /**< @brief Whether gscore_output2 was given.  */
    unsigned int pcl_list_given;    /**< @brief Whether pcl_list was given.  */
    unsigned int binarize_given;    /**< @brief Whether binarize was given.  */
    unsigned int output_pcl_given;    /**< @brief Whether output_pcl was given.  */
    unsigned int dabinput_given;    /**< @brief Whether dabinput was given.  */
    unsigned int hub_dab_output_given;    /**< @brief Whether hub_dab_output was given.  */
    unsigned int aracne_file_given;    /**< @brief Whether aracne_file was given.  */
    unsigned int output_dab_file_given;    /**< @brief Whether output_dab_file was given.  */
    unsigned int dab_file_given;    /**< @brief Whether dab_file was given.  */
    unsigned int output_matrix_given;    /**< @brief Whether output_matrix was given.  */
    unsigned int dweight_dir_given;    /**< @brief Whether dweight_dir was given.  */
    unsigned int dweight_num_given;    /**< @brief Whether dweight_num was given.  */
    unsigned int dweight_map_given;    /**< @brief Whether dweight_map was given.  */
    unsigned int dweight_test_dir_given;    /**< @brief Whether dweight_test_dir was given.  */
    unsigned int dweight_test_num_given;    /**< @brief Whether dweight_test_num was given.  */
    unsigned int gscore_dir1_given;    /**< @brief Whether gscore_dir1 was given.  */
    unsigned int gscore_dir2_given;    /**< @brief Whether gscore_dir2 was given.  */
    unsigned int gscore_num1_given;    /**< @brief Whether gscore_num1 was given.  */
    unsigned int order_stat_single_gene_query_given;    /**< @brief Whether order_stat_single_gene_query was given.  */
    unsigned int db_given;    /**< @brief Whether db was given.  */
    unsigned int dset_list_given;    /**< @brief Whether dset_list was given.  */
    unsigned int input_given;    /**< @brief Whether input was given.  */
    unsigned int single_query_given;    /**< @brief Whether single_query was given.  */
    unsigned int dir_in_given;    /**< @brief Whether dir_in was given.  */
    unsigned int dir_prep_in_given;    /**< @brief Whether dir_prep_in was given.  */
    unsigned int dir_gvar_in_given;    /**< @brief Whether dir_gvar_in was given.  */
    unsigned int dir_sinfo_in_given;    /**< @brief Whether dir_sinfo_in was given.  */
    unsigned int is_nibble_given;    /**< @brief Whether is_nibble was given.  */
    unsigned int platform_dir_given;    /**< @brief Whether platform_dir was given.  */
    unsigned int gvar_cutoff_given;    /**< @brief Whether gvar_cutoff was given.  */
    unsigned int multi_query_given;    /**< @brief Whether multi_query was given.  */
    unsigned int output_file_given;    /**< @brief Whether output_file was given.  */

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


#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /* CMDLINE_H */
