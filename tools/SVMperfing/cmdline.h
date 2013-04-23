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
#define CMDLINE_PARSER_PACKAGE "SVMperfer"
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
  char * labels_arg;	/**< @brief Labels file.  */
  char * labels_orig;	/**< @brief Labels file original value given at command line.  */
  const char *labels_help; /**< @brief Labels file help description.  */
  char * output_arg;	/**< @brief Output file .  */
  char * output_orig;	/**< @brief Output file  original value given at command line.  */
  const char *output_help; /**< @brief Output file  help description.  */
  char * directory_arg;	/**< @brief input directory (must only contain input files).  */
  char * directory_orig;	/**< @brief input directory (must only contain input files) original value given at command line.  */
  const char *directory_help; /**< @brief input directory (must only contain input files) help description.  */
  char * model_arg;	/**< @brief input Model file.  */
  char * model_orig;	/**< @brief input Model file original value given at command line.  */
  const char *model_help; /**< @brief input Model file help description.  */
  char * modelPrefix_arg;	/**< @brief input Model prefix.  */
  char * modelPrefix_orig;	/**< @brief input Model prefix original value given at command line.  */
  const char *modelPrefix_help; /**< @brief input Model prefix help description.  */
  int slack_flag;	/**< @brief Use slack rescaling (not implemented for ROC loss) (default=off).  */
  const char *slack_help; /**< @brief Use slack rescaling (not implemented for ROC loss) help description.  */
  int verbosity_arg;	/**< @brief Sets the svm_struct verbosity (default='0').  */
  char * verbosity_orig;	/**< @brief Sets the svm_struct verbosity original value given at command line.  */
  const char *verbosity_help; /**< @brief Sets the svm_struct verbosity help description.  */
  int cross_validation_arg;	/**< @brief Number of cross-validation sets ( arg of 1 will turn off cross-validation ) (default='5').  */
  char * cross_validation_orig;	/**< @brief Number of cross-validation sets ( arg of 1 will turn off cross-validation ) original value given at command line.  */
  const char *cross_validation_help; /**< @brief Number of cross-validation sets ( arg of 1 will turn off cross-validation ) help description.  */
  int error_function_arg;	/**< @brief Sets the loss function for SVM learning: Choice of:
  0\tZero/one loss: 1 if vector of predictions contains error, 0 otherwise.
  1\tF1: 100 minus the F1-score in percent.
  2\tErrorrate: Percentage of errors in prediction vector.
  3\tPrec/Rec Breakeven: 100 minus PRBEP in percent.
  4\tPrec@k: 100 minus precision at k in percent.
  5\tRec@k: 100 minus recall at k in percent.
  10\tROCArea: Percentage of swapped pos/neg pairs (i.e. 100 - ROCArea).\n (default='10').  */
  char * error_function_orig;	/**< @brief Sets the loss function for SVM learning: Choice of:
  0\tZero/one loss: 1 if vector of predictions contains error, 0 otherwise.
  1\tF1: 100 minus the F1-score in percent.
  2\tErrorrate: Percentage of errors in prediction vector.
  3\tPrec/Rec Breakeven: 100 minus PRBEP in percent.
  4\tPrec@k: 100 minus precision at k in percent.
  5\tRec@k: 100 minus recall at k in percent.
  10\tROCArea: Percentage of swapped pos/neg pairs (i.e. 100 - ROCArea).\n original value given at command line.  */
  const char *error_function_help; /**< @brief Sets the loss function for SVM learning: Choice of:
  0\tZero/one loss: 1 if vector of predictions contains error, 0 otherwise.
  1\tF1: 100 minus the F1-score in percent.
  2\tErrorrate: Percentage of errors in prediction vector.
  3\tPrec/Rec Breakeven: 100 minus PRBEP in percent.
  4\tPrec@k: 100 minus precision at k in percent.
  5\tRec@k: 100 minus recall at k in percent.
  10\tROCArea: Percentage of swapped pos/neg pairs (i.e. 100 - ROCArea).\n help description.  */
  float k_value_arg;	/**< @brief Value of k parameter used for Prec@k and Rec@k in (0,1) (default='0.5').  */
  char * k_value_orig;	/**< @brief Value of k parameter used for Prec@k and Rec@k in (0,1) original value given at command line.  */
  const char *k_value_help; /**< @brief Value of k parameter used for Prec@k and Rec@k in (0,1) help description.  */
  float tradeoff_arg;	/**< @brief SVM tradeoff constant C (default='1').  */
  char * tradeoff_orig;	/**< @brief SVM tradeoff constant C original value given at command line.  */
  const char *tradeoff_help; /**< @brief SVM tradeoff constant C help description.  */
  char * allgenes_arg;	/**< @brief Gene list that list all genes to make predictions.  */
  char * allgenes_orig;	/**< @brief Gene list that list all genes to make predictions original value given at command line.  */
  const char *allgenes_help; /**< @brief Gene list that list all genes to make predictions help description.  */
  char * params_arg;	/**< @brief NOT IMPLEMENTED YET: Parameter file.  */
  char * params_orig;	/**< @brief NOT IMPLEMENTED YET: Parameter file original value given at command line.  */
  const char *params_help; /**< @brief NOT IMPLEMENTED YET: Parameter file help description.  */
  int nan2neg_flag;	/**< @brief set missing values(NaN in dab file) from labels file as negative examples (default=off).  */
  const char *nan2neg_help; /**< @brief set missing values(NaN in dab file) from labels file as negative examples help description.  */
  int mmap_flag;	/**< @brief Memory map binary input (default=off).  */
  const char *mmap_help; /**< @brief Memory map binary input help description.  */
  int random_arg;	/**< @brief Seed random generator (default -1 uses current time) (default='-1').  */
  char * random_orig;	/**< @brief Seed random generator (default -1 uses current time) original value given at command line.  */
  const char *random_help; /**< @brief Seed random generator (default -1 uses current time) help description.  */
  char * tgene_arg;	/**< @brief Target gene list, use this gene list as gene holdout cross-validation and also filter labels that only have one gene in given target gene list.  */
  char * tgene_orig;	/**< @brief Target gene list, use this gene list as gene holdout cross-validation and also filter labels that only have one gene in given target gene list original value given at command line.  */
  const char *tgene_help; /**< @brief Target gene list, use this gene list as gene holdout cross-validation and also filter labels that only have one gene in given target gene list help description.  */
  int balance_flag;	/**< @brief DEBUG: check before usage, Balance the training gene ratios (default=off).  */
  const char *balance_help; /**< @brief DEBUG: check before usage, Balance the training gene ratios help description.  */
  float bfactor_arg;	/**< @brief DEBUG: only for < 500, When balancing neg and pos counts exmaples for training what factor to increase. default is 1..  */
  char * bfactor_orig;	/**< @brief DEBUG: only for < 500, When balancing neg and pos counts exmaples for training what factor to increase. default is 1. original value given at command line.  */
  const char *bfactor_help; /**< @brief DEBUG: only for < 500, When balancing neg and pos counts exmaples for training what factor to increase. default is 1. help description.  */
  int prob_flag;	/**< @brief Output prediction values as estimated probablity (Platt method) (default=off).  */
  const char *prob_help; /**< @brief Output prediction values as estimated probablity (Platt method) help description.  */
  int probCross_flag;	/**< @brief Cross-validation setting for output prediction values as estimated probablity (Platt method) (default=off).  */
  const char *probCross_help; /**< @brief Cross-validation setting for output prediction values as estimated probablity (Platt method) help description.  */
  int normalizeZero_flag;	/**< @brief Normalize input data to the range [0, 1] (default=off).  */
  const char *normalizeZero_help; /**< @brief Normalize input data to the range [0, 1] help description.  */
  int normalizeNPone_flag;	/**< @brief Normalize input data to the range [-1, 1] (default=off).  */
  const char *normalizeNPone_help; /**< @brief Normalize input data to the range [-1, 1] help description.  */
  int edgeholdout_flag;	/**< @brief For cross-validation perform edge holdout (Default is gene holdout) (default=off).  */
  const char *edgeholdout_help; /**< @brief For cross-validation perform edge holdout (Default is gene holdout) help description.  */
  int skipSVM_flag;	/**< @brief If given this flag, skip training SVM models when file already exist. Often used when cluster runs timeout/error and need to re-run jobs. (default=off).  */
  const char *skipSVM_help; /**< @brief If given this flag, skip training SVM models when file already exist. Often used when cluster runs timeout/error and need to re-run jobs. help description.  */
  int aggregateMax_flag;	/**< @brief If given this flag, when predicting for all gene pairs with multiple SVM models(bagging) aggregate using the maximum prediction value (Default: average) (default=off).  */
  const char *aggregateMax_help; /**< @brief If given this flag, when predicting for all gene pairs with multiple SVM models(bagging) aggregate using the maximum prediction value (Default: average) help description.  */
  int NoCrossPredict_flag;	/**< @brief Don't use the cross-validated prediction values for gene pairs that have labels in the final output. This flag will basically let SVM models make prediction on pairs that were also used for training. (default=off).  */
  const char *NoCrossPredict_help; /**< @brief Don't use the cross-validated prediction values for gene pairs that have labels in the final output. This flag will basically let SVM models make prediction on pairs that were also used for training. help description.  */
  char * CrossResult_arg;	/**< @brief Cross-validation prediction results, if given when prediction mode these values are replaced into final prediction values.  */
  char * CrossResult_orig;	/**< @brief Cross-validation prediction results, if given when prediction mode these values are replaced into final prediction values original value given at command line.  */
  const char *CrossResult_help; /**< @brief Cross-validation prediction results, if given when prediction mode these values are replaced into final prediction values help description.  */
  char * SampledLabels_arg;	/**< @brief Save the sampled final training labels to this file.  */
  char * SampledLabels_orig;	/**< @brief Save the sampled final training labels to this file original value given at command line.  */
  const char *SampledLabels_help; /**< @brief Save the sampled final training labels to this file help description.  */
  float subsample_arg;	/**< @brief Sample the labels to the following rate.  */
  char * subsample_orig;	/**< @brief Sample the labels to the following rate original value given at command line.  */
  const char *subsample_help; /**< @brief Sample the labels to the following rate help description.  */
  char * OutLabels_arg;	/**< @brief Save the sampled labels to the file and exit.  */
  char * OutLabels_orig;	/**< @brief Save the sampled labels to the file and exit original value given at command line.  */
  const char *OutLabels_help; /**< @brief Save the sampled labels to the file and exit help description.  */
  int onetgene_flag;	/**< @brief Only keep edges from lables that have one gene in the target gene list (default=off).  */
  const char *onetgene_help; /**< @brief Only keep edges from lables that have one gene in the target gene list help description.  */
  float prior_arg;	/**< @brief Randomly sub-sample the negative labels to reach target prior. If cannot reach target prior, set to closest prior..  */
  char * prior_orig;	/**< @brief Randomly sub-sample the negative labels to reach target prior. If cannot reach target prior, set to closest prior. original value given at command line.  */
  const char *prior_help; /**< @brief Randomly sub-sample the negative labels to reach target prior. If cannot reach target prior, set to closest prior. help description.  */
  int savemodel_flag;	/**< @brief Save model to file (default=off).  */
  const char *savemodel_help; /**< @brief Save model to file help description.  */
  float mintrain_arg;	/**< @brief Minimum number of total positive examples to allow training, if not met exit.  */
  char * mintrain_orig;	/**< @brief Minimum number of total positive examples to allow training, if not met exit original value given at command line.  */
  const char *mintrain_help; /**< @brief Minimum number of total positive examples to allow training, if not met exit help description.  */
  char * context_arg;	/**< @brief Context gene list.  */
  char * context_orig;	/**< @brief Context gene list original value given at command line.  */
  const char *context_help; /**< @brief Context gene list help description.  */
  
  unsigned int help_given ;	/**< @brief Whether help was given.  */
  unsigned int version_given ;	/**< @brief Whether version was given.  */
  unsigned int labels_given ;	/**< @brief Whether labels was given.  */
  unsigned int output_given ;	/**< @brief Whether output was given.  */
  unsigned int directory_given ;	/**< @brief Whether directory was given.  */
  unsigned int model_given ;	/**< @brief Whether model was given.  */
  unsigned int modelPrefix_given ;	/**< @brief Whether modelPrefix was given.  */
  unsigned int slack_given ;	/**< @brief Whether slack was given.  */
  unsigned int verbosity_given ;	/**< @brief Whether verbosity was given.  */
  unsigned int cross_validation_given ;	/**< @brief Whether cross_validation was given.  */
  unsigned int error_function_given ;	/**< @brief Whether error_function was given.  */
  unsigned int k_value_given ;	/**< @brief Whether k_value was given.  */
  unsigned int tradeoff_given ;	/**< @brief Whether tradeoff was given.  */
  unsigned int allgenes_given ;	/**< @brief Whether allgenes was given.  */
  unsigned int params_given ;	/**< @brief Whether params was given.  */
  unsigned int nan2neg_given ;	/**< @brief Whether nan2neg was given.  */
  unsigned int mmap_given ;	/**< @brief Whether mmap was given.  */
  unsigned int random_given ;	/**< @brief Whether random was given.  */
  unsigned int tgene_given ;	/**< @brief Whether tgene was given.  */
  unsigned int balance_given ;	/**< @brief Whether balance was given.  */
  unsigned int bfactor_given ;	/**< @brief Whether bfactor was given.  */
  unsigned int prob_given ;	/**< @brief Whether prob was given.  */
  unsigned int probCross_given ;	/**< @brief Whether probCross was given.  */
  unsigned int normalizeZero_given ;	/**< @brief Whether normalizeZero was given.  */
  unsigned int normalizeNPone_given ;	/**< @brief Whether normalizeNPone was given.  */
  unsigned int edgeholdout_given ;	/**< @brief Whether edgeholdout was given.  */
  unsigned int skipSVM_given ;	/**< @brief Whether skipSVM was given.  */
  unsigned int aggregateMax_given ;	/**< @brief Whether aggregateMax was given.  */
  unsigned int NoCrossPredict_given ;	/**< @brief Whether NoCrossPredict was given.  */
  unsigned int CrossResult_given ;	/**< @brief Whether CrossResult was given.  */
  unsigned int SampledLabels_given ;	/**< @brief Whether SampledLabels was given.  */
  unsigned int subsample_given ;	/**< @brief Whether subsample was given.  */
  unsigned int OutLabels_given ;	/**< @brief Whether OutLabels was given.  */
  unsigned int onetgene_given ;	/**< @brief Whether onetgene was given.  */
  unsigned int prior_given ;	/**< @brief Whether prior was given.  */
  unsigned int savemodel_given ;	/**< @brief Whether savemodel was given.  */
  unsigned int mintrain_given ;	/**< @brief Whether mintrain was given.  */
  unsigned int context_given ;	/**< @brief Whether context was given.  */

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
