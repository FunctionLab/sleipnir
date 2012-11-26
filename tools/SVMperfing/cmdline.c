/*
  File autogenerated by gengetopt version 2.22
  generated with the following command:
  /Genomics/ogtr03/cypark/sleipnir/../sleipnir-extlib/gengetopt-2.22/src/gengetopt -iSVMperfing.ggo --default-optional -u -N -e 

  The developers of gengetopt consider the fixed text that goes in all
  gengetopt output files to be in the public domain:
  we make no copyright claims on it.
*/

/* If we use autoconf.  */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "getopt.h"

#include "cmdline.h"

const char *gengetopt_args_info_purpose = "Wrapper for SVM perf";

const char *gengetopt_args_info_usage = "Usage: SVMperfer [OPTIONS]... [FILES]...";

const char *gengetopt_args_info_description = "";

const char *gengetopt_args_info_help[] = {
  "  -h, --help                  Print help and exit",
  "  -V, --version               Print version and exit",
  "\nMain:",
  "  -l, --labels=filename       Labels file",
  "  -o, --output=filename       Output file ",
  "  -d, --directory=directory   input directory (must only contain input files)",
  "  -m, --model=filename        input Model file",
  "  -S, --slack                 Use slack rescaling (not implemented for ROC \n                                loss)  (default=off)",
  "\nOptions:",
  "  -v, --verbosity=INT         Sets the svm_struct verbosity  (default=`0')",
  "  -c, --cross_validation=INT  Number of cross-validation sets ( arg of 1 will \n                                turn off cross-validation )  (default=`5')",
  "  -e, --error_function=INT    Sets the loss function for SVM learning: Choice \n                                of:\n\n                                0\tZero/one loss: 1 if vector of predictions \n                                contains error, 0 otherwise.\n\n                                1\tF1: 100 minus the F1-score in percent.\n\n                                2\tErrorrate: Percentage of errors in \n                                prediction vector.\n\n                                3\tPrec/Rec Breakeven: 100 minus PRBEP in \n                                percent.\n\n                                4\tPrec@k: 100 minus precision at k in percent.\n\n                                5\tRec@k: 100 minus recall at k in percent.\n\n                                10\tROCArea: Percentage of swapped pos/neg \n                                pairs (i.e. 100 - ROCArea).\n                                  (default=`10')",
  "  -k, --k_value=FLOAT         Value of k parameter used for Prec@k and Rec@k in \n                                (0,1)  (default=`0.5')",
  "  -t, --tradeoff=FLOAT        SVM tradeoff constant C  (default=`1')",
  "  -p, --params=filename       NOT IMPLEMENTED YET: Parameter file",
  "  -n, --nan2neg               set missing values(NaN in dab file) from labels \n                                file as negative examples  (default=off)",
  "  -M, --mmap                  Memory map binary input  (default=off)",
  "  -R, --random=INT            Seed random generator (default -1 uses current \n                                time)  (default=`-1')",
  "  -T, --tgene=filename        Target gene list, use this gene list as gene \n                                holdout cross-validation and also filter labels \n                                that only have one gene in given target gene \n                                list",
  "  -b, --balance               Balance the training gene ratios  (default=off)",
  "  -F, --bfactor=FLOAT         DEBUG: only for < 500, When balancing neg and pos \n                                counts exmaples for training what factor to \n                                increase. default is 1.",
  "  -B, --prob                  Output prediction values as estimated probablity \n                                (Platt method)  (default=off)",
  "  -z, --normalize             Normalize to the range [0,1]  (default=off)",
  "\nFiltering:",
  "  -q, --geneq=filename        Only keep edges that have one gene in this given \n                                list",
  "  -P, --prior=FLOAT           Randomly sub-sample the negative labels to reach \n                                target prior. If cannot reach target prior, set \n                                to closest prior.",
  "  -s, --savemodel             Save model to file  (default=off)",
  "  -E, --mintrain=FLOAT        Minimum number of total positive examples to \n                                allow training, if not met exit",
    0
};

typedef enum {ARG_NO
  , ARG_FLAG
  , ARG_STRING
  , ARG_INT
  , ARG_FLOAT
} cmdline_parser_arg_type;

static
void clear_given (struct gengetopt_args_info *args_info);
static
void clear_args (struct gengetopt_args_info *args_info);

static int
cmdline_parser_internal (int argc, char * const *argv, struct gengetopt_args_info *args_info,
                        struct cmdline_parser_params *params, const char *additional_error);

static int
cmdline_parser_required2 (struct gengetopt_args_info *args_info, const char *prog_name, const char *additional_error);

static char *
gengetopt_strdup (const char *s);

static
void clear_given (struct gengetopt_args_info *args_info)
{
  args_info->help_given = 0 ;
  args_info->version_given = 0 ;
  args_info->labels_given = 0 ;
  args_info->output_given = 0 ;
  args_info->directory_given = 0 ;
  args_info->model_given = 0 ;
  args_info->slack_given = 0 ;
  args_info->verbosity_given = 0 ;
  args_info->cross_validation_given = 0 ;
  args_info->error_function_given = 0 ;
  args_info->k_value_given = 0 ;
  args_info->tradeoff_given = 0 ;
  args_info->params_given = 0 ;
  args_info->nan2neg_given = 0 ;
  args_info->mmap_given = 0 ;
  args_info->random_given = 0 ;
  args_info->tgene_given = 0 ;
  args_info->balance_given = 0 ;
  args_info->bfactor_given = 0 ;
  args_info->prob_given = 0 ;
  args_info->normalize_given = 0 ;
  args_info->geneq_given = 0 ;
  args_info->prior_given = 0 ;
  args_info->savemodel_given = 0 ;
  args_info->mintrain_given = 0 ;
}

static
void clear_args (struct gengetopt_args_info *args_info)
{
  args_info->labels_arg = NULL;
  args_info->labels_orig = NULL;
  args_info->output_arg = NULL;
  args_info->output_orig = NULL;
  args_info->directory_arg = NULL;
  args_info->directory_orig = NULL;
  args_info->model_arg = NULL;
  args_info->model_orig = NULL;
  args_info->slack_flag = 0;
  args_info->verbosity_arg = 0;
  args_info->verbosity_orig = NULL;
  args_info->cross_validation_arg = 5;
  args_info->cross_validation_orig = NULL;
  args_info->error_function_arg = 10;
  args_info->error_function_orig = NULL;
  args_info->k_value_arg = 0.5;
  args_info->k_value_orig = NULL;
  args_info->tradeoff_arg = 1;
  args_info->tradeoff_orig = NULL;
  args_info->params_arg = NULL;
  args_info->params_orig = NULL;
  args_info->nan2neg_flag = 0;
  args_info->mmap_flag = 0;
  args_info->random_arg = -1;
  args_info->random_orig = NULL;
  args_info->tgene_arg = NULL;
  args_info->tgene_orig = NULL;
  args_info->balance_flag = 0;
  args_info->bfactor_orig = NULL;
  args_info->prob_flag = 0;
  args_info->normalize_flag = 0;
  args_info->geneq_arg = NULL;
  args_info->geneq_orig = NULL;
  args_info->prior_orig = NULL;
  args_info->savemodel_flag = 0;
  args_info->mintrain_orig = NULL;
  
}

static
void init_args_info(struct gengetopt_args_info *args_info)
{


  args_info->help_help = gengetopt_args_info_help[0] ;
  args_info->version_help = gengetopt_args_info_help[1] ;
  args_info->labels_help = gengetopt_args_info_help[3] ;
  args_info->output_help = gengetopt_args_info_help[4] ;
  args_info->directory_help = gengetopt_args_info_help[5] ;
  args_info->model_help = gengetopt_args_info_help[6] ;
  args_info->slack_help = gengetopt_args_info_help[7] ;
  args_info->verbosity_help = gengetopt_args_info_help[9] ;
  args_info->cross_validation_help = gengetopt_args_info_help[10] ;
  args_info->error_function_help = gengetopt_args_info_help[11] ;
  args_info->k_value_help = gengetopt_args_info_help[12] ;
  args_info->tradeoff_help = gengetopt_args_info_help[13] ;
  args_info->params_help = gengetopt_args_info_help[14] ;
  args_info->nan2neg_help = gengetopt_args_info_help[15] ;
  args_info->mmap_help = gengetopt_args_info_help[16] ;
  args_info->random_help = gengetopt_args_info_help[17] ;
  args_info->tgene_help = gengetopt_args_info_help[18] ;
  args_info->balance_help = gengetopt_args_info_help[19] ;
  args_info->bfactor_help = gengetopt_args_info_help[20] ;
  args_info->prob_help = gengetopt_args_info_help[21] ;
  args_info->normalize_help = gengetopt_args_info_help[22] ;
  args_info->geneq_help = gengetopt_args_info_help[24] ;
  args_info->prior_help = gengetopt_args_info_help[25] ;
  args_info->savemodel_help = gengetopt_args_info_help[26] ;
  args_info->mintrain_help = gengetopt_args_info_help[27] ;
  
}

void
cmdline_parser_print_version (void)
{
  printf ("%s %s\n", CMDLINE_PARSER_PACKAGE, CMDLINE_PARSER_VERSION);
}

static void print_help_common(void) {
  cmdline_parser_print_version ();

  if (strlen(gengetopt_args_info_purpose) > 0)
    printf("\n%s\n", gengetopt_args_info_purpose);

  if (strlen(gengetopt_args_info_usage) > 0)
    printf("\n%s\n", gengetopt_args_info_usage);

  printf("\n");

  if (strlen(gengetopt_args_info_description) > 0)
    printf("%s\n", gengetopt_args_info_description);
}

void
cmdline_parser_print_help (void)
{
  int i = 0;
  print_help_common();
  while (gengetopt_args_info_help[i])
    printf("%s\n", gengetopt_args_info_help[i++]);
}

void
cmdline_parser_init (struct gengetopt_args_info *args_info)
{
  clear_given (args_info);
  clear_args (args_info);
  init_args_info (args_info);

  args_info->inputs = NULL;
  args_info->inputs_num = 0;
}

void
cmdline_parser_params_init(struct cmdline_parser_params *params)
{
  if (params)
    { 
      params->override = 0;
      params->initialize = 1;
      params->check_required = 1;
      params->check_ambiguity = 0;
      params->print_errors = 1;
    }
}

struct cmdline_parser_params *
cmdline_parser_params_create(void)
{
  struct cmdline_parser_params *params = 
    (struct cmdline_parser_params *)malloc(sizeof(struct cmdline_parser_params));
  cmdline_parser_params_init(params);  
  return params;
}

static void
free_string_field (char **s)
{
  if (*s)
    {
      free (*s);
      *s = 0;
    }
}


static void
cmdline_parser_release (struct gengetopt_args_info *args_info)
{
  unsigned int i;
  free_string_field (&(args_info->labels_arg));
  free_string_field (&(args_info->labels_orig));
  free_string_field (&(args_info->output_arg));
  free_string_field (&(args_info->output_orig));
  free_string_field (&(args_info->directory_arg));
  free_string_field (&(args_info->directory_orig));
  free_string_field (&(args_info->model_arg));
  free_string_field (&(args_info->model_orig));
  free_string_field (&(args_info->verbosity_orig));
  free_string_field (&(args_info->cross_validation_orig));
  free_string_field (&(args_info->error_function_orig));
  free_string_field (&(args_info->k_value_orig));
  free_string_field (&(args_info->tradeoff_orig));
  free_string_field (&(args_info->params_arg));
  free_string_field (&(args_info->params_orig));
  free_string_field (&(args_info->random_orig));
  free_string_field (&(args_info->tgene_arg));
  free_string_field (&(args_info->tgene_orig));
  free_string_field (&(args_info->bfactor_orig));
  free_string_field (&(args_info->geneq_arg));
  free_string_field (&(args_info->geneq_orig));
  free_string_field (&(args_info->prior_orig));
  free_string_field (&(args_info->mintrain_orig));
  
  
  for (i = 0; i < args_info->inputs_num; ++i)
    free (args_info->inputs [i]);

  if (args_info->inputs_num)
    free (args_info->inputs);

  clear_given (args_info);
}


static void
write_into_file(FILE *outfile, const char *opt, const char *arg, char *values[])
{
  if (arg) {
    fprintf(outfile, "%s=\"%s\"\n", opt, arg);
  } else {
    fprintf(outfile, "%s\n", opt);
  }
}


int
cmdline_parser_dump(FILE *outfile, struct gengetopt_args_info *args_info)
{
  int i = 0;

  if (!outfile)
    {
      fprintf (stderr, "%s: cannot dump options to stream\n", CMDLINE_PARSER_PACKAGE);
      return EXIT_FAILURE;
    }

  if (args_info->help_given)
    write_into_file(outfile, "help", 0, 0 );
  if (args_info->version_given)
    write_into_file(outfile, "version", 0, 0 );
  if (args_info->labels_given)
    write_into_file(outfile, "labels", args_info->labels_orig, 0);
  if (args_info->output_given)
    write_into_file(outfile, "output", args_info->output_orig, 0);
  if (args_info->directory_given)
    write_into_file(outfile, "directory", args_info->directory_orig, 0);
  if (args_info->model_given)
    write_into_file(outfile, "model", args_info->model_orig, 0);
  if (args_info->slack_given)
    write_into_file(outfile, "slack", 0, 0 );
  if (args_info->verbosity_given)
    write_into_file(outfile, "verbosity", args_info->verbosity_orig, 0);
  if (args_info->cross_validation_given)
    write_into_file(outfile, "cross_validation", args_info->cross_validation_orig, 0);
  if (args_info->error_function_given)
    write_into_file(outfile, "error_function", args_info->error_function_orig, 0);
  if (args_info->k_value_given)
    write_into_file(outfile, "k_value", args_info->k_value_orig, 0);
  if (args_info->tradeoff_given)
    write_into_file(outfile, "tradeoff", args_info->tradeoff_orig, 0);
  if (args_info->params_given)
    write_into_file(outfile, "params", args_info->params_orig, 0);
  if (args_info->nan2neg_given)
    write_into_file(outfile, "nan2neg", 0, 0 );
  if (args_info->mmap_given)
    write_into_file(outfile, "mmap", 0, 0 );
  if (args_info->random_given)
    write_into_file(outfile, "random", args_info->random_orig, 0);
  if (args_info->tgene_given)
    write_into_file(outfile, "tgene", args_info->tgene_orig, 0);
  if (args_info->balance_given)
    write_into_file(outfile, "balance", 0, 0 );
  if (args_info->bfactor_given)
    write_into_file(outfile, "bfactor", args_info->bfactor_orig, 0);
  if (args_info->prob_given)
    write_into_file(outfile, "prob", 0, 0 );
  if (args_info->normalize_given)
    write_into_file(outfile, "normalize", 0, 0 );
  if (args_info->geneq_given)
    write_into_file(outfile, "geneq", args_info->geneq_orig, 0);
  if (args_info->prior_given)
    write_into_file(outfile, "prior", args_info->prior_orig, 0);
  if (args_info->savemodel_given)
    write_into_file(outfile, "savemodel", 0, 0 );
  if (args_info->mintrain_given)
    write_into_file(outfile, "mintrain", args_info->mintrain_orig, 0);
  

  i = EXIT_SUCCESS;
  return i;
}

int
cmdline_parser_file_save(const char *filename, struct gengetopt_args_info *args_info)
{
  FILE *outfile;
  int i = 0;

  outfile = fopen(filename, "w");

  if (!outfile)
    {
      fprintf (stderr, "%s: cannot open file for writing: %s\n", CMDLINE_PARSER_PACKAGE, filename);
      return EXIT_FAILURE;
    }

  i = cmdline_parser_dump(outfile, args_info);
  fclose (outfile);

  return i;
}

void
cmdline_parser_free (struct gengetopt_args_info *args_info)
{
  cmdline_parser_release (args_info);
}

/** @brief replacement of strdup, which is not standard */
char *
gengetopt_strdup (const char *s)
{
  char *result = NULL;
  if (!s)
    return result;

  result = (char*)malloc(strlen(s) + 1);
  if (result == (char*)0)
    return (char*)0;
  strcpy(result, s);
  return result;
}

int
cmdline_parser (int argc, char * const *argv, struct gengetopt_args_info *args_info)
{
  return cmdline_parser2 (argc, argv, args_info, 0, 1, 1);
}

int
cmdline_parser_ext (int argc, char * const *argv, struct gengetopt_args_info *args_info,
                   struct cmdline_parser_params *params)
{
  int result;
  result = cmdline_parser_internal (argc, argv, args_info, params, NULL);

  return result;
}

int
cmdline_parser2 (int argc, char * const *argv, struct gengetopt_args_info *args_info, int override, int initialize, int check_required)
{
  int result;
  struct cmdline_parser_params params;
  
  params.override = override;
  params.initialize = initialize;
  params.check_required = check_required;
  params.check_ambiguity = 0;
  params.print_errors = 1;

  result = cmdline_parser_internal (argc, argv, args_info, &params, NULL);

  return result;
}

int
cmdline_parser_required (struct gengetopt_args_info *args_info, const char *prog_name)
{
  int result = EXIT_SUCCESS;

  if (cmdline_parser_required2(args_info, prog_name, NULL) > 0)
    result = EXIT_FAILURE;

  return result;
}

int
cmdline_parser_required2 (struct gengetopt_args_info *args_info, const char *prog_name, const char *additional_error)
{
  int error = 0;

  /* checks for required options */
  if (! args_info->directory_given)
    {
      fprintf (stderr, "%s: '--directory' ('-d') option required%s\n", prog_name, (additional_error ? additional_error : ""));
      error = 1;
    }
  
  
  /* checks for dependences among options */

  return error;
}


static char *package_name = 0;

/**
 * @brief updates an option
 * @param field the generic pointer to the field to update
 * @param orig_field the pointer to the orig field
 * @param field_given the pointer to the number of occurrence of this option
 * @param prev_given the pointer to the number of occurrence already seen
 * @param value the argument for this option (if null no arg was specified)
 * @param possible_values the possible values for this option (if specified)
 * @param default_value the default value (in case the option only accepts fixed values)
 * @param arg_type the type of this option
 * @param check_ambiguity @see cmdline_parser_params.check_ambiguity
 * @param override @see cmdline_parser_params.override
 * @param no_free whether to free a possible previous value
 * @param multiple_option whether this is a multiple option
 * @param long_opt the corresponding long option
 * @param short_opt the corresponding short option (or '-' if none)
 * @param additional_error possible further error specification
 */
static
int update_arg(void *field, char **orig_field,
               unsigned int *field_given, unsigned int *prev_given, 
               char *value, char *possible_values[], const char *default_value,
               cmdline_parser_arg_type arg_type,
               int check_ambiguity, int override,
               int no_free, int multiple_option,
               const char *long_opt, char short_opt,
               const char *additional_error)
{
  char *stop_char = 0;
  const char *val = value;
  int found;
  char **string_field;

  stop_char = 0;
  found = 0;

  if (!multiple_option && prev_given && (*prev_given || (check_ambiguity && *field_given)))
    {
      if (short_opt != '-')
        fprintf (stderr, "%s: `--%s' (`-%c') option given more than once%s\n", 
               package_name, long_opt, short_opt,
               (additional_error ? additional_error : ""));
      else
        fprintf (stderr, "%s: `--%s' option given more than once%s\n", 
               package_name, long_opt,
               (additional_error ? additional_error : ""));
      return 1; /* failure */
    }

    
  if (field_given && *field_given && ! override)
    return 0;
  if (prev_given)
    (*prev_given)++;
  if (field_given)
    (*field_given)++;
  if (possible_values)
    val = possible_values[found];

  switch(arg_type) {
  case ARG_FLAG:
    *((int *)field) = !*((int *)field);
    break;
  case ARG_INT:
    if (val) *((int *)field) = strtol (val, &stop_char, 0);
    break;
  case ARG_FLOAT:
    if (val) *((float *)field) = (float)strtod (val, &stop_char);
    break;
  case ARG_STRING:
    if (val) {
      string_field = (char **)field;
      if (!no_free && *string_field)
        free (*string_field); /* free previous string */
      *string_field = gengetopt_strdup (val);
    }
    break;
  default:
    break;
  };

  /* check numeric conversion */
  switch(arg_type) {
  case ARG_INT:
  case ARG_FLOAT:
    if (val && !(stop_char && *stop_char == '\0')) {
      fprintf(stderr, "%s: invalid numeric value: %s\n", package_name, val);
      return 1; /* failure */
    }
    break;
  default:
    ;
  };

  /* store the original value */
  switch(arg_type) {
  case ARG_NO:
  case ARG_FLAG:
    break;
  default:
    if (value && orig_field) {
      if (no_free) {
        *orig_field = value;
      } else {
        if (*orig_field)
          free (*orig_field); /* free previous string */
        *orig_field = gengetopt_strdup (value);
      }
    }
  };

  return 0; /* OK */
}


int
cmdline_parser_internal (int argc, char * const *argv, struct gengetopt_args_info *args_info,
                        struct cmdline_parser_params *params, const char *additional_error)
{
  int c;	/* Character of the parsed option.  */

  int error = 0;
  struct gengetopt_args_info local_args_info;
  
  int override;
  int initialize;
  int check_required;
  int check_ambiguity;
  
  package_name = argv[0];
  
  override = params->override;
  initialize = params->initialize;
  check_required = params->check_required;
  check_ambiguity = params->check_ambiguity;

  if (initialize)
    cmdline_parser_init (args_info);

  cmdline_parser_init (&local_args_info);

  optarg = 0;
  optind = 0;
  opterr = params->print_errors;
  optopt = '?';

  while (1)
    {
      int option_index = 0;

      static struct option long_options[] = {
        { "help",	0, NULL, 'h' },
        { "version",	0, NULL, 'V' },
        { "labels",	1, NULL, 'l' },
        { "output",	1, NULL, 'o' },
        { "directory",	1, NULL, 'd' },
        { "model",	1, NULL, 'm' },
        { "slack",	0, NULL, 'S' },
        { "verbosity",	1, NULL, 'v' },
        { "cross_validation",	1, NULL, 'c' },
        { "error_function",	1, NULL, 'e' },
        { "k_value",	1, NULL, 'k' },
        { "tradeoff",	1, NULL, 't' },
        { "params",	1, NULL, 'p' },
        { "nan2neg",	0, NULL, 'n' },
        { "mmap",	0, NULL, 'M' },
        { "random",	1, NULL, 'R' },
        { "tgene",	1, NULL, 'T' },
        { "balance",	0, NULL, 'b' },
        { "bfactor",	1, NULL, 'F' },
        { "prob",	0, NULL, 'B' },
        { "normalize",	0, NULL, 'z' },
        { "geneq",	1, NULL, 'q' },
        { "prior",	1, NULL, 'P' },
        { "savemodel",	0, NULL, 's' },
        { "mintrain",	1, NULL, 'E' },
        { NULL,	0, NULL, 0 }
      };

      c = getopt_long (argc, argv, "hVl:o:d:m:Sv:c:e:k:t:p:nMR:T:bF:Bzq:P:sE:", long_options, &option_index);

      if (c == -1) break;	/* Exit from `while (1)' loop.  */

      switch (c)
        {
        case 'h':	/* Print help and exit.  */
          cmdline_parser_print_help ();
          cmdline_parser_free (&local_args_info);
          exit (EXIT_SUCCESS);

        case 'V':	/* Print version and exit.  */
        
        
          if (update_arg( 0 , 
               0 , &(args_info->version_given),
              &(local_args_info.version_given), optarg, 0, 0, ARG_NO,
              check_ambiguity, override, 0, 0,
              "version", 'V',
              additional_error))
            goto failure;
          cmdline_parser_free (&local_args_info);
          return 0;
        
          break;
        case 'l':	/* Labels file.  */
        
        
          if (update_arg( (void *)&(args_info->labels_arg), 
               &(args_info->labels_orig), &(args_info->labels_given),
              &(local_args_info.labels_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "labels", 'l',
              additional_error))
            goto failure;
        
          break;
        case 'o':	/* Output file .  */
        
        
          if (update_arg( (void *)&(args_info->output_arg), 
               &(args_info->output_orig), &(args_info->output_given),
              &(local_args_info.output_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "output", 'o',
              additional_error))
            goto failure;
        
          break;
        case 'd':	/* input directory (must only contain input files).  */
        
        
          if (update_arg( (void *)&(args_info->directory_arg), 
               &(args_info->directory_orig), &(args_info->directory_given),
              &(local_args_info.directory_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "directory", 'd',
              additional_error))
            goto failure;
        
          break;
        case 'm':	/* input Model file.  */
        
        
          if (update_arg( (void *)&(args_info->model_arg), 
               &(args_info->model_orig), &(args_info->model_given),
              &(local_args_info.model_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "model", 'm',
              additional_error))
            goto failure;
        
          break;
        case 'S':	/* Use slack rescaling (not implemented for ROC loss).  */
        
        
          if (update_arg((void *)&(args_info->slack_flag), 0, &(args_info->slack_given),
              &(local_args_info.slack_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "slack", 'S',
              additional_error))
            goto failure;
        
          break;
        case 'v':	/* Sets the svm_struct verbosity.  */
        
        
          if (update_arg( (void *)&(args_info->verbosity_arg), 
               &(args_info->verbosity_orig), &(args_info->verbosity_given),
              &(local_args_info.verbosity_given), optarg, 0, "0", ARG_INT,
              check_ambiguity, override, 0, 0,
              "verbosity", 'v',
              additional_error))
            goto failure;
        
          break;
        case 'c':	/* Number of cross-validation sets ( arg of 1 will turn off cross-validation ).  */
        
        
          if (update_arg( (void *)&(args_info->cross_validation_arg), 
               &(args_info->cross_validation_orig), &(args_info->cross_validation_given),
              &(local_args_info.cross_validation_given), optarg, 0, "5", ARG_INT,
              check_ambiguity, override, 0, 0,
              "cross_validation", 'c',
              additional_error))
            goto failure;
        
          break;
        case 'e':	/* Sets the loss function for SVM learning: Choice of:
        0\tZero/one loss: 1 if vector of predictions contains error, 0 otherwise.
        1\tF1: 100 minus the F1-score in percent.
        2\tErrorrate: Percentage of errors in prediction vector.
        3\tPrec/Rec Breakeven: 100 minus PRBEP in percent.
        4\tPrec@k: 100 minus precision at k in percent.
        5\tRec@k: 100 minus recall at k in percent.
        10\tROCArea: Percentage of swapped pos/neg pairs (i.e. 100 - ROCArea).\n.  */
        
        
          if (update_arg( (void *)&(args_info->error_function_arg), 
               &(args_info->error_function_orig), &(args_info->error_function_given),
              &(local_args_info.error_function_given), optarg, 0, "10", ARG_INT,
              check_ambiguity, override, 0, 0,
              "error_function", 'e',
              additional_error))
            goto failure;
        
          break;
        case 'k':	/* Value of k parameter used for Prec@k and Rec@k in (0,1).  */
        
        
          if (update_arg( (void *)&(args_info->k_value_arg), 
               &(args_info->k_value_orig), &(args_info->k_value_given),
              &(local_args_info.k_value_given), optarg, 0, "0.5", ARG_FLOAT,
              check_ambiguity, override, 0, 0,
              "k_value", 'k',
              additional_error))
            goto failure;
        
          break;
        case 't':	/* SVM tradeoff constant C.  */
        
        
          if (update_arg( (void *)&(args_info->tradeoff_arg), 
               &(args_info->tradeoff_orig), &(args_info->tradeoff_given),
              &(local_args_info.tradeoff_given), optarg, 0, "1", ARG_FLOAT,
              check_ambiguity, override, 0, 0,
              "tradeoff", 't',
              additional_error))
            goto failure;
        
          break;
        case 'p':	/* NOT IMPLEMENTED YET: Parameter file.  */
        
        
          if (update_arg( (void *)&(args_info->params_arg), 
               &(args_info->params_orig), &(args_info->params_given),
              &(local_args_info.params_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "params", 'p',
              additional_error))
            goto failure;
        
          break;
        case 'n':	/* set missing values(NaN in dab file) from labels file as negative examples.  */
        
        
          if (update_arg((void *)&(args_info->nan2neg_flag), 0, &(args_info->nan2neg_given),
              &(local_args_info.nan2neg_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "nan2neg", 'n',
              additional_error))
            goto failure;
        
          break;
        case 'M':	/* Memory map binary input.  */
        
        
          if (update_arg((void *)&(args_info->mmap_flag), 0, &(args_info->mmap_given),
              &(local_args_info.mmap_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "mmap", 'M',
              additional_error))
            goto failure;
        
          break;
        case 'R':	/* Seed random generator (default -1 uses current time).  */
        
        
          if (update_arg( (void *)&(args_info->random_arg), 
               &(args_info->random_orig), &(args_info->random_given),
              &(local_args_info.random_given), optarg, 0, "-1", ARG_INT,
              check_ambiguity, override, 0, 0,
              "random", 'R',
              additional_error))
            goto failure;
        
          break;
        case 'T':	/* Target gene list, use this gene list as gene holdout cross-validation and also filter labels that only have one gene in given target gene list.  */
        
        
          if (update_arg( (void *)&(args_info->tgene_arg), 
               &(args_info->tgene_orig), &(args_info->tgene_given),
              &(local_args_info.tgene_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "tgene", 'T',
              additional_error))
            goto failure;
        
          break;
        case 'b':	/* Balance the training gene ratios.  */
        
        
          if (update_arg((void *)&(args_info->balance_flag), 0, &(args_info->balance_given),
              &(local_args_info.balance_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "balance", 'b',
              additional_error))
            goto failure;
        
          break;
        case 'F':	/* DEBUG: only for < 500, When balancing neg and pos counts exmaples for training what factor to increase. default is 1..  */
        
        
          if (update_arg( (void *)&(args_info->bfactor_arg), 
               &(args_info->bfactor_orig), &(args_info->bfactor_given),
              &(local_args_info.bfactor_given), optarg, 0, 0, ARG_FLOAT,
              check_ambiguity, override, 0, 0,
              "bfactor", 'F',
              additional_error))
            goto failure;
        
          break;
        case 'B':	/* Output prediction values as estimated probablity (Platt method).  */
        
        
          if (update_arg((void *)&(args_info->prob_flag), 0, &(args_info->prob_given),
              &(local_args_info.prob_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "prob", 'B',
              additional_error))
            goto failure;
        
          break;
        case 'z':	/* Normalize to the range [0,1].  */
        
        
          if (update_arg((void *)&(args_info->normalize_flag), 0, &(args_info->normalize_given),
              &(local_args_info.normalize_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "normalize", 'z',
              additional_error))
            goto failure;
        
          break;
        case 'q':	/* Only keep edges that have one gene in this given list.  */
        
        
          if (update_arg( (void *)&(args_info->geneq_arg), 
               &(args_info->geneq_orig), &(args_info->geneq_given),
              &(local_args_info.geneq_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "geneq", 'q',
              additional_error))
            goto failure;
        
          break;
        case 'P':	/* Randomly sub-sample the negative labels to reach target prior. If cannot reach target prior, set to closest prior..  */
        
        
          if (update_arg( (void *)&(args_info->prior_arg), 
               &(args_info->prior_orig), &(args_info->prior_given),
              &(local_args_info.prior_given), optarg, 0, 0, ARG_FLOAT,
              check_ambiguity, override, 0, 0,
              "prior", 'P',
              additional_error))
            goto failure;
        
          break;
        case 's':	/* Save model to file.  */
        
        
          if (update_arg((void *)&(args_info->savemodel_flag), 0, &(args_info->savemodel_given),
              &(local_args_info.savemodel_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "savemodel", 's',
              additional_error))
            goto failure;
        
          break;
        case 'E':	/* Minimum number of total positive examples to allow training, if not met exit.  */
        
        
          if (update_arg( (void *)&(args_info->mintrain_arg), 
               &(args_info->mintrain_orig), &(args_info->mintrain_given),
              &(local_args_info.mintrain_given), optarg, 0, 0, ARG_FLOAT,
              check_ambiguity, override, 0, 0,
              "mintrain", 'E',
              additional_error))
            goto failure;
        
          break;

        case 0:	/* Long option with no short option */
        case '?':	/* Invalid option.  */
          /* `getopt_long' already printed an error message.  */
          goto failure;

        default:	/* bug: option not considered.  */
          fprintf (stderr, "%s: option unknown: %c%s\n", CMDLINE_PARSER_PACKAGE, c, (additional_error ? additional_error : ""));
          abort ();
        } /* switch */
    } /* while */



  if (check_required)
    {
      error += cmdline_parser_required2 (args_info, argv[0], additional_error);
    }

  cmdline_parser_release (&local_args_info);

  if ( error )
    return (EXIT_FAILURE);

  if (optind < argc)
    {
      int i = 0 ;
      int found_prog_name = 0;
      /* whether program name, i.e., argv[0], is in the remaining args
         (this may happen with some implementations of getopt,
          but surely not with the one included by gengetopt) */

      i = optind;
      while (i < argc)
        if (argv[i++] == argv[0]) {
          found_prog_name = 1;
          break;
        }
      i = 0;

      args_info->inputs_num = argc - optind - found_prog_name;
      args_info->inputs =
        (char **)(malloc ((args_info->inputs_num)*sizeof(char *))) ;
      while (optind < argc)
        if (argv[optind++] != argv[0])
          args_info->inputs[ i++ ] = gengetopt_strdup (argv[optind-1]) ;
    }

  return 0;

failure:
  
  cmdline_parser_release (&local_args_info);
  return (EXIT_FAILURE);
}
