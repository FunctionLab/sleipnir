/*
  File autogenerated by gengetopt version 2.22
  generated with the following command:
  /home/chuttenh/hg/sleipnir/trunk/../extlib/gengetopt-2.22/bin/gengetopt -iMEFIT.ggo --default-optional -u -N -e 

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

const char *gengetopt_args_info_purpose = "Microarray Expression Functional Integration Technique (Huttenhower et al,\n\n	Bioinformatics 2006)\n\n\n\nMEFIT takes as input:\n\n1. A collection of microarray data sets (PCL files provided on the command \nline)\n\n2. A collection of known biological functions (lists of related genes provided\n\n   using the -r flag)\n\n3. A collection of known unrelated gene pairs (provided using the -u flag)\n\n\n\nIt produces as output:\n\n1. A global Bayesian network learned by considered all of the data sets\n\n   independently of biological function (specified using the -O flag)\n\n2. One Bayesian network per biological function (placed in the directory\n\n   specified by the -o flag)\n\n3. Predicted probabilities of functional relationships within each biological\n\n   function of interest (placed in the directory specified by the -p flag)\n\n4. Trust scores for each input data set and function indicating how\n\n   predictive a data set is within a function (specified by the -t flag)";

const char *gengetopt_args_info_usage = "Usage: MEFIT [OPTIONS]... [FILES]...";

const char *gengetopt_args_info_description = "";

const char *gengetopt_args_info_help[] = {
  "  -h, --help                   Print help and exit",
  "  -V, --version                Print version and exit",
  "\nInputs:",
  "  -r, --related=directory      Directory containing lists of known related \n                                 genes",
  "  -u, --unrelated=filename     List of known unrelated gene pairs",
  "  -d, --distance=STRING        Similarity measure  (possible \n                                 values=\"pearson\", \"euclidean\", \n                                 \"kendalls\", \"kolm-smir\", \"spearman\", \n                                 \"pearnorm\" default=`pearnorm')",
  "  -b, --bins=filename          Tab separated QUANT bin cutoffs",
  "\nOutputs:",
  "  -o, --output=directory       Directory to contain learned per-function \n                                 Bayesian networks",
  "  -O, --global=filename        Global learned Bayesian network",
  "  -p, --predictions=directory  Directory to contain predicted probabilities of \n                                 functional relationship",
  "  -t, --trusts=filename        Trust scores learned per data set and function",
  "\nLearning/Evaluation/Features:",
  "  -g, --genes=filename         Subset of genes to include in evaluation",
  "  -G, --genex=filename         Subset of genes to exclude from evaluation",
  "  -z, --zero                   Zero missing values  (default=off)",
  "  -c, --cutoff=DOUBLE          Include only confidences above cutoff  \n                                 (default=`0')",
  "  -s, --skip=INT               Additional columns to skip in input PCLs  \n                                 (default=`2')",
  "  -x, --xdsl                   Output XDSL files in place of DSLs  (default=on)",
  "  -a, --dab                    Output DAB files in place of DATs  (default=on)",
  "  -R, --random=INT             Seed random generator  (default=`0')",
  "  -v, --verbosity=INT          Message verbosity  (default=`5')",
    0
};

typedef enum {ARG_NO
  , ARG_FLAG
  , ARG_STRING
  , ARG_INT
  , ARG_DOUBLE
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

char *cmdline_parser_distance_values[] = {"pearson", "euclidean", "kendalls", "kolm-smir", "spearman", "pearnorm", 0} ;	/* Possible values for distance.  */

static char *
gengetopt_strdup (const char *s);

static
void clear_given (struct gengetopt_args_info *args_info)
{
  args_info->help_given = 0 ;
  args_info->version_given = 0 ;
  args_info->related_given = 0 ;
  args_info->unrelated_given = 0 ;
  args_info->distance_given = 0 ;
  args_info->bins_given = 0 ;
  args_info->output_given = 0 ;
  args_info->global_given = 0 ;
  args_info->predictions_given = 0 ;
  args_info->trusts_given = 0 ;
  args_info->genes_given = 0 ;
  args_info->genex_given = 0 ;
  args_info->zero_given = 0 ;
  args_info->cutoff_given = 0 ;
  args_info->skip_given = 0 ;
  args_info->xdsl_given = 0 ;
  args_info->dab_given = 0 ;
  args_info->random_given = 0 ;
  args_info->verbosity_given = 0 ;
}

static
void clear_args (struct gengetopt_args_info *args_info)
{
  args_info->related_arg = NULL;
  args_info->related_orig = NULL;
  args_info->unrelated_arg = NULL;
  args_info->unrelated_orig = NULL;
  args_info->distance_arg = gengetopt_strdup ("pearnorm");
  args_info->distance_orig = NULL;
  args_info->bins_arg = NULL;
  args_info->bins_orig = NULL;
  args_info->output_arg = NULL;
  args_info->output_orig = NULL;
  args_info->global_arg = NULL;
  args_info->global_orig = NULL;
  args_info->predictions_arg = NULL;
  args_info->predictions_orig = NULL;
  args_info->trusts_arg = NULL;
  args_info->trusts_orig = NULL;
  args_info->genes_arg = NULL;
  args_info->genes_orig = NULL;
  args_info->genex_arg = NULL;
  args_info->genex_orig = NULL;
  args_info->zero_flag = 0;
  args_info->cutoff_arg = 0;
  args_info->cutoff_orig = NULL;
  args_info->skip_arg = 2;
  args_info->skip_orig = NULL;
  args_info->xdsl_flag = 1;
  args_info->dab_flag = 1;
  args_info->random_arg = 0;
  args_info->random_orig = NULL;
  args_info->verbosity_arg = 5;
  args_info->verbosity_orig = NULL;
  
}

static
void init_args_info(struct gengetopt_args_info *args_info)
{


  args_info->help_help = gengetopt_args_info_help[0] ;
  args_info->version_help = gengetopt_args_info_help[1] ;
  args_info->related_help = gengetopt_args_info_help[3] ;
  args_info->unrelated_help = gengetopt_args_info_help[4] ;
  args_info->distance_help = gengetopt_args_info_help[5] ;
  args_info->bins_help = gengetopt_args_info_help[6] ;
  args_info->output_help = gengetopt_args_info_help[8] ;
  args_info->global_help = gengetopt_args_info_help[9] ;
  args_info->predictions_help = gengetopt_args_info_help[10] ;
  args_info->trusts_help = gengetopt_args_info_help[11] ;
  args_info->genes_help = gengetopt_args_info_help[13] ;
  args_info->genex_help = gengetopt_args_info_help[14] ;
  args_info->zero_help = gengetopt_args_info_help[15] ;
  args_info->cutoff_help = gengetopt_args_info_help[16] ;
  args_info->skip_help = gengetopt_args_info_help[17] ;
  args_info->xdsl_help = gengetopt_args_info_help[18] ;
  args_info->dab_help = gengetopt_args_info_help[19] ;
  args_info->random_help = gengetopt_args_info_help[20] ;
  args_info->verbosity_help = gengetopt_args_info_help[21] ;
  
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
  free_string_field (&(args_info->related_arg));
  free_string_field (&(args_info->related_orig));
  free_string_field (&(args_info->unrelated_arg));
  free_string_field (&(args_info->unrelated_orig));
  free_string_field (&(args_info->distance_arg));
  free_string_field (&(args_info->distance_orig));
  free_string_field (&(args_info->bins_arg));
  free_string_field (&(args_info->bins_orig));
  free_string_field (&(args_info->output_arg));
  free_string_field (&(args_info->output_orig));
  free_string_field (&(args_info->global_arg));
  free_string_field (&(args_info->global_orig));
  free_string_field (&(args_info->predictions_arg));
  free_string_field (&(args_info->predictions_orig));
  free_string_field (&(args_info->trusts_arg));
  free_string_field (&(args_info->trusts_orig));
  free_string_field (&(args_info->genes_arg));
  free_string_field (&(args_info->genes_orig));
  free_string_field (&(args_info->genex_arg));
  free_string_field (&(args_info->genex_orig));
  free_string_field (&(args_info->cutoff_orig));
  free_string_field (&(args_info->skip_orig));
  free_string_field (&(args_info->random_orig));
  free_string_field (&(args_info->verbosity_orig));
  
  
  for (i = 0; i < args_info->inputs_num; ++i)
    free (args_info->inputs [i]);

  if (args_info->inputs_num)
    free (args_info->inputs);

  clear_given (args_info);
}

/**
 * @param val the value to check
 * @param values the possible values
 * @return the index of the matched value:
 * -1 if no value matched,
 * -2 if more than one value has matched
 */
static int
check_possible_values(const char *val, char *values[])
{
  int i, found, last;
  size_t len;

  if (!val)   /* otherwise strlen() crashes below */
    return -1; /* -1 means no argument for the option */

  found = last = 0;

  for (i = 0, len = strlen(val); values[i]; ++i)
    {
      if (strncmp(val, values[i], len) == 0)
        {
          ++found;
          last = i;
          if (strlen(values[i]) == len)
            return i; /* exact macth no need to check more */
        }
    }

  if (found == 1) /* one match: OK */
    return last;

  return (found ? -2 : -1); /* return many values or none matched */
}


static void
write_into_file(FILE *outfile, const char *opt, const char *arg, char *values[])
{
  int found = -1;
  if (arg) {
    if (values) {
      found = check_possible_values(arg, values);      
    }
    if (found >= 0)
      fprintf(outfile, "%s=\"%s\" # %s\n", opt, arg, values[found]);
    else
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
  if (args_info->related_given)
    write_into_file(outfile, "related", args_info->related_orig, 0);
  if (args_info->unrelated_given)
    write_into_file(outfile, "unrelated", args_info->unrelated_orig, 0);
  if (args_info->distance_given)
    write_into_file(outfile, "distance", args_info->distance_orig, cmdline_parser_distance_values);
  if (args_info->bins_given)
    write_into_file(outfile, "bins", args_info->bins_orig, 0);
  if (args_info->output_given)
    write_into_file(outfile, "output", args_info->output_orig, 0);
  if (args_info->global_given)
    write_into_file(outfile, "global", args_info->global_orig, 0);
  if (args_info->predictions_given)
    write_into_file(outfile, "predictions", args_info->predictions_orig, 0);
  if (args_info->trusts_given)
    write_into_file(outfile, "trusts", args_info->trusts_orig, 0);
  if (args_info->genes_given)
    write_into_file(outfile, "genes", args_info->genes_orig, 0);
  if (args_info->genex_given)
    write_into_file(outfile, "genex", args_info->genex_orig, 0);
  if (args_info->zero_given)
    write_into_file(outfile, "zero", 0, 0 );
  if (args_info->cutoff_given)
    write_into_file(outfile, "cutoff", args_info->cutoff_orig, 0);
  if (args_info->skip_given)
    write_into_file(outfile, "skip", args_info->skip_orig, 0);
  if (args_info->xdsl_given)
    write_into_file(outfile, "xdsl", 0, 0 );
  if (args_info->dab_given)
    write_into_file(outfile, "dab", 0, 0 );
  if (args_info->random_given)
    write_into_file(outfile, "random", args_info->random_orig, 0);
  if (args_info->verbosity_given)
    write_into_file(outfile, "verbosity", args_info->verbosity_orig, 0);
  

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
  if (! args_info->related_given)
    {
      fprintf (stderr, "%s: '--related' ('-r') option required%s\n", prog_name, (additional_error ? additional_error : ""));
      error = 1;
    }
  
  if (! args_info->unrelated_given)
    {
      fprintf (stderr, "%s: '--unrelated' ('-u') option required%s\n", prog_name, (additional_error ? additional_error : ""));
      error = 1;
    }
  
  if (! args_info->output_given)
    {
      fprintf (stderr, "%s: '--output' ('-o') option required%s\n", prog_name, (additional_error ? additional_error : ""));
      error = 1;
    }
  
  if (! args_info->global_given)
    {
      fprintf (stderr, "%s: '--global' ('-O') option required%s\n", prog_name, (additional_error ? additional_error : ""));
      error = 1;
    }
  
  if (! args_info->predictions_given)
    {
      fprintf (stderr, "%s: '--predictions' ('-p') option required%s\n", prog_name, (additional_error ? additional_error : ""));
      error = 1;
    }
  
  if (! args_info->trusts_given)
    {
      fprintf (stderr, "%s: '--trusts' ('-t') option required%s\n", prog_name, (additional_error ? additional_error : ""));
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

  if (possible_values && (found = check_possible_values((value ? value : default_value), possible_values)) < 0)
    {
      if (short_opt != '-')
        fprintf (stderr, "%s: %s argument, \"%s\", for option `--%s' (`-%c')%s\n", 
          package_name, (found == -2) ? "ambiguous" : "invalid", value, long_opt, short_opt,
          (additional_error ? additional_error : ""));
      else
        fprintf (stderr, "%s: %s argument, \"%s\", for option `--%s'%s\n", 
          package_name, (found == -2) ? "ambiguous" : "invalid", value, long_opt,
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
  case ARG_DOUBLE:
    if (val) *((double *)field) = strtod (val, &stop_char);
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
  case ARG_DOUBLE:
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
        { "related",	1, NULL, 'r' },
        { "unrelated",	1, NULL, 'u' },
        { "distance",	1, NULL, 'd' },
        { "bins",	1, NULL, 'b' },
        { "output",	1, NULL, 'o' },
        { "global",	1, NULL, 'O' },
        { "predictions",	1, NULL, 'p' },
        { "trusts",	1, NULL, 't' },
        { "genes",	1, NULL, 'g' },
        { "genex",	1, NULL, 'G' },
        { "zero",	0, NULL, 'z' },
        { "cutoff",	1, NULL, 'c' },
        { "skip",	1, NULL, 's' },
        { "xdsl",	0, NULL, 'x' },
        { "dab",	0, NULL, 'a' },
        { "random",	1, NULL, 'R' },
        { "verbosity",	1, NULL, 'v' },
        { NULL,	0, NULL, 0 }
      };

      c = getopt_long (argc, argv, "hVr:u:d:b:o:O:p:t:g:G:zc:s:xaR:v:", long_options, &option_index);

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
        case 'r':	/* Directory containing lists of known related genes.  */
        
        
          if (update_arg( (void *)&(args_info->related_arg), 
               &(args_info->related_orig), &(args_info->related_given),
              &(local_args_info.related_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "related", 'r',
              additional_error))
            goto failure;
        
          break;
        case 'u':	/* List of known unrelated gene pairs.  */
        
        
          if (update_arg( (void *)&(args_info->unrelated_arg), 
               &(args_info->unrelated_orig), &(args_info->unrelated_given),
              &(local_args_info.unrelated_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "unrelated", 'u',
              additional_error))
            goto failure;
        
          break;
        case 'd':	/* Similarity measure.  */
        
        
          if (update_arg( (void *)&(args_info->distance_arg), 
               &(args_info->distance_orig), &(args_info->distance_given),
              &(local_args_info.distance_given), optarg, cmdline_parser_distance_values, "pearnorm", ARG_STRING,
              check_ambiguity, override, 0, 0,
              "distance", 'd',
              additional_error))
            goto failure;
        
          break;
        case 'b':	/* Tab separated QUANT bin cutoffs.  */
        
        
          if (update_arg( (void *)&(args_info->bins_arg), 
               &(args_info->bins_orig), &(args_info->bins_given),
              &(local_args_info.bins_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "bins", 'b',
              additional_error))
            goto failure;
        
          break;
        case 'o':	/* Directory to contain learned per-function Bayesian networks.  */
        
        
          if (update_arg( (void *)&(args_info->output_arg), 
               &(args_info->output_orig), &(args_info->output_given),
              &(local_args_info.output_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "output", 'o',
              additional_error))
            goto failure;
        
          break;
        case 'O':	/* Global learned Bayesian network.  */
        
        
          if (update_arg( (void *)&(args_info->global_arg), 
               &(args_info->global_orig), &(args_info->global_given),
              &(local_args_info.global_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "global", 'O',
              additional_error))
            goto failure;
        
          break;
        case 'p':	/* Directory to contain predicted probabilities of functional relationship.  */
        
        
          if (update_arg( (void *)&(args_info->predictions_arg), 
               &(args_info->predictions_orig), &(args_info->predictions_given),
              &(local_args_info.predictions_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "predictions", 'p',
              additional_error))
            goto failure;
        
          break;
        case 't':	/* Trust scores learned per data set and function.  */
        
        
          if (update_arg( (void *)&(args_info->trusts_arg), 
               &(args_info->trusts_orig), &(args_info->trusts_given),
              &(local_args_info.trusts_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "trusts", 't',
              additional_error))
            goto failure;
        
          break;
        case 'g':	/* Subset of genes to include in evaluation.  */
        
        
          if (update_arg( (void *)&(args_info->genes_arg), 
               &(args_info->genes_orig), &(args_info->genes_given),
              &(local_args_info.genes_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "genes", 'g',
              additional_error))
            goto failure;
        
          break;
        case 'G':	/* Subset of genes to exclude from evaluation.  */
        
        
          if (update_arg( (void *)&(args_info->genex_arg), 
               &(args_info->genex_orig), &(args_info->genex_given),
              &(local_args_info.genex_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "genex", 'G',
              additional_error))
            goto failure;
        
          break;
        case 'z':	/* Zero missing values.  */
        
        
          if (update_arg((void *)&(args_info->zero_flag), 0, &(args_info->zero_given),
              &(local_args_info.zero_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "zero", 'z',
              additional_error))
            goto failure;
        
          break;
        case 'c':	/* Include only confidences above cutoff.  */
        
        
          if (update_arg( (void *)&(args_info->cutoff_arg), 
               &(args_info->cutoff_orig), &(args_info->cutoff_given),
              &(local_args_info.cutoff_given), optarg, 0, "0", ARG_DOUBLE,
              check_ambiguity, override, 0, 0,
              "cutoff", 'c',
              additional_error))
            goto failure;
        
          break;
        case 's':	/* Additional columns to skip in input PCLs.  */
        
        
          if (update_arg( (void *)&(args_info->skip_arg), 
               &(args_info->skip_orig), &(args_info->skip_given),
              &(local_args_info.skip_given), optarg, 0, "2", ARG_INT,
              check_ambiguity, override, 0, 0,
              "skip", 's',
              additional_error))
            goto failure;
        
          break;
        case 'x':	/* Output XDSL files in place of DSLs.  */
        
        
          if (update_arg((void *)&(args_info->xdsl_flag), 0, &(args_info->xdsl_given),
              &(local_args_info.xdsl_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "xdsl", 'x',
              additional_error))
            goto failure;
        
          break;
        case 'a':	/* Output DAB files in place of DATs.  */
        
        
          if (update_arg((void *)&(args_info->dab_flag), 0, &(args_info->dab_given),
              &(local_args_info.dab_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "dab", 'a',
              additional_error))
            goto failure;
        
          break;
        case 'R':	/* Seed random generator.  */
        
        
          if (update_arg( (void *)&(args_info->random_arg), 
               &(args_info->random_orig), &(args_info->random_given),
              &(local_args_info.random_given), optarg, 0, "0", ARG_INT,
              check_ambiguity, override, 0, 0,
              "random", 'R',
              additional_error))
            goto failure;
        
          break;
        case 'v':	/* Message verbosity.  */
        
        
          if (update_arg( (void *)&(args_info->verbosity_arg), 
               &(args_info->verbosity_orig), &(args_info->verbosity_given),
              &(local_args_info.verbosity_given), optarg, 0, "5", ARG_INT,
              check_ambiguity, override, 0, 0,
              "verbosity", 'v',
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
