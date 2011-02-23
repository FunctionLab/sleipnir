/*
  File autogenerated by gengetopt version 2.22.4
  generated with the following command:
  /usr/bin/gengetopt -iCounter.ggo --default-optional -u -N -e 

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

#ifndef FIX_UNUSED
#define FIX_UNUSED(X) (void) (X) /* avoid warnings for unused params */
#endif

#include <getopt.h>

#include "cmdline.h"

const char *gengetopt_args_info_purpose = "Pre-Bayesian learning tool; counts distributions of values in data";

const char *gengetopt_args_info_usage = "Usage: Counter [OPTIONS]... [FILES]...";

const char *gengetopt_args_info_description = "";

const char *gengetopt_args_info_help[] = {
  "  -h, --help                    Print help and exit",
  "  -V, --version                 Print version and exit",
  "\n Group: Mode",
  "  -w, --answers=filename        Answer file",
  "  -k, --counts=directory        Directory containing count files",
  "  -n, --networks=filename       Bayes nets",
  "\nMain:",
  "  -o, --output=filename or directory\n                                Output count directory, Bayes nets, or \n                                  inferences",
  "  -d, --directory=directory     Data directory  (default=`.')",
  "  -s, --datasets=filename       Dataset ID text file",
  "  -e, --genome=filename         Gene ID text file",
  "  -X, --contexts=filename       Context ID text file",
  "\nLearning/Evaluation:",
  "  -g, --genes=filename          Gene inclusion file",
  "  -G, --genex=filename          Gene exclusion file",
  "  -c, --genet=filename          Term inclusion file",
  "  -C, --genee=filename          Edge inclusion file",
  "\nNetwork Features:",
  "  -b, --default=filename        Count file containing defaults for cases with \n                                  missing data",
  "  -Z, --zeros=filename          Read zeroed node IDs/outputs from the given \n                                  file",
  "  -S, --genewise                Evaluate networks assuming genewise contexts  \n                                  (default=off)",
  "\nBayesian Regularization:",
  "  -p, --pseudocounts=FLOAT      Effective number of pseudocounts to use  \n                                  (default=`-1')",
  "  -a, --alphas=filename         File containing equivalent sample sizes \n                                  (alphas) for each node",
  "  -r, --regularize              Automatically regularize based on similarity  \n                                  (default=off)",
  "  -R, --reggroups=filename      Automatically regularize based on given groups",
  "\nOptional:",
  "  -y, --temporary=directory     Directory for temporary files  (default=`.')",
  "  -l, --smile                   Output SMILE (X)DSL files rather than minimal \n                                  networks  (default=off)",
  "  -x, --xdsl                    Generate XDSL output rather than DSL  \n                                  (default=on)",
  "  -m, --memmap                  Memory map input files  (default=off)",
  "  -M, --memmapout               Memory map output files (only for inference \n                                  mode)  (default=off)",
  "  -t, --threads=INT             Maximum number of threads to spawn  \n                                  (default=`-1')",
  "  -v, --verbosity=INT           Message verbosity  (default=`5')",
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
cmdline_parser_internal (int argc, char **argv, struct gengetopt_args_info *args_info,
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
  args_info->answers_given = 0 ;
  args_info->counts_given = 0 ;
  args_info->networks_given = 0 ;
  args_info->output_given = 0 ;
  args_info->directory_given = 0 ;
  args_info->datasets_given = 0 ;
  args_info->genome_given = 0 ;
  args_info->contexts_given = 0 ;
  args_info->genes_given = 0 ;
  args_info->genex_given = 0 ;
  args_info->genet_given = 0 ;
  args_info->genee_given = 0 ;
  args_info->default_given = 0 ;
  args_info->zeros_given = 0 ;
  args_info->genewise_given = 0 ;
  args_info->pseudocounts_given = 0 ;
  args_info->alphas_given = 0 ;
  args_info->regularize_given = 0 ;
  args_info->reggroups_given = 0 ;
  args_info->temporary_given = 0 ;
  args_info->smile_given = 0 ;
  args_info->xdsl_given = 0 ;
  args_info->memmap_given = 0 ;
  args_info->memmapout_given = 0 ;
  args_info->threads_given = 0 ;
  args_info->verbosity_given = 0 ;
  args_info->Mode_group_counter = 0 ;
}

static
void clear_args (struct gengetopt_args_info *args_info)
{
  FIX_UNUSED (args_info);
  args_info->answers_arg = NULL;
  args_info->answers_orig = NULL;
  args_info->counts_arg = NULL;
  args_info->counts_orig = NULL;
  args_info->networks_arg = NULL;
  args_info->networks_orig = NULL;
  args_info->output_arg = NULL;
  args_info->output_orig = NULL;
  args_info->directory_arg = gengetopt_strdup (".");
  args_info->directory_orig = NULL;
  args_info->datasets_arg = NULL;
  args_info->datasets_orig = NULL;
  args_info->genome_arg = NULL;
  args_info->genome_orig = NULL;
  args_info->contexts_arg = NULL;
  args_info->contexts_orig = NULL;
  args_info->genes_arg = NULL;
  args_info->genes_orig = NULL;
  args_info->genex_arg = NULL;
  args_info->genex_orig = NULL;
  args_info->genet_arg = NULL;
  args_info->genet_orig = NULL;
  args_info->genee_arg = NULL;
  args_info->genee_orig = NULL;
  args_info->default_arg = NULL;
  args_info->default_orig = NULL;
  args_info->zeros_arg = NULL;
  args_info->zeros_orig = NULL;
  args_info->genewise_flag = 0;
  args_info->pseudocounts_arg = -1;
  args_info->pseudocounts_orig = NULL;
  args_info->alphas_arg = NULL;
  args_info->alphas_orig = NULL;
  args_info->regularize_flag = 0;
  args_info->reggroups_arg = NULL;
  args_info->reggroups_orig = NULL;
  args_info->temporary_arg = gengetopt_strdup (".");
  args_info->temporary_orig = NULL;
  args_info->smile_flag = 0;
  args_info->xdsl_flag = 1;
  args_info->memmap_flag = 0;
  args_info->memmapout_flag = 0;
  args_info->threads_arg = -1;
  args_info->threads_orig = NULL;
  args_info->verbosity_arg = 5;
  args_info->verbosity_orig = NULL;
  
}

static
void init_args_info(struct gengetopt_args_info *args_info)
{


  args_info->help_help = gengetopt_args_info_help[0] ;
  args_info->version_help = gengetopt_args_info_help[1] ;
  args_info->answers_help = gengetopt_args_info_help[3] ;
  args_info->counts_help = gengetopt_args_info_help[4] ;
  args_info->networks_help = gengetopt_args_info_help[5] ;
  args_info->output_help = gengetopt_args_info_help[7] ;
  args_info->directory_help = gengetopt_args_info_help[8] ;
  args_info->datasets_help = gengetopt_args_info_help[9] ;
  args_info->genome_help = gengetopt_args_info_help[10] ;
  args_info->contexts_help = gengetopt_args_info_help[11] ;
  args_info->genes_help = gengetopt_args_info_help[13] ;
  args_info->genex_help = gengetopt_args_info_help[14] ;
  args_info->genet_help = gengetopt_args_info_help[15] ;
  args_info->genee_help = gengetopt_args_info_help[16] ;
  args_info->default_help = gengetopt_args_info_help[18] ;
  args_info->zeros_help = gengetopt_args_info_help[19] ;
  args_info->genewise_help = gengetopt_args_info_help[20] ;
  args_info->pseudocounts_help = gengetopt_args_info_help[22] ;
  args_info->alphas_help = gengetopt_args_info_help[23] ;
  args_info->regularize_help = gengetopt_args_info_help[24] ;
  args_info->reggroups_help = gengetopt_args_info_help[25] ;
  args_info->temporary_help = gengetopt_args_info_help[27] ;
  args_info->smile_help = gengetopt_args_info_help[28] ;
  args_info->xdsl_help = gengetopt_args_info_help[29] ;
  args_info->memmap_help = gengetopt_args_info_help[30] ;
  args_info->memmapout_help = gengetopt_args_info_help[31] ;
  args_info->threads_help = gengetopt_args_info_help[32] ;
  args_info->verbosity_help = gengetopt_args_info_help[33] ;
  
}

void
cmdline_parser_print_version (void)
{
  printf ("%s %s\n",
     (strlen(CMDLINE_PARSER_PACKAGE_NAME) ? CMDLINE_PARSER_PACKAGE_NAME : CMDLINE_PARSER_PACKAGE),
     CMDLINE_PARSER_VERSION);
}

static void print_help_common(void) {
  cmdline_parser_print_version ();

  if (strlen(gengetopt_args_info_purpose) > 0)
    printf("\n%s\n", gengetopt_args_info_purpose);

  if (strlen(gengetopt_args_info_usage) > 0)
    printf("\n%s\n", gengetopt_args_info_usage);

  printf("\n");

  if (strlen(gengetopt_args_info_description) > 0)
    printf("%s\n\n", gengetopt_args_info_description);
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

  args_info->inputs = 0;
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
  free_string_field (&(args_info->answers_arg));
  free_string_field (&(args_info->answers_orig));
  free_string_field (&(args_info->counts_arg));
  free_string_field (&(args_info->counts_orig));
  free_string_field (&(args_info->networks_arg));
  free_string_field (&(args_info->networks_orig));
  free_string_field (&(args_info->output_arg));
  free_string_field (&(args_info->output_orig));
  free_string_field (&(args_info->directory_arg));
  free_string_field (&(args_info->directory_orig));
  free_string_field (&(args_info->datasets_arg));
  free_string_field (&(args_info->datasets_orig));
  free_string_field (&(args_info->genome_arg));
  free_string_field (&(args_info->genome_orig));
  free_string_field (&(args_info->contexts_arg));
  free_string_field (&(args_info->contexts_orig));
  free_string_field (&(args_info->genes_arg));
  free_string_field (&(args_info->genes_orig));
  free_string_field (&(args_info->genex_arg));
  free_string_field (&(args_info->genex_orig));
  free_string_field (&(args_info->genet_arg));
  free_string_field (&(args_info->genet_orig));
  free_string_field (&(args_info->genee_arg));
  free_string_field (&(args_info->genee_orig));
  free_string_field (&(args_info->default_arg));
  free_string_field (&(args_info->default_orig));
  free_string_field (&(args_info->zeros_arg));
  free_string_field (&(args_info->zeros_orig));
  free_string_field (&(args_info->pseudocounts_orig));
  free_string_field (&(args_info->alphas_arg));
  free_string_field (&(args_info->alphas_orig));
  free_string_field (&(args_info->reggroups_arg));
  free_string_field (&(args_info->reggroups_orig));
  free_string_field (&(args_info->temporary_arg));
  free_string_field (&(args_info->temporary_orig));
  free_string_field (&(args_info->threads_orig));
  free_string_field (&(args_info->verbosity_orig));
  
  
  for (i = 0; i < args_info->inputs_num; ++i)
    free (args_info->inputs [i]);

  if (args_info->inputs_num)
    free (args_info->inputs);

  clear_given (args_info);
}


static void
write_into_file(FILE *outfile, const char *opt, const char *arg, const char *values[])
{
  FIX_UNUSED (values);
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
  if (args_info->answers_given)
    write_into_file(outfile, "answers", args_info->answers_orig, 0);
  if (args_info->counts_given)
    write_into_file(outfile, "counts", args_info->counts_orig, 0);
  if (args_info->networks_given)
    write_into_file(outfile, "networks", args_info->networks_orig, 0);
  if (args_info->output_given)
    write_into_file(outfile, "output", args_info->output_orig, 0);
  if (args_info->directory_given)
    write_into_file(outfile, "directory", args_info->directory_orig, 0);
  if (args_info->datasets_given)
    write_into_file(outfile, "datasets", args_info->datasets_orig, 0);
  if (args_info->genome_given)
    write_into_file(outfile, "genome", args_info->genome_orig, 0);
  if (args_info->contexts_given)
    write_into_file(outfile, "contexts", args_info->contexts_orig, 0);
  if (args_info->genes_given)
    write_into_file(outfile, "genes", args_info->genes_orig, 0);
  if (args_info->genex_given)
    write_into_file(outfile, "genex", args_info->genex_orig, 0);
  if (args_info->genet_given)
    write_into_file(outfile, "genet", args_info->genet_orig, 0);
  if (args_info->genee_given)
    write_into_file(outfile, "genee", args_info->genee_orig, 0);
  if (args_info->default_given)
    write_into_file(outfile, "default", args_info->default_orig, 0);
  if (args_info->zeros_given)
    write_into_file(outfile, "zeros", args_info->zeros_orig, 0);
  if (args_info->genewise_given)
    write_into_file(outfile, "genewise", 0, 0 );
  if (args_info->pseudocounts_given)
    write_into_file(outfile, "pseudocounts", args_info->pseudocounts_orig, 0);
  if (args_info->alphas_given)
    write_into_file(outfile, "alphas", args_info->alphas_orig, 0);
  if (args_info->regularize_given)
    write_into_file(outfile, "regularize", 0, 0 );
  if (args_info->reggroups_given)
    write_into_file(outfile, "reggroups", args_info->reggroups_orig, 0);
  if (args_info->temporary_given)
    write_into_file(outfile, "temporary", args_info->temporary_orig, 0);
  if (args_info->smile_given)
    write_into_file(outfile, "smile", 0, 0 );
  if (args_info->xdsl_given)
    write_into_file(outfile, "xdsl", 0, 0 );
  if (args_info->memmap_given)
    write_into_file(outfile, "memmap", 0, 0 );
  if (args_info->memmapout_given)
    write_into_file(outfile, "memmapout", 0, 0 );
  if (args_info->threads_given)
    write_into_file(outfile, "threads", args_info->threads_orig, 0);
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
  char *result = 0;
  if (!s)
    return result;

  result = (char*)malloc(strlen(s) + 1);
  if (result == (char*)0)
    return (char*)0;
  strcpy(result, s);
  return result;
}

static void
reset_group_Mode(struct gengetopt_args_info *args_info)
{
  if (! args_info->Mode_group_counter)
    return;
  
  args_info->answers_given = 0 ;
  free_string_field (&(args_info->answers_arg));
  free_string_field (&(args_info->answers_orig));
  args_info->counts_given = 0 ;
  free_string_field (&(args_info->counts_arg));
  free_string_field (&(args_info->counts_orig));
  args_info->networks_given = 0 ;
  free_string_field (&(args_info->networks_arg));
  free_string_field (&(args_info->networks_orig));

  args_info->Mode_group_counter = 0;
}

int
cmdline_parser (int argc, char **argv, struct gengetopt_args_info *args_info)
{
  return cmdline_parser2 (argc, argv, args_info, 0, 1, 1);
}

int
cmdline_parser_ext (int argc, char **argv, struct gengetopt_args_info *args_info,
                   struct cmdline_parser_params *params)
{
  int result;
  result = cmdline_parser_internal (argc, argv, args_info, params, 0);

  return result;
}

int
cmdline_parser2 (int argc, char **argv, struct gengetopt_args_info *args_info, int override, int initialize, int check_required)
{
  int result;
  struct cmdline_parser_params params;
  
  params.override = override;
  params.initialize = initialize;
  params.check_required = check_required;
  params.check_ambiguity = 0;
  params.print_errors = 1;

  result = cmdline_parser_internal (argc, argv, args_info, &params, 0);

  return result;
}

int
cmdline_parser_required (struct gengetopt_args_info *args_info, const char *prog_name)
{
  int result = EXIT_SUCCESS;

  if (cmdline_parser_required2(args_info, prog_name, 0) > 0)
    result = EXIT_FAILURE;

  return result;
}

int
cmdline_parser_required2 (struct gengetopt_args_info *args_info, const char *prog_name, const char *additional_error)
{
  int error = 0;
  FIX_UNUSED (additional_error);

  /* checks for required options */
  if (! args_info->output_given)
    {
      fprintf (stderr, "%s: '--output' ('-o') option required%s\n", prog_name, (additional_error ? additional_error : ""));
      error = 1;
    }
  
  if (args_info->Mode_group_counter == 0)
    {
      fprintf (stderr, "%s: %d options of group Mode were given. One is required%s.\n", prog_name, args_info->Mode_group_counter, (additional_error ? additional_error : ""));
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
               char *value, const char *possible_values[],
               const char *default_value,
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
  FIX_UNUSED (field);

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

  FIX_UNUSED (default_value);
    
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
cmdline_parser_internal (
  int argc, char **argv, struct gengetopt_args_info *args_info,
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
        { "answers",	1, NULL, 'w' },
        { "counts",	1, NULL, 'k' },
        { "networks",	1, NULL, 'n' },
        { "output",	1, NULL, 'o' },
        { "directory",	1, NULL, 'd' },
        { "datasets",	1, NULL, 's' },
        { "genome",	1, NULL, 'e' },
        { "contexts",	1, NULL, 'X' },
        { "genes",	1, NULL, 'g' },
        { "genex",	1, NULL, 'G' },
        { "genet",	1, NULL, 'c' },
        { "genee",	1, NULL, 'C' },
        { "default",	1, NULL, 'b' },
        { "zeros",	1, NULL, 'Z' },
        { "genewise",	0, NULL, 'S' },
        { "pseudocounts",	1, NULL, 'p' },
        { "alphas",	1, NULL, 'a' },
        { "regularize",	0, NULL, 'r' },
        { "reggroups",	1, NULL, 'R' },
        { "temporary",	1, NULL, 'y' },
        { "smile",	0, NULL, 'l' },
        { "xdsl",	0, NULL, 'x' },
        { "memmap",	0, NULL, 'm' },
        { "memmapout",	0, NULL, 'M' },
        { "threads",	1, NULL, 't' },
        { "verbosity",	1, NULL, 'v' },
        { 0,  0, 0, 0 }
      };

      c = getopt_long (argc, argv, "hVw:k:n:o:d:s:e:X:g:G:c:C:b:Z:Sp:a:rR:y:lxmMt:v:", long_options, &option_index);

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
        case 'w':	/* Answer file.  */
        
          if (args_info->Mode_group_counter && override)
            reset_group_Mode (args_info);
          args_info->Mode_group_counter += 1;
        
          if (update_arg( (void *)&(args_info->answers_arg), 
               &(args_info->answers_orig), &(args_info->answers_given),
              &(local_args_info.answers_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "answers", 'w',
              additional_error))
            goto failure;
        
          break;
        case 'k':	/* Directory containing count files.  */
        
          if (args_info->Mode_group_counter && override)
            reset_group_Mode (args_info);
          args_info->Mode_group_counter += 1;
        
          if (update_arg( (void *)&(args_info->counts_arg), 
               &(args_info->counts_orig), &(args_info->counts_given),
              &(local_args_info.counts_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "counts", 'k',
              additional_error))
            goto failure;
        
          break;
        case 'n':	/* Bayes nets.  */
        
          if (args_info->Mode_group_counter && override)
            reset_group_Mode (args_info);
          args_info->Mode_group_counter += 1;
        
          if (update_arg( (void *)&(args_info->networks_arg), 
               &(args_info->networks_orig), &(args_info->networks_given),
              &(local_args_info.networks_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "networks", 'n',
              additional_error))
            goto failure;
        
          break;
        case 'o':	/* Output count directory, Bayes nets, or inferences.  */
        
        
          if (update_arg( (void *)&(args_info->output_arg), 
               &(args_info->output_orig), &(args_info->output_given),
              &(local_args_info.output_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "output", 'o',
              additional_error))
            goto failure;
        
          break;
        case 'd':	/* Data directory.  */
        
        
          if (update_arg( (void *)&(args_info->directory_arg), 
               &(args_info->directory_orig), &(args_info->directory_given),
              &(local_args_info.directory_given), optarg, 0, ".", ARG_STRING,
              check_ambiguity, override, 0, 0,
              "directory", 'd',
              additional_error))
            goto failure;
        
          break;
        case 's':	/* Dataset ID text file.  */
        
        
          if (update_arg( (void *)&(args_info->datasets_arg), 
               &(args_info->datasets_orig), &(args_info->datasets_given),
              &(local_args_info.datasets_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "datasets", 's',
              additional_error))
            goto failure;
        
          break;
        case 'e':	/* Gene ID text file.  */
        
        
          if (update_arg( (void *)&(args_info->genome_arg), 
               &(args_info->genome_orig), &(args_info->genome_given),
              &(local_args_info.genome_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "genome", 'e',
              additional_error))
            goto failure;
        
          break;
        case 'X':	/* Context ID text file.  */
        
        
          if (update_arg( (void *)&(args_info->contexts_arg), 
               &(args_info->contexts_orig), &(args_info->contexts_given),
              &(local_args_info.contexts_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "contexts", 'X',
              additional_error))
            goto failure;
        
          break;
        case 'g':	/* Gene inclusion file.  */
        
        
          if (update_arg( (void *)&(args_info->genes_arg), 
               &(args_info->genes_orig), &(args_info->genes_given),
              &(local_args_info.genes_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "genes", 'g',
              additional_error))
            goto failure;
        
          break;
        case 'G':	/* Gene exclusion file.  */
        
        
          if (update_arg( (void *)&(args_info->genex_arg), 
               &(args_info->genex_orig), &(args_info->genex_given),
              &(local_args_info.genex_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "genex", 'G',
              additional_error))
            goto failure;
        
          break;
        case 'c':	/* Term inclusion file.  */
        
        
          if (update_arg( (void *)&(args_info->genet_arg), 
               &(args_info->genet_orig), &(args_info->genet_given),
              &(local_args_info.genet_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "genet", 'c',
              additional_error))
            goto failure;
        
          break;
        case 'C':	/* Edge inclusion file.  */
        
        
          if (update_arg( (void *)&(args_info->genee_arg), 
               &(args_info->genee_orig), &(args_info->genee_given),
              &(local_args_info.genee_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "genee", 'C',
              additional_error))
            goto failure;
        
          break;
        case 'b':	/* Count file containing defaults for cases with missing data.  */
        
        
          if (update_arg( (void *)&(args_info->default_arg), 
               &(args_info->default_orig), &(args_info->default_given),
              &(local_args_info.default_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "default", 'b',
              additional_error))
            goto failure;
        
          break;
        case 'Z':	/* Read zeroed node IDs/outputs from the given file.  */
        
        
          if (update_arg( (void *)&(args_info->zeros_arg), 
               &(args_info->zeros_orig), &(args_info->zeros_given),
              &(local_args_info.zeros_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "zeros", 'Z',
              additional_error))
            goto failure;
        
          break;
        case 'S':	/* Evaluate networks assuming genewise contexts.  */
        
        
          if (update_arg((void *)&(args_info->genewise_flag), 0, &(args_info->genewise_given),
              &(local_args_info.genewise_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "genewise", 'S',
              additional_error))
            goto failure;
        
          break;
        case 'p':	/* Effective number of pseudocounts to use.  */
        
        
          if (update_arg( (void *)&(args_info->pseudocounts_arg), 
               &(args_info->pseudocounts_orig), &(args_info->pseudocounts_given),
              &(local_args_info.pseudocounts_given), optarg, 0, "-1", ARG_FLOAT,
              check_ambiguity, override, 0, 0,
              "pseudocounts", 'p',
              additional_error))
            goto failure;
        
          break;
        case 'a':	/* File containing equivalent sample sizes (alphas) for each node.  */
        
        
          if (update_arg( (void *)&(args_info->alphas_arg), 
               &(args_info->alphas_orig), &(args_info->alphas_given),
              &(local_args_info.alphas_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "alphas", 'a',
              additional_error))
            goto failure;
        
          break;
        case 'r':	/* Automatically regularize based on similarity.  */
        
        
          if (update_arg((void *)&(args_info->regularize_flag), 0, &(args_info->regularize_given),
              &(local_args_info.regularize_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "regularize", 'r',
              additional_error))
            goto failure;
        
          break;
        case 'R':	/* Automatically regularize based on given groups.  */
        
        
          if (update_arg( (void *)&(args_info->reggroups_arg), 
               &(args_info->reggroups_orig), &(args_info->reggroups_given),
              &(local_args_info.reggroups_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "reggroups", 'R',
              additional_error))
            goto failure;
        
          break;
        case 'y':	/* Directory for temporary files.  */
        
        
          if (update_arg( (void *)&(args_info->temporary_arg), 
               &(args_info->temporary_orig), &(args_info->temporary_given),
              &(local_args_info.temporary_given), optarg, 0, ".", ARG_STRING,
              check_ambiguity, override, 0, 0,
              "temporary", 'y',
              additional_error))
            goto failure;
        
          break;
        case 'l':	/* Output SMILE (X)DSL files rather than minimal networks.  */
        
        
          if (update_arg((void *)&(args_info->smile_flag), 0, &(args_info->smile_given),
              &(local_args_info.smile_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "smile", 'l',
              additional_error))
            goto failure;
        
          break;
        case 'x':	/* Generate XDSL output rather than DSL.  */
        
        
          if (update_arg((void *)&(args_info->xdsl_flag), 0, &(args_info->xdsl_given),
              &(local_args_info.xdsl_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "xdsl", 'x',
              additional_error))
            goto failure;
        
          break;
        case 'm':	/* Memory map input files.  */
        
        
          if (update_arg((void *)&(args_info->memmap_flag), 0, &(args_info->memmap_given),
              &(local_args_info.memmap_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "memmap", 'm',
              additional_error))
            goto failure;
        
          break;
        case 'M':	/* Memory map output files (only for inference mode).  */
        
        
          if (update_arg((void *)&(args_info->memmapout_flag), 0, &(args_info->memmapout_given),
              &(local_args_info.memmapout_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "memmapout", 'M',
              additional_error))
            goto failure;
        
          break;
        case 't':	/* Maximum number of threads to spawn.  */
        
        
          if (update_arg( (void *)&(args_info->threads_arg), 
               &(args_info->threads_orig), &(args_info->threads_given),
              &(local_args_info.threads_given), optarg, 0, "-1", ARG_INT,
              check_ambiguity, override, 0, 0,
              "threads", 't',
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

  if (args_info->Mode_group_counter > 1)
    {
      fprintf (stderr, "%s: %d options of group Mode were given. One is required%s.\n", argv[0], args_info->Mode_group_counter, (additional_error ? additional_error : ""));
      error = 1;
    }
  


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
