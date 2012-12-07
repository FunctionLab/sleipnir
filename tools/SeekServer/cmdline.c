/*
  File autogenerated by gengetopt version 2.22.5
  generated with the following command:
  /usr/bin/gengetopt -iSeekServer.ggo --default-optional -u -N -e 

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

const char *gengetopt_args_info_purpose = "Performs cross-platform microarray query-guided search in server mode";

const char *gengetopt_args_info_usage = "Usage: SeekServer [OPTIONS]... [FILES]...";

const char *gengetopt_args_info_description = "";

const char *gengetopt_args_info_help[] = {
  "  -h, --help                    Print help and exit",
  "  -V, --version                 Print version and exit",
  "\nMain:",
  "  -t, --port=STRING             Port  (default=`9000')",
  "  -x, --dset=filename           Input a set of datasets",
  "  -i, --input=filename          Input gene mapping",
  "  -d, --dir_in=directory        Database directory",
  "  -p, --dir_prep_in=directory   Prep directory (containing .gavg, .gpres files)",
  "  -P, --dir_platform=directory  Platform directory (containing .gplatavg, \n                                  .gplatstdev, .gplatorder files)",
  "  -u, --dir_sinfo=directory     Sinfo Directory (containing .sinfo files)  \n                                  (default=`NA')",
  "  -U, --dir_gvar=directory      Gene variance directory (containing .gexpvar \n                                  files)  (default=`NA')",
  "  -Q, --quant=filename          quant file (assuming all datasets use the same \n                                  quantization)",
  "  -n, --num_db=INT              Number of databaselets in database  \n                                  (default=`1000')",
  "\nOptional - Parameter tweaking:",
  "  -T, --correlation             Use Pearson correlation values, instead of \n                                  z-score. -m, -M, -r do not apply  \n                                  (default=off)",
  "  -m, --norm_subavg             Per dataset, normalize z-scores by subtracting \n                                  average of result gene  (default=off)",
  "  -M, --norm_platsubavg         Per platform, normalize z-scores by subtracting \n                                  average of query gene across platform  \n                                  (default=off)",
  "  -r, --norm_platstdev          Per platform, normalize z-scores by dividing \n                                  stdev of query gene across platform  \n                                  (default=off)",
  "  -c, --score_cutoff=FLOAT      Cutoff on the gene-gene score before adding, \n                                  default: no cutoff  (default=`-9999')",
  "  -C, --per_q_required=FLOAT    Fraction (max 1.0) of query required to \n                                  correlate with a gene, in order to count the \n                                  gene's query score. A gene may not correlate \n                                  with a query gene if it is absent, or its \n                                  correlation with query does not pass cut-off \n                                  (specified by --score_cutoff). Use this with \n                                  caution. Be careful if using with \n                                  --score_cutoff.  (default=`0.0')",
  "  -e, --square_z                If using z-score, square-transform z-scores. \n                                  Usually used in conjunction with \n                                  --score-cutoff  (default=off)",
  "\nMISC:",
  "  -N, --is_nibble               Whether the input DB is nibble type  \n                                  (default=off)",
  "  -b, --buffer=INT              Number of Databaselets to store in memory  \n                                  (default=`20')",
  "  -O, --output_text             Output results (gene list and dataset weights) \n                                  as text  (default=off)",
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
  args_info->port_given = 0 ;
  args_info->dset_given = 0 ;
  args_info->input_given = 0 ;
  args_info->dir_in_given = 0 ;
  args_info->dir_prep_in_given = 0 ;
  args_info->dir_platform_given = 0 ;
  args_info->dir_sinfo_given = 0 ;
  args_info->dir_gvar_given = 0 ;
  args_info->quant_given = 0 ;
  args_info->num_db_given = 0 ;
  args_info->correlation_given = 0 ;
  args_info->norm_subavg_given = 0 ;
  args_info->norm_platsubavg_given = 0 ;
  args_info->norm_platstdev_given = 0 ;
  args_info->score_cutoff_given = 0 ;
  args_info->per_q_required_given = 0 ;
  args_info->square_z_given = 0 ;
  args_info->is_nibble_given = 0 ;
  args_info->buffer_given = 0 ;
  args_info->output_text_given = 0 ;
}

static
void clear_args (struct gengetopt_args_info *args_info)
{
  FIX_UNUSED (args_info);
  args_info->port_arg = gengetopt_strdup ("9000");
  args_info->port_orig = NULL;
  args_info->dset_arg = NULL;
  args_info->dset_orig = NULL;
  args_info->input_arg = NULL;
  args_info->input_orig = NULL;
  args_info->dir_in_arg = NULL;
  args_info->dir_in_orig = NULL;
  args_info->dir_prep_in_arg = NULL;
  args_info->dir_prep_in_orig = NULL;
  args_info->dir_platform_arg = NULL;
  args_info->dir_platform_orig = NULL;
  args_info->dir_sinfo_arg = gengetopt_strdup ("NA");
  args_info->dir_sinfo_orig = NULL;
  args_info->dir_gvar_arg = gengetopt_strdup ("NA");
  args_info->dir_gvar_orig = NULL;
  args_info->quant_arg = NULL;
  args_info->quant_orig = NULL;
  args_info->num_db_arg = 1000;
  args_info->num_db_orig = NULL;
  args_info->correlation_flag = 0;
  args_info->norm_subavg_flag = 0;
  args_info->norm_platsubavg_flag = 0;
  args_info->norm_platstdev_flag = 0;
  args_info->score_cutoff_arg = -9999;
  args_info->score_cutoff_orig = NULL;
  args_info->per_q_required_arg = 0.0;
  args_info->per_q_required_orig = NULL;
  args_info->square_z_flag = 0;
  args_info->is_nibble_flag = 0;
  args_info->buffer_arg = 20;
  args_info->buffer_orig = NULL;
  args_info->output_text_flag = 0;
  
}

static
void init_args_info(struct gengetopt_args_info *args_info)
{


  args_info->help_help = gengetopt_args_info_help[0] ;
  args_info->version_help = gengetopt_args_info_help[1] ;
  args_info->port_help = gengetopt_args_info_help[3] ;
  args_info->dset_help = gengetopt_args_info_help[4] ;
  args_info->input_help = gengetopt_args_info_help[5] ;
  args_info->dir_in_help = gengetopt_args_info_help[6] ;
  args_info->dir_prep_in_help = gengetopt_args_info_help[7] ;
  args_info->dir_platform_help = gengetopt_args_info_help[8] ;
  args_info->dir_sinfo_help = gengetopt_args_info_help[9] ;
  args_info->dir_gvar_help = gengetopt_args_info_help[10] ;
  args_info->quant_help = gengetopt_args_info_help[11] ;
  args_info->num_db_help = gengetopt_args_info_help[12] ;
  args_info->correlation_help = gengetopt_args_info_help[14] ;
  args_info->norm_subavg_help = gengetopt_args_info_help[15] ;
  args_info->norm_platsubavg_help = gengetopt_args_info_help[16] ;
  args_info->norm_platstdev_help = gengetopt_args_info_help[17] ;
  args_info->score_cutoff_help = gengetopt_args_info_help[18] ;
  args_info->per_q_required_help = gengetopt_args_info_help[19] ;
  args_info->square_z_help = gengetopt_args_info_help[20] ;
  args_info->is_nibble_help = gengetopt_args_info_help[22] ;
  args_info->buffer_help = gengetopt_args_info_help[23] ;
  args_info->output_text_help = gengetopt_args_info_help[24] ;
  
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
  free_string_field (&(args_info->port_arg));
  free_string_field (&(args_info->port_orig));
  free_string_field (&(args_info->dset_arg));
  free_string_field (&(args_info->dset_orig));
  free_string_field (&(args_info->input_arg));
  free_string_field (&(args_info->input_orig));
  free_string_field (&(args_info->dir_in_arg));
  free_string_field (&(args_info->dir_in_orig));
  free_string_field (&(args_info->dir_prep_in_arg));
  free_string_field (&(args_info->dir_prep_in_orig));
  free_string_field (&(args_info->dir_platform_arg));
  free_string_field (&(args_info->dir_platform_orig));
  free_string_field (&(args_info->dir_sinfo_arg));
  free_string_field (&(args_info->dir_sinfo_orig));
  free_string_field (&(args_info->dir_gvar_arg));
  free_string_field (&(args_info->dir_gvar_orig));
  free_string_field (&(args_info->quant_arg));
  free_string_field (&(args_info->quant_orig));
  free_string_field (&(args_info->num_db_orig));
  free_string_field (&(args_info->score_cutoff_orig));
  free_string_field (&(args_info->per_q_required_orig));
  free_string_field (&(args_info->buffer_orig));
  
  
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
  if (args_info->port_given)
    write_into_file(outfile, "port", args_info->port_orig, 0);
  if (args_info->dset_given)
    write_into_file(outfile, "dset", args_info->dset_orig, 0);
  if (args_info->input_given)
    write_into_file(outfile, "input", args_info->input_orig, 0);
  if (args_info->dir_in_given)
    write_into_file(outfile, "dir_in", args_info->dir_in_orig, 0);
  if (args_info->dir_prep_in_given)
    write_into_file(outfile, "dir_prep_in", args_info->dir_prep_in_orig, 0);
  if (args_info->dir_platform_given)
    write_into_file(outfile, "dir_platform", args_info->dir_platform_orig, 0);
  if (args_info->dir_sinfo_given)
    write_into_file(outfile, "dir_sinfo", args_info->dir_sinfo_orig, 0);
  if (args_info->dir_gvar_given)
    write_into_file(outfile, "dir_gvar", args_info->dir_gvar_orig, 0);
  if (args_info->quant_given)
    write_into_file(outfile, "quant", args_info->quant_orig, 0);
  if (args_info->num_db_given)
    write_into_file(outfile, "num_db", args_info->num_db_orig, 0);
  if (args_info->correlation_given)
    write_into_file(outfile, "correlation", 0, 0 );
  if (args_info->norm_subavg_given)
    write_into_file(outfile, "norm_subavg", 0, 0 );
  if (args_info->norm_platsubavg_given)
    write_into_file(outfile, "norm_platsubavg", 0, 0 );
  if (args_info->norm_platstdev_given)
    write_into_file(outfile, "norm_platstdev", 0, 0 );
  if (args_info->score_cutoff_given)
    write_into_file(outfile, "score_cutoff", args_info->score_cutoff_orig, 0);
  if (args_info->per_q_required_given)
    write_into_file(outfile, "per_q_required", args_info->per_q_required_orig, 0);
  if (args_info->square_z_given)
    write_into_file(outfile, "square_z", 0, 0 );
  if (args_info->is_nibble_given)
    write_into_file(outfile, "is_nibble", 0, 0 );
  if (args_info->buffer_given)
    write_into_file(outfile, "buffer", args_info->buffer_orig, 0);
  if (args_info->output_text_given)
    write_into_file(outfile, "output_text", 0, 0 );
  

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
  if (! args_info->dset_given)
    {
      fprintf (stderr, "%s: '--dset' ('-x') option required%s\n", prog_name, (additional_error ? additional_error : ""));
      error = 1;
    }
  
  if (! args_info->input_given)
    {
      fprintf (stderr, "%s: '--input' ('-i') option required%s\n", prog_name, (additional_error ? additional_error : ""));
      error = 1;
    }
  
  if (! args_info->dir_in_given)
    {
      fprintf (stderr, "%s: '--dir_in' ('-d') option required%s\n", prog_name, (additional_error ? additional_error : ""));
      error = 1;
    }
  
  if (! args_info->dir_prep_in_given)
    {
      fprintf (stderr, "%s: '--dir_prep_in' ('-p') option required%s\n", prog_name, (additional_error ? additional_error : ""));
      error = 1;
    }
  
  if (! args_info->dir_platform_given)
    {
      fprintf (stderr, "%s: '--dir_platform' ('-P') option required%s\n", prog_name, (additional_error ? additional_error : ""));
      error = 1;
    }
  
  if (! args_info->quant_given)
    {
      fprintf (stderr, "%s: '--quant' ('-Q') option required%s\n", prog_name, (additional_error ? additional_error : ""));
      error = 1;
    }
  
  if (! args_info->num_db_given)
    {
      fprintf (stderr, "%s: '--num_db' ('-n') option required%s\n", prog_name, (additional_error ? additional_error : ""));
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
        { "port",	1, NULL, 't' },
        { "dset",	1, NULL, 'x' },
        { "input",	1, NULL, 'i' },
        { "dir_in",	1, NULL, 'd' },
        { "dir_prep_in",	1, NULL, 'p' },
        { "dir_platform",	1, NULL, 'P' },
        { "dir_sinfo",	1, NULL, 'u' },
        { "dir_gvar",	1, NULL, 'U' },
        { "quant",	1, NULL, 'Q' },
        { "num_db",	1, NULL, 'n' },
        { "correlation",	0, NULL, 'T' },
        { "norm_subavg",	0, NULL, 'm' },
        { "norm_platsubavg",	0, NULL, 'M' },
        { "norm_platstdev",	0, NULL, 'r' },
        { "score_cutoff",	1, NULL, 'c' },
        { "per_q_required",	1, NULL, 'C' },
        { "square_z",	0, NULL, 'e' },
        { "is_nibble",	0, NULL, 'N' },
        { "buffer",	1, NULL, 'b' },
        { "output_text",	0, NULL, 'O' },
        { 0,  0, 0, 0 }
      };

      c = getopt_long (argc, argv, "hVt:x:i:d:p:P:u:U:Q:n:TmMrc:C:eNb:O", long_options, &option_index);

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
        case 't':	/* Port.  */
        
        
          if (update_arg( (void *)&(args_info->port_arg), 
               &(args_info->port_orig), &(args_info->port_given),
              &(local_args_info.port_given), optarg, 0, "9000", ARG_STRING,
              check_ambiguity, override, 0, 0,
              "port", 't',
              additional_error))
            goto failure;
        
          break;
        case 'x':	/* Input a set of datasets.  */
        
        
          if (update_arg( (void *)&(args_info->dset_arg), 
               &(args_info->dset_orig), &(args_info->dset_given),
              &(local_args_info.dset_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "dset", 'x',
              additional_error))
            goto failure;
        
          break;
        case 'i':	/* Input gene mapping.  */
        
        
          if (update_arg( (void *)&(args_info->input_arg), 
               &(args_info->input_orig), &(args_info->input_given),
              &(local_args_info.input_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "input", 'i',
              additional_error))
            goto failure;
        
          break;
        case 'd':	/* Database directory.  */
        
        
          if (update_arg( (void *)&(args_info->dir_in_arg), 
               &(args_info->dir_in_orig), &(args_info->dir_in_given),
              &(local_args_info.dir_in_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "dir_in", 'd',
              additional_error))
            goto failure;
        
          break;
        case 'p':	/* Prep directory (containing .gavg, .gpres files).  */
        
        
          if (update_arg( (void *)&(args_info->dir_prep_in_arg), 
               &(args_info->dir_prep_in_orig), &(args_info->dir_prep_in_given),
              &(local_args_info.dir_prep_in_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "dir_prep_in", 'p',
              additional_error))
            goto failure;
        
          break;
        case 'P':	/* Platform directory (containing .gplatavg, .gplatstdev, .gplatorder files).  */
        
        
          if (update_arg( (void *)&(args_info->dir_platform_arg), 
               &(args_info->dir_platform_orig), &(args_info->dir_platform_given),
              &(local_args_info.dir_platform_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "dir_platform", 'P',
              additional_error))
            goto failure;
        
          break;
        case 'u':	/* Sinfo Directory (containing .sinfo files).  */
        
        
          if (update_arg( (void *)&(args_info->dir_sinfo_arg), 
               &(args_info->dir_sinfo_orig), &(args_info->dir_sinfo_given),
              &(local_args_info.dir_sinfo_given), optarg, 0, "NA", ARG_STRING,
              check_ambiguity, override, 0, 0,
              "dir_sinfo", 'u',
              additional_error))
            goto failure;
        
          break;
        case 'U':	/* Gene variance directory (containing .gexpvar files).  */
        
        
          if (update_arg( (void *)&(args_info->dir_gvar_arg), 
               &(args_info->dir_gvar_orig), &(args_info->dir_gvar_given),
              &(local_args_info.dir_gvar_given), optarg, 0, "NA", ARG_STRING,
              check_ambiguity, override, 0, 0,
              "dir_gvar", 'U',
              additional_error))
            goto failure;
        
          break;
        case 'Q':	/* quant file (assuming all datasets use the same quantization).  */
        
        
          if (update_arg( (void *)&(args_info->quant_arg), 
               &(args_info->quant_orig), &(args_info->quant_given),
              &(local_args_info.quant_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "quant", 'Q',
              additional_error))
            goto failure;
        
          break;
        case 'n':	/* Number of databaselets in database.  */
        
        
          if (update_arg( (void *)&(args_info->num_db_arg), 
               &(args_info->num_db_orig), &(args_info->num_db_given),
              &(local_args_info.num_db_given), optarg, 0, "1000", ARG_INT,
              check_ambiguity, override, 0, 0,
              "num_db", 'n',
              additional_error))
            goto failure;
        
          break;
        case 'T':	/* Use Pearson correlation values, instead of z-score. -m, -M, -r do not apply.  */
        
        
          if (update_arg((void *)&(args_info->correlation_flag), 0, &(args_info->correlation_given),
              &(local_args_info.correlation_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "correlation", 'T',
              additional_error))
            goto failure;
        
          break;
        case 'm':	/* Per dataset, normalize z-scores by subtracting average of result gene.  */
        
        
          if (update_arg((void *)&(args_info->norm_subavg_flag), 0, &(args_info->norm_subavg_given),
              &(local_args_info.norm_subavg_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "norm_subavg", 'm',
              additional_error))
            goto failure;
        
          break;
        case 'M':	/* Per platform, normalize z-scores by subtracting average of query gene across platform.  */
        
        
          if (update_arg((void *)&(args_info->norm_platsubavg_flag), 0, &(args_info->norm_platsubavg_given),
              &(local_args_info.norm_platsubavg_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "norm_platsubavg", 'M',
              additional_error))
            goto failure;
        
          break;
        case 'r':	/* Per platform, normalize z-scores by dividing stdev of query gene across platform.  */
        
        
          if (update_arg((void *)&(args_info->norm_platstdev_flag), 0, &(args_info->norm_platstdev_given),
              &(local_args_info.norm_platstdev_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "norm_platstdev", 'r',
              additional_error))
            goto failure;
        
          break;
        case 'c':	/* Cutoff on the gene-gene score before adding, default: no cutoff.  */
        
        
          if (update_arg( (void *)&(args_info->score_cutoff_arg), 
               &(args_info->score_cutoff_orig), &(args_info->score_cutoff_given),
              &(local_args_info.score_cutoff_given), optarg, 0, "-9999", ARG_FLOAT,
              check_ambiguity, override, 0, 0,
              "score_cutoff", 'c',
              additional_error))
            goto failure;
        
          break;
        case 'C':	/* Fraction (max 1.0) of query required to correlate with a gene, in order to count the gene's query score. A gene may not correlate with a query gene if it is absent, or its correlation with query does not pass cut-off (specified by --score_cutoff). Use this with caution. Be careful if using with --score_cutoff..  */
        
        
          if (update_arg( (void *)&(args_info->per_q_required_arg), 
               &(args_info->per_q_required_orig), &(args_info->per_q_required_given),
              &(local_args_info.per_q_required_given), optarg, 0, "0.0", ARG_FLOAT,
              check_ambiguity, override, 0, 0,
              "per_q_required", 'C',
              additional_error))
            goto failure;
        
          break;
        case 'e':	/* If using z-score, square-transform z-scores. Usually used in conjunction with --score-cutoff.  */
        
        
          if (update_arg((void *)&(args_info->square_z_flag), 0, &(args_info->square_z_given),
              &(local_args_info.square_z_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "square_z", 'e',
              additional_error))
            goto failure;
        
          break;
        case 'N':	/* Whether the input DB is nibble type.  */
        
        
          if (update_arg((void *)&(args_info->is_nibble_flag), 0, &(args_info->is_nibble_given),
              &(local_args_info.is_nibble_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "is_nibble", 'N',
              additional_error))
            goto failure;
        
          break;
        case 'b':	/* Number of Databaselets to store in memory.  */
        
        
          if (update_arg( (void *)&(args_info->buffer_arg), 
               &(args_info->buffer_orig), &(args_info->buffer_given),
              &(local_args_info.buffer_given), optarg, 0, "20", ARG_INT,
              check_ambiguity, override, 0, 0,
              "buffer", 'b',
              additional_error))
            goto failure;
        
          break;
        case 'O':	/* Output results (gene list and dataset weights) as text.  */
        
        
          if (update_arg((void *)&(args_info->output_text_flag), 0, &(args_info->output_text_given),
              &(local_args_info.output_text_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "output_text", 'O',
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
