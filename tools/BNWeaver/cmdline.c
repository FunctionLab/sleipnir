/*
  File autogenerated by gengetopt version 2.22
  generated with the following command:
  /r01/tergeo/chuttenh/sleipnir/trunk/../extlib/gengetopt-2.22/bin/gengetopt -iBNWeaver.ggo --default-optional -u -N -e 

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

const char *gengetopt_args_info_purpose = "Bayes net construction and training from data";

const char *gengetopt_args_info_usage = "Usage: BNWeaver [OPTIONS]... [FILES]...";

const char *gengetopt_args_info_description = "";

const char *gengetopt_args_info_help[] = {
  "  -h, --help                 Print help and exit",
  "  -V, --version              Print version and exit",
  "\nMain:",
  "  -w, --answers=filename     Answer file",
  "  -o, --output=directory     Output directory  (default=`.')",
  "  -d, --directory=directory  Data directory  (default=`.')",
  "\nLearning/Evaluation:",
  "  -G, --genex=filename       Gene exclusion file",
  "  -n, --negatives=filename   Gene set for negative pairs",
  "  -a, --randomize            Randomize data before training  (default=off)",
  "\nNetwork Features:",
  "  -b, --default=filename     Bayes net containing defaults for cases with \n                               missing data",
  "  -z, --zero                 Zero missing values  (default=off)",
  "  -Z, --zeros=filename       Read zeroed node IDs/outputs from the given file",
  "\nOptional:",
  "  -m, --memmap               Memory map input files  (default=off)",
  "  -t, --threads=INT          Maximum number of threads to spawn  (default=`-1')",
  "  -x, --xdsl                 Generate XDSL output rather than DSL  (default=on)",
  "  -u, --group                Group identical inputs  (default=on)",
  "  -r, --random=INT           Seed random generator  (default=`0')",
  "  -v, --verbosity=INT        Message verbosity  (default=`5')",
    0
};

typedef enum {ARG_NO
  , ARG_FLAG
  , ARG_STRING
  , ARG_INT
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
  args_info->answers_given = 0 ;
  args_info->output_given = 0 ;
  args_info->directory_given = 0 ;
  args_info->genex_given = 0 ;
  args_info->negatives_given = 0 ;
  args_info->randomize_given = 0 ;
  args_info->default_given = 0 ;
  args_info->zero_given = 0 ;
  args_info->zeros_given = 0 ;
  args_info->memmap_given = 0 ;
  args_info->threads_given = 0 ;
  args_info->xdsl_given = 0 ;
  args_info->group_given = 0 ;
  args_info->random_given = 0 ;
  args_info->verbosity_given = 0 ;
}

static
void clear_args (struct gengetopt_args_info *args_info)
{
  args_info->answers_arg = NULL;
  args_info->answers_orig = NULL;
  args_info->output_arg = gengetopt_strdup (".");
  args_info->output_orig = NULL;
  args_info->directory_arg = gengetopt_strdup (".");
  args_info->directory_orig = NULL;
  args_info->genex_arg = NULL;
  args_info->genex_orig = NULL;
  args_info->negatives_arg = NULL;
  args_info->negatives_orig = NULL;
  args_info->randomize_flag = 0;
  args_info->default_arg = NULL;
  args_info->default_orig = NULL;
  args_info->zero_flag = 0;
  args_info->zeros_arg = NULL;
  args_info->zeros_orig = NULL;
  args_info->memmap_flag = 0;
  args_info->threads_arg = -1;
  args_info->threads_orig = NULL;
  args_info->xdsl_flag = 1;
  args_info->group_flag = 1;
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
  args_info->answers_help = gengetopt_args_info_help[3] ;
  args_info->output_help = gengetopt_args_info_help[4] ;
  args_info->directory_help = gengetopt_args_info_help[5] ;
  args_info->genex_help = gengetopt_args_info_help[7] ;
  args_info->negatives_help = gengetopt_args_info_help[8] ;
  args_info->randomize_help = gengetopt_args_info_help[9] ;
  args_info->default_help = gengetopt_args_info_help[11] ;
  args_info->zero_help = gengetopt_args_info_help[12] ;
  args_info->zeros_help = gengetopt_args_info_help[13] ;
  args_info->memmap_help = gengetopt_args_info_help[15] ;
  args_info->threads_help = gengetopt_args_info_help[16] ;
  args_info->xdsl_help = gengetopt_args_info_help[17] ;
  args_info->group_help = gengetopt_args_info_help[18] ;
  args_info->random_help = gengetopt_args_info_help[19] ;
  args_info->verbosity_help = gengetopt_args_info_help[20] ;
  
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
  free_string_field (&(args_info->answers_arg));
  free_string_field (&(args_info->answers_orig));
  free_string_field (&(args_info->output_arg));
  free_string_field (&(args_info->output_orig));
  free_string_field (&(args_info->directory_arg));
  free_string_field (&(args_info->directory_orig));
  free_string_field (&(args_info->genex_arg));
  free_string_field (&(args_info->genex_orig));
  free_string_field (&(args_info->negatives_arg));
  free_string_field (&(args_info->negatives_orig));
  free_string_field (&(args_info->default_arg));
  free_string_field (&(args_info->default_orig));
  free_string_field (&(args_info->zeros_arg));
  free_string_field (&(args_info->zeros_orig));
  free_string_field (&(args_info->threads_orig));
  free_string_field (&(args_info->random_orig));
  free_string_field (&(args_info->verbosity_orig));
  
  
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
  if (args_info->answers_given)
    write_into_file(outfile, "answers", args_info->answers_orig, 0);
  if (args_info->output_given)
    write_into_file(outfile, "output", args_info->output_orig, 0);
  if (args_info->directory_given)
    write_into_file(outfile, "directory", args_info->directory_orig, 0);
  if (args_info->genex_given)
    write_into_file(outfile, "genex", args_info->genex_orig, 0);
  if (args_info->negatives_given)
    write_into_file(outfile, "negatives", args_info->negatives_orig, 0);
  if (args_info->randomize_given)
    write_into_file(outfile, "randomize", 0, 0 );
  if (args_info->default_given)
    write_into_file(outfile, "default", args_info->default_orig, 0);
  if (args_info->zero_given)
    write_into_file(outfile, "zero", 0, 0 );
  if (args_info->zeros_given)
    write_into_file(outfile, "zeros", args_info->zeros_orig, 0);
  if (args_info->memmap_given)
    write_into_file(outfile, "memmap", 0, 0 );
  if (args_info->threads_given)
    write_into_file(outfile, "threads", args_info->threads_orig, 0);
  if (args_info->xdsl_given)
    write_into_file(outfile, "xdsl", 0, 0 );
  if (args_info->group_given)
    write_into_file(outfile, "group", 0, 0 );
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
  if (! args_info->answers_given)
    {
      fprintf (stderr, "%s: '--answers' ('-w') option required%s\n", prog_name, (additional_error ? additional_error : ""));
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
        { "answers",	1, NULL, 'w' },
        { "output",	1, NULL, 'o' },
        { "directory",	1, NULL, 'd' },
        { "genex",	1, NULL, 'G' },
        { "negatives",	1, NULL, 'n' },
        { "randomize",	0, NULL, 'a' },
        { "default",	1, NULL, 'b' },
        { "zero",	0, NULL, 'z' },
        { "zeros",	1, NULL, 'Z' },
        { "memmap",	0, NULL, 'm' },
        { "threads",	1, NULL, 't' },
        { "xdsl",	0, NULL, 'x' },
        { "group",	0, NULL, 'u' },
        { "random",	1, NULL, 'r' },
        { "verbosity",	1, NULL, 'v' },
        { NULL,	0, NULL, 0 }
      };

      c = getopt_long (argc, argv, "hVw:o:d:G:n:ab:zZ:mt:xur:v:", long_options, &option_index);

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
        
        
          if (update_arg( (void *)&(args_info->answers_arg), 
               &(args_info->answers_orig), &(args_info->answers_given),
              &(local_args_info.answers_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "answers", 'w',
              additional_error))
            goto failure;
        
          break;
        case 'o':	/* Output directory.  */
        
        
          if (update_arg( (void *)&(args_info->output_arg), 
               &(args_info->output_orig), &(args_info->output_given),
              &(local_args_info.output_given), optarg, 0, ".", ARG_STRING,
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
        case 'G':	/* Gene exclusion file.  */
        
        
          if (update_arg( (void *)&(args_info->genex_arg), 
               &(args_info->genex_orig), &(args_info->genex_given),
              &(local_args_info.genex_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "genex", 'G',
              additional_error))
            goto failure;
        
          break;
        case 'n':	/* Gene set for negative pairs.  */
        
        
          if (update_arg( (void *)&(args_info->negatives_arg), 
               &(args_info->negatives_orig), &(args_info->negatives_given),
              &(local_args_info.negatives_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "negatives", 'n',
              additional_error))
            goto failure;
        
          break;
        case 'a':	/* Randomize data before training.  */
        
        
          if (update_arg((void *)&(args_info->randomize_flag), 0, &(args_info->randomize_given),
              &(local_args_info.randomize_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "randomize", 'a',
              additional_error))
            goto failure;
        
          break;
        case 'b':	/* Bayes net containing defaults for cases with missing data.  */
        
        
          if (update_arg( (void *)&(args_info->default_arg), 
               &(args_info->default_orig), &(args_info->default_given),
              &(local_args_info.default_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "default", 'b',
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
        case 'Z':	/* Read zeroed node IDs/outputs from the given file.  */
        
        
          if (update_arg( (void *)&(args_info->zeros_arg), 
               &(args_info->zeros_orig), &(args_info->zeros_given),
              &(local_args_info.zeros_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "zeros", 'Z',
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
        case 't':	/* Maximum number of threads to spawn.  */
        
        
          if (update_arg( (void *)&(args_info->threads_arg), 
               &(args_info->threads_orig), &(args_info->threads_given),
              &(local_args_info.threads_given), optarg, 0, "-1", ARG_INT,
              check_ambiguity, override, 0, 0,
              "threads", 't',
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
        case 'u':	/* Group identical inputs.  */
        
        
          if (update_arg((void *)&(args_info->group_flag), 0, &(args_info->group_given),
              &(local_args_info.group_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "group", 'u',
              additional_error))
            goto failure;
        
          break;
        case 'r':	/* Seed random generator.  */
        
        
          if (update_arg( (void *)&(args_info->random_arg), 
               &(args_info->random_orig), &(args_info->random_given),
              &(local_args_info.random_given), optarg, 0, "0", ARG_INT,
              check_ambiguity, override, 0, 0,
              "random", 'r',
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
