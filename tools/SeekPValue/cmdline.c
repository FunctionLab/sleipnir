/*
  File autogenerated by gengetopt version 2.22.5
  generated with the following command:
  /usr/bin/gengetopt -iSeekPValue.ggo --default-optional -u -N -e 

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

const char *gengetopt_args_info_purpose = "Estimates P-Value of the retrieved genes based on a background of random \nqueries";

const char *gengetopt_args_info_usage = "Usage: SeekPValue [OPTIONS]... [FILES]...";

const char *gengetopt_args_info_description = "";

const char *gengetopt_args_info_help[] = {
  "  -h, --help                    Print help and exit",
  "  -V, --version                 Print version and exit",
  "\nMain:",
  "  -m, --mode=STRING             Mode (datasets or genes)  (possible \n                                  values=\"datasets\", \"genes\" \n                                  default=`genes')",
  "  -t, --port=STRING             Port  (default=`9005')",
  "\nDataset mode:",
  "  -d, --dset_platform=filename  Dataset platform file",
  "  -p, --param_dir=directory     Parameter directory",
  "\nGene mode:",
  "  -R, --random_dir=directory    Random directory",
  "  -N, --random_num=INT          Number of random trials  (default=`100')",
  "  -i, --input=filename          Gene mapping file",
  "  -n, --nan=FLOAT               Define NaN score  (default=`-320')",
    0
};

typedef enum {ARG_NO
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


const char *cmdline_parser_mode_values[] = {"datasets", "genes", 0}; /*< Possible values for mode. */

static char *
gengetopt_strdup (const char *s);

static
void clear_given (struct gengetopt_args_info *args_info)
{
  args_info->help_given = 0 ;
  args_info->version_given = 0 ;
  args_info->mode_given = 0 ;
  args_info->port_given = 0 ;
  args_info->dset_platform_given = 0 ;
  args_info->param_dir_given = 0 ;
  args_info->random_dir_given = 0 ;
  args_info->random_num_given = 0 ;
  args_info->input_given = 0 ;
  args_info->nan_given = 0 ;
}

static
void clear_args (struct gengetopt_args_info *args_info)
{
  FIX_UNUSED (args_info);
  args_info->mode_arg = gengetopt_strdup ("genes");
  args_info->mode_orig = NULL;
  args_info->port_arg = gengetopt_strdup ("9005");
  args_info->port_orig = NULL;
  args_info->dset_platform_arg = NULL;
  args_info->dset_platform_orig = NULL;
  args_info->param_dir_arg = NULL;
  args_info->param_dir_orig = NULL;
  args_info->random_dir_arg = NULL;
  args_info->random_dir_orig = NULL;
  args_info->random_num_arg = 100;
  args_info->random_num_orig = NULL;
  args_info->input_arg = NULL;
  args_info->input_orig = NULL;
  args_info->nan_arg = -320;
  args_info->nan_orig = NULL;
  
}

static
void init_args_info(struct gengetopt_args_info *args_info)
{


  args_info->help_help = gengetopt_args_info_help[0] ;
  args_info->version_help = gengetopt_args_info_help[1] ;
  args_info->mode_help = gengetopt_args_info_help[3] ;
  args_info->port_help = gengetopt_args_info_help[4] ;
  args_info->dset_platform_help = gengetopt_args_info_help[6] ;
  args_info->param_dir_help = gengetopt_args_info_help[7] ;
  args_info->random_dir_help = gengetopt_args_info_help[9] ;
  args_info->random_num_help = gengetopt_args_info_help[10] ;
  args_info->input_help = gengetopt_args_info_help[11] ;
  args_info->nan_help = gengetopt_args_info_help[12] ;
  
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
  free_string_field (&(args_info->mode_arg));
  free_string_field (&(args_info->mode_orig));
  free_string_field (&(args_info->port_arg));
  free_string_field (&(args_info->port_orig));
  free_string_field (&(args_info->dset_platform_arg));
  free_string_field (&(args_info->dset_platform_orig));
  free_string_field (&(args_info->param_dir_arg));
  free_string_field (&(args_info->param_dir_orig));
  free_string_field (&(args_info->random_dir_arg));
  free_string_field (&(args_info->random_dir_orig));
  free_string_field (&(args_info->random_num_orig));
  free_string_field (&(args_info->input_arg));
  free_string_field (&(args_info->input_orig));
  free_string_field (&(args_info->nan_orig));
  
  
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
check_possible_values(const char *val, const char *values[])
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
write_into_file(FILE *outfile, const char *opt, const char *arg, const char *values[])
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
  if (args_info->mode_given)
    write_into_file(outfile, "mode", args_info->mode_orig, cmdline_parser_mode_values);
  if (args_info->port_given)
    write_into_file(outfile, "port", args_info->port_orig, 0);
  if (args_info->dset_platform_given)
    write_into_file(outfile, "dset_platform", args_info->dset_platform_orig, 0);
  if (args_info->param_dir_given)
    write_into_file(outfile, "param_dir", args_info->param_dir_orig, 0);
  if (args_info->random_dir_given)
    write_into_file(outfile, "random_dir", args_info->random_dir_orig, 0);
  if (args_info->random_num_given)
    write_into_file(outfile, "random_num", args_info->random_num_orig, 0);
  if (args_info->input_given)
    write_into_file(outfile, "input", args_info->input_orig, 0);
  if (args_info->nan_given)
    write_into_file(outfile, "nan", args_info->nan_orig, 0);
  

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
  FIX_UNUSED (args_info);
  FIX_UNUSED (prog_name);
  return EXIT_SUCCESS;
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
        { "mode",	1, NULL, 'm' },
        { "port",	1, NULL, 't' },
        { "dset_platform",	1, NULL, 'd' },
        { "param_dir",	1, NULL, 'p' },
        { "random_dir",	1, NULL, 'R' },
        { "random_num",	1, NULL, 'N' },
        { "input",	1, NULL, 'i' },
        { "nan",	1, NULL, 'n' },
        { 0,  0, 0, 0 }
      };

      c = getopt_long (argc, argv, "hVm:t:d:p:R:N:i:n:", long_options, &option_index);

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
        case 'm':	/* Mode (datasets or genes).  */
        
        
          if (update_arg( (void *)&(args_info->mode_arg), 
               &(args_info->mode_orig), &(args_info->mode_given),
              &(local_args_info.mode_given), optarg, cmdline_parser_mode_values, "genes", ARG_STRING,
              check_ambiguity, override, 0, 0,
              "mode", 'm',
              additional_error))
            goto failure;
        
          break;
        case 't':	/* Port.  */
        
        
          if (update_arg( (void *)&(args_info->port_arg), 
               &(args_info->port_orig), &(args_info->port_given),
              &(local_args_info.port_given), optarg, 0, "9005", ARG_STRING,
              check_ambiguity, override, 0, 0,
              "port", 't',
              additional_error))
            goto failure;
        
          break;
        case 'd':	/* Dataset platform file.  */
        
        
          if (update_arg( (void *)&(args_info->dset_platform_arg), 
               &(args_info->dset_platform_orig), &(args_info->dset_platform_given),
              &(local_args_info.dset_platform_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "dset_platform", 'd',
              additional_error))
            goto failure;
        
          break;
        case 'p':	/* Parameter directory.  */
        
        
          if (update_arg( (void *)&(args_info->param_dir_arg), 
               &(args_info->param_dir_orig), &(args_info->param_dir_given),
              &(local_args_info.param_dir_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "param_dir", 'p',
              additional_error))
            goto failure;
        
          break;
        case 'R':	/* Random directory.  */
        
        
          if (update_arg( (void *)&(args_info->random_dir_arg), 
               &(args_info->random_dir_orig), &(args_info->random_dir_given),
              &(local_args_info.random_dir_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "random_dir", 'R',
              additional_error))
            goto failure;
        
          break;
        case 'N':	/* Number of random trials.  */
        
        
          if (update_arg( (void *)&(args_info->random_num_arg), 
               &(args_info->random_num_orig), &(args_info->random_num_given),
              &(local_args_info.random_num_given), optarg, 0, "100", ARG_INT,
              check_ambiguity, override, 0, 0,
              "random_num", 'N',
              additional_error))
            goto failure;
        
          break;
        case 'i':	/* Gene mapping file.  */
        
        
          if (update_arg( (void *)&(args_info->input_arg), 
               &(args_info->input_orig), &(args_info->input_given),
              &(local_args_info.input_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "input", 'i',
              additional_error))
            goto failure;
        
          break;
        case 'n':	/* Define NaN score.  */
        
        
          if (update_arg( (void *)&(args_info->nan_arg), 
               &(args_info->nan_orig), &(args_info->nan_given),
              &(local_args_info.nan_given), optarg, 0, "-320", ARG_FLOAT,
              check_ambiguity, override, 0, 0,
              "nan", 'n',
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
