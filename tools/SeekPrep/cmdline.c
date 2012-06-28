/*
  File autogenerated by gengetopt version 2.22.5
  generated with the following command:
  /usr/local/bin/gengetopt -iSeekPrep.ggo --default-optional -u -N -e 

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

const char *gengetopt_args_info_purpose = "Preprocess datasets for Seek";

const char *gengetopt_args_info_usage = "Usage: SeekPrep [OPTIONS]... [FILES]...";

const char *gengetopt_args_info_description = "";

const char *gengetopt_args_info_help[] = {
  "  -h, --help               Print help and exit",
  "  -V, --version            Print version and exit",
  "\nMain:",
  "  -x, --dab=filename       Input dataset (.dab) or databaselet (.db)",
  "  -s, --sinfo              Generates sinfo file (with dataset wide mean and \n                             stdev)  (default=off)",
  "  -a, --gavg               Generates gene-average file  (default=off)",
  "  -P, --gplat              Generates gene-platform average + stdev file \n                             (requires .db input)  (default=off)",
  "  -p, --gpres              Generates gene-presence file  (default=off)",
  "  -v, --gvar               Generates gene-variance file  (default=off)",
  "  -D, --dir_out=directory  Output directory",
  "  -i, --input=filename     Gene mapping file",
  "  -A, --dset=filename      Dataset ordering file (with platform info) (required \n                             for -P)",
  "  -N, --useNibble          Is DB file nibble? (required for -P)  (default=off)",
    0
};

typedef enum {ARG_NO
  , ARG_FLAG
  , ARG_STRING
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
  args_info->dab_given = 0 ;
  args_info->sinfo_given = 0 ;
  args_info->gavg_given = 0 ;
  args_info->gplat_given = 0 ;
  args_info->gpres_given = 0 ;
  args_info->gvar_given = 0 ;
  args_info->dir_out_given = 0 ;
  args_info->input_given = 0 ;
  args_info->dset_given = 0 ;
  args_info->useNibble_given = 0 ;
}

static
void clear_args (struct gengetopt_args_info *args_info)
{
  FIX_UNUSED (args_info);
  args_info->dab_arg = NULL;
  args_info->dab_orig = NULL;
  args_info->sinfo_flag = 0;
  args_info->gavg_flag = 0;
  args_info->gplat_flag = 0;
  args_info->gpres_flag = 0;
  args_info->gvar_flag = 0;
  args_info->dir_out_arg = NULL;
  args_info->dir_out_orig = NULL;
  args_info->input_arg = NULL;
  args_info->input_orig = NULL;
  args_info->dset_arg = NULL;
  args_info->dset_orig = NULL;
  args_info->useNibble_flag = 0;
  
}

static
void init_args_info(struct gengetopt_args_info *args_info)
{


  args_info->help_help = gengetopt_args_info_help[0] ;
  args_info->version_help = gengetopt_args_info_help[1] ;
  args_info->dab_help = gengetopt_args_info_help[3] ;
  args_info->sinfo_help = gengetopt_args_info_help[4] ;
  args_info->gavg_help = gengetopt_args_info_help[5] ;
  args_info->gplat_help = gengetopt_args_info_help[6] ;
  args_info->gpres_help = gengetopt_args_info_help[7] ;
  args_info->gvar_help = gengetopt_args_info_help[8] ;
  args_info->dir_out_help = gengetopt_args_info_help[9] ;
  args_info->input_help = gengetopt_args_info_help[10] ;
  args_info->dset_help = gengetopt_args_info_help[11] ;
  args_info->useNibble_help = gengetopt_args_info_help[12] ;
  
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
  free_string_field (&(args_info->dab_arg));
  free_string_field (&(args_info->dab_orig));
  free_string_field (&(args_info->dir_out_arg));
  free_string_field (&(args_info->dir_out_orig));
  free_string_field (&(args_info->input_arg));
  free_string_field (&(args_info->input_orig));
  free_string_field (&(args_info->dset_arg));
  free_string_field (&(args_info->dset_orig));
  
  
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
  if (args_info->dab_given)
    write_into_file(outfile, "dab", args_info->dab_orig, 0);
  if (args_info->sinfo_given)
    write_into_file(outfile, "sinfo", 0, 0 );
  if (args_info->gavg_given)
    write_into_file(outfile, "gavg", 0, 0 );
  if (args_info->gplat_given)
    write_into_file(outfile, "gplat", 0, 0 );
  if (args_info->gpres_given)
    write_into_file(outfile, "gpres", 0, 0 );
  if (args_info->gvar_given)
    write_into_file(outfile, "gvar", 0, 0 );
  if (args_info->dir_out_given)
    write_into_file(outfile, "dir_out", args_info->dir_out_orig, 0);
  if (args_info->input_given)
    write_into_file(outfile, "input", args_info->input_orig, 0);
  if (args_info->dset_given)
    write_into_file(outfile, "dset", args_info->dset_orig, 0);
  if (args_info->useNibble_given)
    write_into_file(outfile, "useNibble", 0, 0 );
  

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
  if (! args_info->dab_given)
    {
      fprintf (stderr, "%s: '--dab' ('-x') option required%s\n", prog_name, (additional_error ? additional_error : ""));
      error = 1;
    }
  
  if (! args_info->dir_out_given)
    {
      fprintf (stderr, "%s: '--dir_out' ('-D') option required%s\n", prog_name, (additional_error ? additional_error : ""));
      error = 1;
    }
  
  if (! args_info->input_given)
    {
      fprintf (stderr, "%s: '--input' ('-i') option required%s\n", prog_name, (additional_error ? additional_error : ""));
      error = 1;
    }
  
  if (! args_info->dset_given)
    {
      fprintf (stderr, "%s: '--dset' ('-A') option required%s\n", prog_name, (additional_error ? additional_error : ""));
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
        { "dab",	1, NULL, 'x' },
        { "sinfo",	0, NULL, 's' },
        { "gavg",	0, NULL, 'a' },
        { "gplat",	0, NULL, 'P' },
        { "gpres",	0, NULL, 'p' },
        { "gvar",	0, NULL, 'v' },
        { "dir_out",	1, NULL, 'D' },
        { "input",	1, NULL, 'i' },
        { "dset",	1, NULL, 'A' },
        { "useNibble",	0, NULL, 'N' },
        { 0,  0, 0, 0 }
      };

      c = getopt_long (argc, argv, "hVx:saPpvD:i:A:N", long_options, &option_index);

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
        case 'x':	/* Input dataset (.dab) or databaselet (.db).  */
        
        
          if (update_arg( (void *)&(args_info->dab_arg), 
               &(args_info->dab_orig), &(args_info->dab_given),
              &(local_args_info.dab_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "dab", 'x',
              additional_error))
            goto failure;
        
          break;
        case 's':	/* Generates sinfo file (with dataset wide mean and stdev).  */
        
        
          if (update_arg((void *)&(args_info->sinfo_flag), 0, &(args_info->sinfo_given),
              &(local_args_info.sinfo_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "sinfo", 's',
              additional_error))
            goto failure;
        
          break;
        case 'a':	/* Generates gene-average file.  */
        
        
          if (update_arg((void *)&(args_info->gavg_flag), 0, &(args_info->gavg_given),
              &(local_args_info.gavg_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "gavg", 'a',
              additional_error))
            goto failure;
        
          break;
        case 'P':	/* Generates gene-platform average + stdev file (requires .db input).  */
        
        
          if (update_arg((void *)&(args_info->gplat_flag), 0, &(args_info->gplat_given),
              &(local_args_info.gplat_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "gplat", 'P',
              additional_error))
            goto failure;
        
          break;
        case 'p':	/* Generates gene-presence file.  */
        
        
          if (update_arg((void *)&(args_info->gpres_flag), 0, &(args_info->gpres_given),
              &(local_args_info.gpres_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "gpres", 'p',
              additional_error))
            goto failure;
        
          break;
        case 'v':	/* Generates gene-variance file.  */
        
        
          if (update_arg((void *)&(args_info->gvar_flag), 0, &(args_info->gvar_given),
              &(local_args_info.gvar_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "gvar", 'v',
              additional_error))
            goto failure;
        
          break;
        case 'D':	/* Output directory.  */
        
        
          if (update_arg( (void *)&(args_info->dir_out_arg), 
               &(args_info->dir_out_orig), &(args_info->dir_out_given),
              &(local_args_info.dir_out_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "dir_out", 'D',
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
        case 'A':	/* Dataset ordering file (with platform info) (required for -P).  */
        
        
          if (update_arg( (void *)&(args_info->dset_arg), 
               &(args_info->dset_orig), &(args_info->dset_given),
              &(local_args_info.dset_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "dset", 'A',
              additional_error))
            goto failure;
        
          break;
        case 'N':	/* Is DB file nibble? (required for -P).  */
        
        
          if (update_arg((void *)&(args_info->useNibble_flag), 0, &(args_info->useNibble_given),
              &(local_args_info.useNibble_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "useNibble", 'N',
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
