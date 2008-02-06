/*
  File autogenerated by gengetopt version 2.13.1
  generated with the following command:
  ..\gengetopt-2.13.1\Release\gengetopt.exe -u -N -e -i Combiner.ggo 

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

static
void clear_given (struct gengetopt_args_info *args_info);
static
void clear_args (struct gengetopt_args_info *args_info);

static int
cmdline_parser_internal (int argc, char * const *argv, struct gengetopt_args_info *args_info, int override, int initialize, int check_required, const char *additional_error);


char *type_values[] = {"pcl", "dat", "dab", 0} ;	/* Possible values for type.  */
char *method_values[] = {"min", "max", "mean", "gmean", "hmean", 0} ;	/* Possible values for method.  */

static char *
gengetopt_strdup (const char *s);

static
void clear_given (struct gengetopt_args_info *args_info)
{
  args_info->help_given = 0 ;
  args_info->version_given = 0 ;
  args_info->type_given = 0 ;
  args_info->method_given = 0 ;
  args_info->output_given = 0 ;
  args_info->subset_given = 0 ;
  args_info->skip_given = 0 ;
  args_info->memmap_given = 0 ;
  args_info->verbosity_given = 0 ;
  args_info->normalize_given = 0 ;
}

static
void clear_args (struct gengetopt_args_info *args_info)
{
  args_info->type_arg = gengetopt_strdup ("pcl");
  args_info->method_arg = gengetopt_strdup ("mean");
  args_info->output_arg = NULL;
  args_info->subset_arg = 0 ;
  args_info->skip_arg = 2 ;
  args_info->memmap_flag = 0;
  args_info->verbosity_arg = 5 ;
  args_info->normalize_flag = 0;
}

void
cmdline_parser_print_version (void)
{
  printf ("%s %s\n", CMDLINE_PARSER_PACKAGE, CMDLINE_PARSER_VERSION);
}

void
cmdline_parser_print_help (void)
{
  cmdline_parser_print_version ();
  printf("\n%s\n", "PCL and data file combination tool");
  printf("\nUsage: Combiner [OPTIONS]... [FILES]...\n\n");
  printf("%s\n","  -h, --help             Print help and exit");
  printf("%s\n","  -V, --version          Print version and exit");
  printf("%s\n","\nOptional:");
  printf("%s\n","  -t, --type=STRING      Data file type  (possible values=\"pcl\", \"dat\", \n                           \"dab\" default=`pcl')");
  printf("%s\n","  -m, --method=STRING    Combination method  (possible values=\"min\", \"max\", \n                           \"mean\", \"gmean\", \"hmean\" default=`mean')");
  printf("%s\n","  -o, --output=filename  Output file");
  printf("%s\n","  -s, --subset=INT       Subset size (none if zero)  (default=`0')");
  printf("%s\n","  -k, --skip=INT         Columns to skip in input PCLs  (default=`2')");
  printf("%s\n","  -p, --memmap           Memory map input files  (default=off)");
  printf("%s\n","  -v, --verbosity=INT    Message verbosity  (default=`5')");
  printf("%s\n","  -n, --normalize        Normalize inputs before combining  (default=off)");
  
}

void
cmdline_parser_init (struct gengetopt_args_info *args_info)
{
  clear_given (args_info);
  clear_args (args_info);

  args_info->inputs = NULL;
  args_info->inputs_num = 0;
}

void
cmdline_parser_free (struct gengetopt_args_info *args_info)
{
  
  unsigned int i;
  if (args_info->type_arg)
    {
      free (args_info->type_arg); /* free previous argument */
      args_info->type_arg = 0;
    }
  if (args_info->method_arg)
    {
      free (args_info->method_arg); /* free previous argument */
      args_info->method_arg = 0;
    }
  if (args_info->output_arg)
    {
      free (args_info->output_arg); /* free previous argument */
      args_info->output_arg = 0;
    }
  
  for (i = 0; i < args_info->inputs_num; ++i)
    free (args_info->inputs [i]);
  
  if (args_info->inputs_num)
    free (args_info->inputs);
  
  clear_given (args_info);
}

static int
check_possible_values(const char *val, char *values[])
{
  int i, found;
  size_t len;

  if (!val)   /* otherwise strlen() crashes below */
    return 0; /* NULL means no argument for the option */

  for (found = i = 0, len = strlen(val); values[i]; ++i)
    {
      if (strncmp(val, values[i], len) == 0)
        {
          found++;
          if (strlen(values[i]) == len)
            return 1; /* exact macth no need to check more */
        }
    }

  return found; /* return how many values are matched */
}

/* gengetopt_strdup() */
/* strdup.c replacement of strdup, which is not standard */
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
cmdline_parser2 (int argc, char * const *argv, struct gengetopt_args_info *args_info, int override, int initialize, int check_required)
{
  int result;

  result = cmdline_parser_internal (argc, argv, args_info, override, initialize, check_required, NULL);

  return result;
}

int
cmdline_parser_required (struct gengetopt_args_info *args_info, const char *prog_name)
{
  return EXIT_SUCCESS;
}

int
cmdline_parser_internal (int argc, char * const *argv, struct gengetopt_args_info *args_info, int override, int initialize, int check_required, const char *additional_error)
{
  int c;	/* Character of the parsed option.  */

  int error = 0;
  struct gengetopt_args_info local_args_info;

  if (initialize)
    cmdline_parser_init (args_info);

  cmdline_parser_init (&local_args_info);

  optarg = 0;
  optind = 1;
  opterr = 1;
  optopt = '?';

  while (1)
    {
      int found = 0;
      int option_index = 0;
      char *stop_char;

      static struct option long_options[] = {
        { "help",	0, NULL, 'h' },
        { "version",	0, NULL, 'V' },
        { "type",	1, NULL, 't' },
        { "method",	1, NULL, 'm' },
        { "output",	1, NULL, 'o' },
        { "subset",	1, NULL, 's' },
        { "skip",	1, NULL, 'k' },
        { "memmap",	0, NULL, 'p' },
        { "verbosity",	1, NULL, 'v' },
        { "normalize",	0, NULL, 'n' },
        { NULL,	0, NULL, 0 }
      };

      stop_char = 0;
      c = getopt_long (argc, argv, "hVt:m:o:s:k:pv:n", long_options, &option_index);

      if (c == -1) break;	/* Exit from `while (1)' loop.  */

      switch (c)
        {
        case 'h':	/* Print help and exit.  */
          cmdline_parser_print_help ();
          cmdline_parser_free (&local_args_info);
          exit (EXIT_SUCCESS);

        case 'V':	/* Print version and exit.  */
          if (local_args_info.version_given)
            {
              fprintf (stderr, "%s: `--version' (`-V') option given more than once%s\n", argv[0], (additional_error ? additional_error : ""));
              goto failure;
            }
          if (args_info->version_given && ! override)
            continue;
          local_args_info.version_given = 1;
          args_info->version_given = 1;
          cmdline_parser_free (&local_args_info);
          return 0;

        case 't':	/* Data file type.  */
          if (local_args_info.type_given)
            {
              fprintf (stderr, "%s: `--type' (`-t') option given more than once%s\n", argv[0], (additional_error ? additional_error : ""));
              goto failure;
            }
          if ((found = check_possible_values(optarg, type_values)) != 1)
            {
              fprintf (stderr, "%s: %s argument, \"%s\", for option `--type' (`-t')%s\n", argv[0], found ? "ambiguous" : "invalid", optarg, (additional_error ? additional_error : ""));
              goto failure;
            }
          if (args_info->type_given && ! override)
            continue;
          local_args_info.type_given = 1;
          args_info->type_given = 1;
          if (args_info->type_arg)
            free (args_info->type_arg); /* free previous string */
          args_info->type_arg = gengetopt_strdup (optarg);
          break;

        case 'm':	/* Combination method.  */
          if (local_args_info.method_given)
            {
              fprintf (stderr, "%s: `--method' (`-m') option given more than once%s\n", argv[0], (additional_error ? additional_error : ""));
              goto failure;
            }
          if ((found = check_possible_values(optarg, method_values)) != 1)
            {
              fprintf (stderr, "%s: %s argument, \"%s\", for option `--method' (`-m')%s\n", argv[0], found ? "ambiguous" : "invalid", optarg, (additional_error ? additional_error : ""));
              goto failure;
            }
          if (args_info->method_given && ! override)
            continue;
          local_args_info.method_given = 1;
          args_info->method_given = 1;
          if (args_info->method_arg)
            free (args_info->method_arg); /* free previous string */
          args_info->method_arg = gengetopt_strdup (optarg);
          break;

        case 'o':	/* Output file.  */
          if (local_args_info.output_given)
            {
              fprintf (stderr, "%s: `--output' (`-o') option given more than once%s\n", argv[0], (additional_error ? additional_error : ""));
              goto failure;
            }
          if (args_info->output_given && ! override)
            continue;
          local_args_info.output_given = 1;
          args_info->output_given = 1;
          if (args_info->output_arg)
            free (args_info->output_arg); /* free previous string */
          args_info->output_arg = gengetopt_strdup (optarg);
          break;

        case 's':	/* Subset size (none if zero).  */
          if (local_args_info.subset_given)
            {
              fprintf (stderr, "%s: `--subset' (`-s') option given more than once%s\n", argv[0], (additional_error ? additional_error : ""));
              goto failure;
            }
          if (args_info->subset_given && ! override)
            continue;
          local_args_info.subset_given = 1;
          args_info->subset_given = 1;
          args_info->subset_arg = strtol (optarg,&stop_char,0);
          break;

        case 'k':	/* Columns to skip in input PCLs.  */
          if (local_args_info.skip_given)
            {
              fprintf (stderr, "%s: `--skip' (`-k') option given more than once%s\n", argv[0], (additional_error ? additional_error : ""));
              goto failure;
            }
          if (args_info->skip_given && ! override)
            continue;
          local_args_info.skip_given = 1;
          args_info->skip_given = 1;
          args_info->skip_arg = strtol (optarg,&stop_char,0);
          break;

        case 'p':	/* Memory map input files.  */
          if (local_args_info.memmap_given)
            {
              fprintf (stderr, "%s: `--memmap' (`-p') option given more than once%s\n", argv[0], (additional_error ? additional_error : ""));
              goto failure;
            }
          if (args_info->memmap_given && ! override)
            continue;
          local_args_info.memmap_given = 1;
          args_info->memmap_given = 1;
          args_info->memmap_flag = !(args_info->memmap_flag);
          break;

        case 'v':	/* Message verbosity.  */
          if (local_args_info.verbosity_given)
            {
              fprintf (stderr, "%s: `--verbosity' (`-v') option given more than once%s\n", argv[0], (additional_error ? additional_error : ""));
              goto failure;
            }
          if (args_info->verbosity_given && ! override)
            continue;
          local_args_info.verbosity_given = 1;
          args_info->verbosity_given = 1;
          args_info->verbosity_arg = strtol (optarg,&stop_char,0);
          break;

        case 'n':	/* Normalize inputs before combining.  */
          if (local_args_info.normalize_given)
            {
              fprintf (stderr, "%s: `--normalize' (`-n') option given more than once%s\n", argv[0], (additional_error ? additional_error : ""));
              goto failure;
            }
          if (args_info->normalize_given && ! override)
            continue;
          local_args_info.normalize_given = 1;
          args_info->normalize_given = 1;
          args_info->normalize_flag = !(args_info->normalize_flag);
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




  cmdline_parser_free (&local_args_info);

  if ( error )
    return (EXIT_FAILURE);

  if (optind < argc)
    {
      int i = 0 ;

      args_info->inputs_num = argc - optind ;
      args_info->inputs =
        (char **)(malloc ((args_info->inputs_num)*sizeof(char *))) ;
      while (optind < argc)
        args_info->inputs[ i++ ] = gengetopt_strdup (argv[optind++]) ;
    }

  return 0;

failure:
  
  cmdline_parser_free (&local_args_info);
  return (EXIT_FAILURE);
}
