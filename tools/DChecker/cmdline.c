/*
  File autogenerated by gengetopt version 2.13.1
  generated with the following command:
  ..\gengetopt-2.13.1\Release\gengetopt.exe -u -N -e -i DChecker.ggo 

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

static int
cmdline_parser_required2 (struct gengetopt_args_info *args_info, const char *prog_name, const char *additional_error);

static char *
gengetopt_strdup (const char *s);

static
void clear_given (struct gengetopt_args_info *args_info)
{
  args_info->help_given = 0 ;
  args_info->version_given = 0 ;
  args_info->input_given = 0 ;
  args_info->answers_given = 0 ;
  args_info->invert_given = 0 ;
  args_info->normalize_given = 0 ;
  args_info->min_given = 0 ;
  args_info->max_given = 0 ;
  args_info->delta_given = 0 ;
  args_info->bins_given = 0 ;
  args_info->finite_given = 0 ;
  args_info->directory_given = 0 ;
  args_info->genes_given = 0 ;
  args_info->genex_given = 0 ;
  args_info->genet_given = 0 ;
  args_info->genee_given = 0 ;
  args_info->sse_given = 0 ;
  args_info->verbosity_given = 0 ;
  args_info->memmap_given = 0 ;
  args_info->unannotated_given = 0 ;
}

static
void clear_args (struct gengetopt_args_info *args_info)
{
  args_info->input_arg = NULL;
  args_info->answers_arg = NULL;
  args_info->invert_flag = 0;
  args_info->normalize_flag = 0;
  args_info->min_arg = 0 ;
  args_info->max_arg = 1 ;
  args_info->delta_arg = 0.01 ;
  args_info->bins_arg = 0 ;
  args_info->finite_flag = 0;
  args_info->directory_arg = gengetopt_strdup (".");
  args_info->genes_arg = NULL;
  args_info->genex_arg = NULL;
  args_info->genet_arg = NULL;
  args_info->genee_arg = NULL;
  args_info->sse_flag = 0;
  args_info->verbosity_arg = 5 ;
  args_info->memmap_flag = 0;
  args_info->unannotated_flag = 0;
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
  printf("\n%s\n", "Correlation to answer file checker");
  printf("\nUsage: DChecker [OPTIONS]... [FILES]...\n\n");
  printf("%s\n","  -h, --help                 Print help and exit");
  printf("%s\n","  -V, --version              Print version and exit");
  printf("%s\n","  -i, --input=filename       Correlation file");
  printf("%s\n","  -w, --answers=filename     Answer file");
  printf("%s\n","\nOptional:");
  printf("%s\n","  -t, --invert               Invert correlations to distances  (default=off)");
  printf("%s\n","  -n, --normalize            Normlize scores before processing  (default=off)");
  printf("%s\n","  -m, --min=FLOAT            Minimum correlation to process  (default=`0')");
  printf("%s\n","  -M, --max=FLOAT            Maximum correlation to process  (default=`1')");
  printf("%s\n","  -e, --delta=DOUBLE         Size of correlation bins  (default=`0.01')");
  printf("%s\n","  -b, --bins=INT             Bins for quantile sorting  (default=`0')");
  printf("%s\n","  -f, --finite               Count finitely many bins  (default=off)");
  printf("%s\n","  -d, --directory=directory  Output directory  (default=`.')");
  printf("%s\n","  -g, --genes=filename       Gene inclusion file");
  printf("%s\n","  -G, --genex=filename       Gene exclusion file");
  printf("%s\n","  -c, --genet=filename       Term inclusion file");
  printf("%s\n","  -C, --genee=filename       Edge inclusion file");
  printf("%s\n","  -s, --sse                  Calculate sum of squared errors  (default=off)");
  printf("%s\n","  -v, --verbosity=INT        Message verbosity  (default=`5')");
  printf("%s\n","  -p, --memmap               Memory map input DABs  (default=off)");
  printf("%s\n","  -u, --unannotated          Include unannotated genes in standard  \n                               (default=off)");
  
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
  if (args_info->input_arg)
    {
      free (args_info->input_arg); /* free previous argument */
      args_info->input_arg = 0;
    }
  if (args_info->answers_arg)
    {
      free (args_info->answers_arg); /* free previous argument */
      args_info->answers_arg = 0;
    }
  if (args_info->directory_arg)
    {
      free (args_info->directory_arg); /* free previous argument */
      args_info->directory_arg = 0;
    }
  if (args_info->genes_arg)
    {
      free (args_info->genes_arg); /* free previous argument */
      args_info->genes_arg = 0;
    }
  if (args_info->genex_arg)
    {
      free (args_info->genex_arg); /* free previous argument */
      args_info->genex_arg = 0;
    }
  if (args_info->genet_arg)
    {
      free (args_info->genet_arg); /* free previous argument */
      args_info->genet_arg = 0;
    }
  if (args_info->genee_arg)
    {
      free (args_info->genee_arg); /* free previous argument */
      args_info->genee_arg = 0;
    }
  
  for (i = 0; i < args_info->inputs_num; ++i)
    free (args_info->inputs [i]);
  
  if (args_info->inputs_num)
    free (args_info->inputs);
  
  clear_given (args_info);
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
  int result = EXIT_SUCCESS;

  if (cmdline_parser_required2(args_info, prog_name, NULL) > 0)
    result = EXIT_FAILURE;

  return result;
}

int
cmdline_parser_required2 (struct gengetopt_args_info *args_info, const char *prog_name, const char *additional_error)
{
  int error = 0;

  if (! args_info->input_given)
    {
      fprintf (stderr, "%s: '--input' ('-i') option required%s\n", prog_name, (additional_error ? additional_error : ""));
      error = 1;
    }
  if (! args_info->answers_given)
    {
      fprintf (stderr, "%s: '--answers' ('-w') option required%s\n", prog_name, (additional_error ? additional_error : ""));
      error = 1;
    }

  return error;
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
      int option_index = 0;
      char *stop_char;

      static struct option long_options[] = {
        { "help",	0, NULL, 'h' },
        { "version",	0, NULL, 'V' },
        { "input",	1, NULL, 'i' },
        { "answers",	1, NULL, 'w' },
        { "invert",	0, NULL, 't' },
        { "normalize",	0, NULL, 'n' },
        { "min",	1, NULL, 'm' },
        { "max",	1, NULL, 'M' },
        { "delta",	1, NULL, 'e' },
        { "bins",	1, NULL, 'b' },
        { "finite",	0, NULL, 'f' },
        { "directory",	1, NULL, 'd' },
        { "genes",	1, NULL, 'g' },
        { "genex",	1, NULL, 'G' },
        { "genet",	1, NULL, 'c' },
        { "genee",	1, NULL, 'C' },
        { "sse",	0, NULL, 's' },
        { "verbosity",	1, NULL, 'v' },
        { "memmap",	0, NULL, 'p' },
        { "unannotated",	0, NULL, 'u' },
        { NULL,	0, NULL, 0 }
      };

      stop_char = 0;
      c = getopt_long (argc, argv, "hVi:w:tnm:M:e:b:fd:g:G:c:C:sv:pu", long_options, &option_index);

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

        case 'i':	/* Correlation file.  */
          if (local_args_info.input_given)
            {
              fprintf (stderr, "%s: `--input' (`-i') option given more than once%s\n", argv[0], (additional_error ? additional_error : ""));
              goto failure;
            }
          if (args_info->input_given && ! override)
            continue;
          local_args_info.input_given = 1;
          args_info->input_given = 1;
          if (args_info->input_arg)
            free (args_info->input_arg); /* free previous string */
          args_info->input_arg = gengetopt_strdup (optarg);
          break;

        case 'w':	/* Answer file.  */
          if (local_args_info.answers_given)
            {
              fprintf (stderr, "%s: `--answers' (`-w') option given more than once%s\n", argv[0], (additional_error ? additional_error : ""));
              goto failure;
            }
          if (args_info->answers_given && ! override)
            continue;
          local_args_info.answers_given = 1;
          args_info->answers_given = 1;
          if (args_info->answers_arg)
            free (args_info->answers_arg); /* free previous string */
          args_info->answers_arg = gengetopt_strdup (optarg);
          break;

        case 't':	/* Invert correlations to distances.  */
          if (local_args_info.invert_given)
            {
              fprintf (stderr, "%s: `--invert' (`-t') option given more than once%s\n", argv[0], (additional_error ? additional_error : ""));
              goto failure;
            }
          if (args_info->invert_given && ! override)
            continue;
          local_args_info.invert_given = 1;
          args_info->invert_given = 1;
          args_info->invert_flag = !(args_info->invert_flag);
          break;

        case 'n':	/* Normlize scores before processing.  */
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

        case 'm':	/* Minimum correlation to process.  */
          if (local_args_info.min_given)
            {
              fprintf (stderr, "%s: `--min' (`-m') option given more than once%s\n", argv[0], (additional_error ? additional_error : ""));
              goto failure;
            }
          if (args_info->min_given && ! override)
            continue;
          local_args_info.min_given = 1;
          args_info->min_given = 1;
          args_info->min_arg = (float)strtod (optarg, NULL);
          break;

        case 'M':	/* Maximum correlation to process.  */
          if (local_args_info.max_given)
            {
              fprintf (stderr, "%s: `--max' (`-M') option given more than once%s\n", argv[0], (additional_error ? additional_error : ""));
              goto failure;
            }
          if (args_info->max_given && ! override)
            continue;
          local_args_info.max_given = 1;
          args_info->max_given = 1;
          args_info->max_arg = (float)strtod (optarg, NULL);
          break;

        case 'e':	/* Size of correlation bins.  */
          if (local_args_info.delta_given)
            {
              fprintf (stderr, "%s: `--delta' (`-e') option given more than once%s\n", argv[0], (additional_error ? additional_error : ""));
              goto failure;
            }
          if (args_info->delta_given && ! override)
            continue;
          local_args_info.delta_given = 1;
          args_info->delta_given = 1;
          args_info->delta_arg = strtod (optarg, NULL);
          break;

        case 'b':	/* Bins for quantile sorting.  */
          if (local_args_info.bins_given)
            {
              fprintf (stderr, "%s: `--bins' (`-b') option given more than once%s\n", argv[0], (additional_error ? additional_error : ""));
              goto failure;
            }
          if (args_info->bins_given && ! override)
            continue;
          local_args_info.bins_given = 1;
          args_info->bins_given = 1;
          args_info->bins_arg = strtol (optarg,&stop_char,0);
          break;

        case 'f':	/* Count finitely many bins.  */
          if (local_args_info.finite_given)
            {
              fprintf (stderr, "%s: `--finite' (`-f') option given more than once%s\n", argv[0], (additional_error ? additional_error : ""));
              goto failure;
            }
          if (args_info->finite_given && ! override)
            continue;
          local_args_info.finite_given = 1;
          args_info->finite_given = 1;
          args_info->finite_flag = !(args_info->finite_flag);
          break;

        case 'd':	/* Output directory.  */
          if (local_args_info.directory_given)
            {
              fprintf (stderr, "%s: `--directory' (`-d') option given more than once%s\n", argv[0], (additional_error ? additional_error : ""));
              goto failure;
            }
          if (args_info->directory_given && ! override)
            continue;
          local_args_info.directory_given = 1;
          args_info->directory_given = 1;
          if (args_info->directory_arg)
            free (args_info->directory_arg); /* free previous string */
          args_info->directory_arg = gengetopt_strdup (optarg);
          break;

        case 'g':	/* Gene inclusion file.  */
          if (local_args_info.genes_given)
            {
              fprintf (stderr, "%s: `--genes' (`-g') option given more than once%s\n", argv[0], (additional_error ? additional_error : ""));
              goto failure;
            }
          if (args_info->genes_given && ! override)
            continue;
          local_args_info.genes_given = 1;
          args_info->genes_given = 1;
          if (args_info->genes_arg)
            free (args_info->genes_arg); /* free previous string */
          args_info->genes_arg = gengetopt_strdup (optarg);
          break;

        case 'G':	/* Gene exclusion file.  */
          if (local_args_info.genex_given)
            {
              fprintf (stderr, "%s: `--genex' (`-G') option given more than once%s\n", argv[0], (additional_error ? additional_error : ""));
              goto failure;
            }
          if (args_info->genex_given && ! override)
            continue;
          local_args_info.genex_given = 1;
          args_info->genex_given = 1;
          if (args_info->genex_arg)
            free (args_info->genex_arg); /* free previous string */
          args_info->genex_arg = gengetopt_strdup (optarg);
          break;

        case 'c':	/* Term inclusion file.  */
          if (local_args_info.genet_given)
            {
              fprintf (stderr, "%s: `--genet' (`-c') option given more than once%s\n", argv[0], (additional_error ? additional_error : ""));
              goto failure;
            }
          if (args_info->genet_given && ! override)
            continue;
          local_args_info.genet_given = 1;
          args_info->genet_given = 1;
          if (args_info->genet_arg)
            free (args_info->genet_arg); /* free previous string */
          args_info->genet_arg = gengetopt_strdup (optarg);
          break;

        case 'C':	/* Edge inclusion file.  */
          if (local_args_info.genee_given)
            {
              fprintf (stderr, "%s: `--genee' (`-C') option given more than once%s\n", argv[0], (additional_error ? additional_error : ""));
              goto failure;
            }
          if (args_info->genee_given && ! override)
            continue;
          local_args_info.genee_given = 1;
          args_info->genee_given = 1;
          if (args_info->genee_arg)
            free (args_info->genee_arg); /* free previous string */
          args_info->genee_arg = gengetopt_strdup (optarg);
          break;

        case 's':	/* Calculate sum of squared errors.  */
          if (local_args_info.sse_given)
            {
              fprintf (stderr, "%s: `--sse' (`-s') option given more than once%s\n", argv[0], (additional_error ? additional_error : ""));
              goto failure;
            }
          if (args_info->sse_given && ! override)
            continue;
          local_args_info.sse_given = 1;
          args_info->sse_given = 1;
          args_info->sse_flag = !(args_info->sse_flag);
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

        case 'p':	/* Memory map input DABs.  */
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

        case 'u':	/* Include unannotated genes in standard.  */
          if (local_args_info.unannotated_given)
            {
              fprintf (stderr, "%s: `--unannotated' (`-u') option given more than once%s\n", argv[0], (additional_error ? additional_error : ""));
              goto failure;
            }
          if (args_info->unannotated_given && ! override)
            continue;
          local_args_info.unannotated_given = 1;
          args_info->unannotated_given = 1;
          args_info->unannotated_flag = !(args_info->unannotated_flag);
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
