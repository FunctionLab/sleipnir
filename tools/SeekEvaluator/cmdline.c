/*
  File autogenerated by gengetopt version 2.22
  generated with the following command:
  /memex/qzhu/usr/bin/gengetopt -iSeekEvaluator.ggo --default-optional -u -N -e 

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

const char *gengetopt_args_info_purpose = "Evaluate results generated from SeekMiner";

const char *gengetopt_args_info_usage = "Usage: SeekEvaluator [OPTIONS]... [FILES]...";

const char *gengetopt_args_info_description = "";

const char *gengetopt_args_info_help[] = {
  "  -h, --help                   Print help and exit",
  "  -V, --version                Print version and exit",
  "\nMode:",
  "  -O, --single                 Evaluate one query's rank list result  \n                                 (default=on)",
  "  -M, --aggregate              Evaluate multiple queries and aggregates results \n                                  (default=off)",
  "\nMetric:",
  "  -r, --rbp                    Rank biased precision (requires parameter p to \n                                 be set)  (default=off)",
  "  -a, --avgp                   Average precision for the top X positives, where \n                                 X = integer, or % of total positives  \n                                 (default=off)",
  "  -t, --pr                     Precision at X-th positive, where X = integer, \n                                 or % of total positives  (default=off)",
  "  -c, --pr_all                 Precision at all positive depths (useful for \n                                 drawing precision-recall curve)  (default=off)",
  "  -u, --auc                    AUC  (default=off)",
  "  -x, --x_int=INT              Parameter X = integer, for --avgp, --pr  \n                                 (default=`-1')",
  "  -e, --x_per=FLOAT            Parameter X = percentage, for --avgp, --pr  \n                                 (default=`0')",
  "  -p, --rbp_p=FLOAT            Parameter p, for --rbp  (default=`0.95')",
  "\nAggregator (for multi-query evaluation):",
  "  -A, --agg_avg                Show the average, standard deviation of the \n                                 metric for all queries  (default=off)",
  "  -B, --agg_quartile           Show the min, max, as well as the 1st, 2nd, 3rd \n                                 quartile of the metric for all queries  \n                                 (default=off)",
  "  -C, --agg_ranksum            Sum up the ranks of genes in all query rankings \n                                 to produce a master list sorted by summed \n                                 rank, and perform metric on this list  \n                                 (default=off)",
  "  -D, --agg_scoresum           Sum up the scores of genes in all query rankings \n                                 to produce a master list sorted by summed \n                                 score, and perform metric on this list  \n                                 (default=off)",
  "  -E, --display_all            Display the metric for all queries  \n                                 (default=off)",
  "\nInput required by all:",
  "  -i, --input=filename         Gene mapping file",
  "\nInput required by single-query:",
  "  -s, --goldstd=filename       Gold standard gene set file (one line, space \n                                 delimited)",
  "  -g, --gscore=filename        Gene score file (.gscore)",
  "  -q, --query=filename         Query gene set file (to be excluded from \n                                 evaluation) (.query)",
  "  -n, --nan=FLOAT              Define NaN score  (default=`-320')",
  "\nInput required by multi-query:",
  "  -S, --goldstd_list=filename  List of gold standard gene set files",
  "  -G, --gscore_list=filename   List of gene score files",
  "  -Q, --query_list=filename    List of query gene set files",
  "\nOutput:",
  "  -d, --dir_out=directory      Output directory",
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
  args_info->single_given = 0 ;
  args_info->aggregate_given = 0 ;
  args_info->rbp_given = 0 ;
  args_info->avgp_given = 0 ;
  args_info->pr_given = 0 ;
  args_info->pr_all_given = 0 ;
  args_info->auc_given = 0 ;
  args_info->x_int_given = 0 ;
  args_info->x_per_given = 0 ;
  args_info->rbp_p_given = 0 ;
  args_info->agg_avg_given = 0 ;
  args_info->agg_quartile_given = 0 ;
  args_info->agg_ranksum_given = 0 ;
  args_info->agg_scoresum_given = 0 ;
  args_info->display_all_given = 0 ;
  args_info->input_given = 0 ;
  args_info->goldstd_given = 0 ;
  args_info->gscore_given = 0 ;
  args_info->query_given = 0 ;
  args_info->nan_given = 0 ;
  args_info->goldstd_list_given = 0 ;
  args_info->gscore_list_given = 0 ;
  args_info->query_list_given = 0 ;
  args_info->dir_out_given = 0 ;
}

static
void clear_args (struct gengetopt_args_info *args_info)
{
  args_info->single_flag = 1;
  args_info->aggregate_flag = 0;
  args_info->rbp_flag = 0;
  args_info->avgp_flag = 0;
  args_info->pr_flag = 0;
  args_info->pr_all_flag = 0;
  args_info->auc_flag = 0;
  args_info->x_int_arg = -1;
  args_info->x_int_orig = NULL;
  args_info->x_per_arg = 0;
  args_info->x_per_orig = NULL;
  args_info->rbp_p_arg = 0.95;
  args_info->rbp_p_orig = NULL;
  args_info->agg_avg_flag = 0;
  args_info->agg_quartile_flag = 0;
  args_info->agg_ranksum_flag = 0;
  args_info->agg_scoresum_flag = 0;
  args_info->display_all_flag = 0;
  args_info->input_arg = NULL;
  args_info->input_orig = NULL;
  args_info->goldstd_arg = NULL;
  args_info->goldstd_orig = NULL;
  args_info->gscore_arg = NULL;
  args_info->gscore_orig = NULL;
  args_info->query_arg = NULL;
  args_info->query_orig = NULL;
  args_info->nan_arg = -320;
  args_info->nan_orig = NULL;
  args_info->goldstd_list_arg = NULL;
  args_info->goldstd_list_orig = NULL;
  args_info->gscore_list_arg = NULL;
  args_info->gscore_list_orig = NULL;
  args_info->query_list_arg = NULL;
  args_info->query_list_orig = NULL;
  args_info->dir_out_arg = NULL;
  args_info->dir_out_orig = NULL;
  
}

static
void init_args_info(struct gengetopt_args_info *args_info)
{


  args_info->help_help = gengetopt_args_info_help[0] ;
  args_info->version_help = gengetopt_args_info_help[1] ;
  args_info->single_help = gengetopt_args_info_help[3] ;
  args_info->aggregate_help = gengetopt_args_info_help[4] ;
  args_info->rbp_help = gengetopt_args_info_help[6] ;
  args_info->avgp_help = gengetopt_args_info_help[7] ;
  args_info->pr_help = gengetopt_args_info_help[8] ;
  args_info->pr_all_help = gengetopt_args_info_help[9] ;
  args_info->auc_help = gengetopt_args_info_help[10] ;
  args_info->x_int_help = gengetopt_args_info_help[11] ;
  args_info->x_per_help = gengetopt_args_info_help[12] ;
  args_info->rbp_p_help = gengetopt_args_info_help[13] ;
  args_info->agg_avg_help = gengetopt_args_info_help[15] ;
  args_info->agg_quartile_help = gengetopt_args_info_help[16] ;
  args_info->agg_ranksum_help = gengetopt_args_info_help[17] ;
  args_info->agg_scoresum_help = gengetopt_args_info_help[18] ;
  args_info->display_all_help = gengetopt_args_info_help[19] ;
  args_info->input_help = gengetopt_args_info_help[21] ;
  args_info->goldstd_help = gengetopt_args_info_help[23] ;
  args_info->gscore_help = gengetopt_args_info_help[24] ;
  args_info->query_help = gengetopt_args_info_help[25] ;
  args_info->nan_help = gengetopt_args_info_help[26] ;
  args_info->goldstd_list_help = gengetopt_args_info_help[28] ;
  args_info->gscore_list_help = gengetopt_args_info_help[29] ;
  args_info->query_list_help = gengetopt_args_info_help[30] ;
  args_info->dir_out_help = gengetopt_args_info_help[32] ;
  
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
  free_string_field (&(args_info->x_int_orig));
  free_string_field (&(args_info->x_per_orig));
  free_string_field (&(args_info->rbp_p_orig));
  free_string_field (&(args_info->input_arg));
  free_string_field (&(args_info->input_orig));
  free_string_field (&(args_info->goldstd_arg));
  free_string_field (&(args_info->goldstd_orig));
  free_string_field (&(args_info->gscore_arg));
  free_string_field (&(args_info->gscore_orig));
  free_string_field (&(args_info->query_arg));
  free_string_field (&(args_info->query_orig));
  free_string_field (&(args_info->nan_orig));
  free_string_field (&(args_info->goldstd_list_arg));
  free_string_field (&(args_info->goldstd_list_orig));
  free_string_field (&(args_info->gscore_list_arg));
  free_string_field (&(args_info->gscore_list_orig));
  free_string_field (&(args_info->query_list_arg));
  free_string_field (&(args_info->query_list_orig));
  free_string_field (&(args_info->dir_out_arg));
  free_string_field (&(args_info->dir_out_orig));
  
  
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
  if (args_info->single_given)
    write_into_file(outfile, "single", 0, 0 );
  if (args_info->aggregate_given)
    write_into_file(outfile, "aggregate", 0, 0 );
  if (args_info->rbp_given)
    write_into_file(outfile, "rbp", 0, 0 );
  if (args_info->avgp_given)
    write_into_file(outfile, "avgp", 0, 0 );
  if (args_info->pr_given)
    write_into_file(outfile, "pr", 0, 0 );
  if (args_info->pr_all_given)
    write_into_file(outfile, "pr_all", 0, 0 );
  if (args_info->auc_given)
    write_into_file(outfile, "auc", 0, 0 );
  if (args_info->x_int_given)
    write_into_file(outfile, "x_int", args_info->x_int_orig, 0);
  if (args_info->x_per_given)
    write_into_file(outfile, "x_per", args_info->x_per_orig, 0);
  if (args_info->rbp_p_given)
    write_into_file(outfile, "rbp_p", args_info->rbp_p_orig, 0);
  if (args_info->agg_avg_given)
    write_into_file(outfile, "agg_avg", 0, 0 );
  if (args_info->agg_quartile_given)
    write_into_file(outfile, "agg_quartile", 0, 0 );
  if (args_info->agg_ranksum_given)
    write_into_file(outfile, "agg_ranksum", 0, 0 );
  if (args_info->agg_scoresum_given)
    write_into_file(outfile, "agg_scoresum", 0, 0 );
  if (args_info->display_all_given)
    write_into_file(outfile, "display_all", 0, 0 );
  if (args_info->input_given)
    write_into_file(outfile, "input", args_info->input_orig, 0);
  if (args_info->goldstd_given)
    write_into_file(outfile, "goldstd", args_info->goldstd_orig, 0);
  if (args_info->gscore_given)
    write_into_file(outfile, "gscore", args_info->gscore_orig, 0);
  if (args_info->query_given)
    write_into_file(outfile, "query", args_info->query_orig, 0);
  if (args_info->nan_given)
    write_into_file(outfile, "nan", args_info->nan_orig, 0);
  if (args_info->goldstd_list_given)
    write_into_file(outfile, "goldstd_list", args_info->goldstd_list_orig, 0);
  if (args_info->gscore_list_given)
    write_into_file(outfile, "gscore_list", args_info->gscore_list_orig, 0);
  if (args_info->query_list_given)
    write_into_file(outfile, "query_list", args_info->query_list_orig, 0);
  if (args_info->dir_out_given)
    write_into_file(outfile, "dir_out", args_info->dir_out_orig, 0);
  

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
  if (! args_info->input_given)
    {
      fprintf (stderr, "%s: '--input' ('-i') option required%s\n", prog_name, (additional_error ? additional_error : ""));
      error = 1;
    }
  
  if (! args_info->dir_out_given)
    {
      fprintf (stderr, "%s: '--dir_out' ('-d') option required%s\n", prog_name, (additional_error ? additional_error : ""));
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
        { "single",	0, NULL, 'O' },
        { "aggregate",	0, NULL, 'M' },
        { "rbp",	0, NULL, 'r' },
        { "avgp",	0, NULL, 'a' },
        { "pr",	0, NULL, 't' },
        { "pr_all",	0, NULL, 'c' },
        { "auc",	0, NULL, 'u' },
        { "x_int",	1, NULL, 'x' },
        { "x_per",	1, NULL, 'e' },
        { "rbp_p",	1, NULL, 'p' },
        { "agg_avg",	0, NULL, 'A' },
        { "agg_quartile",	0, NULL, 'B' },
        { "agg_ranksum",	0, NULL, 'C' },
        { "agg_scoresum",	0, NULL, 'D' },
        { "display_all",	0, NULL, 'E' },
        { "input",	1, NULL, 'i' },
        { "goldstd",	1, NULL, 's' },
        { "gscore",	1, NULL, 'g' },
        { "query",	1, NULL, 'q' },
        { "nan",	1, NULL, 'n' },
        { "goldstd_list",	1, NULL, 'S' },
        { "gscore_list",	1, NULL, 'G' },
        { "query_list",	1, NULL, 'Q' },
        { "dir_out",	1, NULL, 'd' },
        { NULL,	0, NULL, 0 }
      };

      c = getopt_long (argc, argv, "hVOMratcux:e:p:ABCDEi:s:g:q:n:S:G:Q:d:", long_options, &option_index);

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
        case 'O':	/* Evaluate one query's rank list result.  */
        
        
          if (update_arg((void *)&(args_info->single_flag), 0, &(args_info->single_given),
              &(local_args_info.single_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "single", 'O',
              additional_error))
            goto failure;
        
          break;
        case 'M':	/* Evaluate multiple queries and aggregates results.  */
        
        
          if (update_arg((void *)&(args_info->aggregate_flag), 0, &(args_info->aggregate_given),
              &(local_args_info.aggregate_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "aggregate", 'M',
              additional_error))
            goto failure;
        
          break;
        case 'r':	/* Rank biased precision (requires parameter p to be set).  */
        
        
          if (update_arg((void *)&(args_info->rbp_flag), 0, &(args_info->rbp_given),
              &(local_args_info.rbp_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "rbp", 'r',
              additional_error))
            goto failure;
        
          break;
        case 'a':	/* Average precision for the top X positives, where X = integer, or % of total positives.  */
        
        
          if (update_arg((void *)&(args_info->avgp_flag), 0, &(args_info->avgp_given),
              &(local_args_info.avgp_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "avgp", 'a',
              additional_error))
            goto failure;
        
          break;
        case 't':	/* Precision at X-th positive, where X = integer, or % of total positives.  */
        
        
          if (update_arg((void *)&(args_info->pr_flag), 0, &(args_info->pr_given),
              &(local_args_info.pr_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "pr", 't',
              additional_error))
            goto failure;
        
          break;
        case 'c':	/* Precision at all positive depths (useful for drawing precision-recall curve).  */
        
        
          if (update_arg((void *)&(args_info->pr_all_flag), 0, &(args_info->pr_all_given),
              &(local_args_info.pr_all_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "pr_all", 'c',
              additional_error))
            goto failure;
        
          break;
        case 'u':	/* AUC.  */
        
        
          if (update_arg((void *)&(args_info->auc_flag), 0, &(args_info->auc_given),
              &(local_args_info.auc_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "auc", 'u',
              additional_error))
            goto failure;
        
          break;
        case 'x':	/* Parameter X = integer, for --avgp, --pr.  */
        
        
          if (update_arg( (void *)&(args_info->x_int_arg), 
               &(args_info->x_int_orig), &(args_info->x_int_given),
              &(local_args_info.x_int_given), optarg, 0, "-1", ARG_INT,
              check_ambiguity, override, 0, 0,
              "x_int", 'x',
              additional_error))
            goto failure;
        
          break;
        case 'e':	/* Parameter X = percentage, for --avgp, --pr.  */
        
        
          if (update_arg( (void *)&(args_info->x_per_arg), 
               &(args_info->x_per_orig), &(args_info->x_per_given),
              &(local_args_info.x_per_given), optarg, 0, "0", ARG_FLOAT,
              check_ambiguity, override, 0, 0,
              "x_per", 'e',
              additional_error))
            goto failure;
        
          break;
        case 'p':	/* Parameter p, for --rbp.  */
        
        
          if (update_arg( (void *)&(args_info->rbp_p_arg), 
               &(args_info->rbp_p_orig), &(args_info->rbp_p_given),
              &(local_args_info.rbp_p_given), optarg, 0, "0.95", ARG_FLOAT,
              check_ambiguity, override, 0, 0,
              "rbp_p", 'p',
              additional_error))
            goto failure;
        
          break;
        case 'A':	/* Show the average, standard deviation of the metric for all queries.  */
        
        
          if (update_arg((void *)&(args_info->agg_avg_flag), 0, &(args_info->agg_avg_given),
              &(local_args_info.agg_avg_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "agg_avg", 'A',
              additional_error))
            goto failure;
        
          break;
        case 'B':	/* Show the min, max, as well as the 1st, 2nd, 3rd quartile of the metric for all queries.  */
        
        
          if (update_arg((void *)&(args_info->agg_quartile_flag), 0, &(args_info->agg_quartile_given),
              &(local_args_info.agg_quartile_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "agg_quartile", 'B',
              additional_error))
            goto failure;
        
          break;
        case 'C':	/* Sum up the ranks of genes in all query rankings to produce a master list sorted by summed rank, and perform metric on this list.  */
        
        
          if (update_arg((void *)&(args_info->agg_ranksum_flag), 0, &(args_info->agg_ranksum_given),
              &(local_args_info.agg_ranksum_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "agg_ranksum", 'C',
              additional_error))
            goto failure;
        
          break;
        case 'D':	/* Sum up the scores of genes in all query rankings to produce a master list sorted by summed score, and perform metric on this list.  */
        
        
          if (update_arg((void *)&(args_info->agg_scoresum_flag), 0, &(args_info->agg_scoresum_given),
              &(local_args_info.agg_scoresum_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "agg_scoresum", 'D',
              additional_error))
            goto failure;
        
          break;
        case 'E':	/* Display the metric for all queries.  */
        
        
          if (update_arg((void *)&(args_info->display_all_flag), 0, &(args_info->display_all_given),
              &(local_args_info.display_all_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "display_all", 'E',
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
        case 's':	/* Gold standard gene set file (one line, space delimited).  */
        
        
          if (update_arg( (void *)&(args_info->goldstd_arg), 
               &(args_info->goldstd_orig), &(args_info->goldstd_given),
              &(local_args_info.goldstd_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "goldstd", 's',
              additional_error))
            goto failure;
        
          break;
        case 'g':	/* Gene score file (.gscore).  */
        
        
          if (update_arg( (void *)&(args_info->gscore_arg), 
               &(args_info->gscore_orig), &(args_info->gscore_given),
              &(local_args_info.gscore_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "gscore", 'g',
              additional_error))
            goto failure;
        
          break;
        case 'q':	/* Query gene set file (to be excluded from evaluation) (.query).  */
        
        
          if (update_arg( (void *)&(args_info->query_arg), 
               &(args_info->query_orig), &(args_info->query_given),
              &(local_args_info.query_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "query", 'q',
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
        case 'S':	/* List of gold standard gene set files.  */
        
        
          if (update_arg( (void *)&(args_info->goldstd_list_arg), 
               &(args_info->goldstd_list_orig), &(args_info->goldstd_list_given),
              &(local_args_info.goldstd_list_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "goldstd_list", 'S',
              additional_error))
            goto failure;
        
          break;
        case 'G':	/* List of gene score files.  */
        
        
          if (update_arg( (void *)&(args_info->gscore_list_arg), 
               &(args_info->gscore_list_orig), &(args_info->gscore_list_given),
              &(local_args_info.gscore_list_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "gscore_list", 'G',
              additional_error))
            goto failure;
        
          break;
        case 'Q':	/* List of query gene set files.  */
        
        
          if (update_arg( (void *)&(args_info->query_list_arg), 
               &(args_info->query_list_orig), &(args_info->query_list_given),
              &(local_args_info.query_list_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "query_list", 'Q',
              additional_error))
            goto failure;
        
          break;
        case 'd':	/* Output directory.  */
        
        
          if (update_arg( (void *)&(args_info->dir_out_arg), 
               &(args_info->dir_out_orig), &(args_info->dir_out_given),
              &(local_args_info.dir_out_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "dir_out", 'd',
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
