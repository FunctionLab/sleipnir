/*
  File autogenerated by gengetopt version 2.22
  generated with the following command:
  ..\..\..\..\extlib\proj\vs2005\release\gengetopt.exe --default-optional -N -e --output-dir=c:\Documents and Settings\eblis\My Documents\Visual Studio 2008\Projects\Sleipnir\trunk\tools\COALESCE\ -i ..\..\..\tools\COALESCE\COALESCE.ggo 

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

const char *gengetopt_args_info_purpose = "Implements a Bayesian Iterative Signature Algorithm for biclustering and TFBS \ndiscovery.";

const char *gengetopt_args_info_usage = "Usage: COALESCE [OPTIONS]...";

const char *gengetopt_args_info_description = "";

const char *gengetopt_args_info_help[] = {
  "  -h, --help                    Print help and exit",
  "  -V, --version                 Print version and exit",
  "\nMain:",
  "  -i, --input=filename          Input PCL file",
  "  -f, --fasta=filename          Input FASTA file",
  "  -d, --datasets=filename       Condition groupings into dataset blocks",
  "  -o, --output=directory        Directory for intermediate output files (PCLs)",
  "\nAlgorithm Parameters:",
  "  -p, --prob_gene=DOUBLE        Probability threshhold for gene inclusion  \n                                  (default=`0.95')",
  "  -P, --pvalue_cond=DOUBLE      P-value threshhold for condition inclusion  \n                                  (default=`0.05')",
  "  -m, --pvalue_motif=DOUBLE     P-value threshhold for motif inclusion  \n                                  (default=`0.05')",
  "\nSequence Parameters:",
  "  -k, --k=INT                   Sequence kmer length  (default=`7')",
  "  -g, --pvalue_merge=DOUBLE     P-value threshhold for motif merging  \n                                  (default=`0.05')",
  "  -G, --cutoff_merge=DOUBLE     Edit distance cutoff for motif merging  \n                                  (default=`2.5')",
  "  -y, --penalty_gap=DOUBLE      Edit distance penalty for gaps  (default=`1')",
  "  -Y, --penalty_mismatch=DOUBLE Edit distance penalty for mismatches  \n                                  (default=`2.1')",
  "\nPerformance Parameters:",
  "  -c, --pvalue_correl=DOUBLE    P-value threshhold for significant correlation  \n                                  (default=`0.05')",
  "  -C, --number_correl=INT       Maximum number of pairs to sample for \n                                  significant correlation  (default=`100000')",
  "  -q, --sequences=STRING        Sequence types to use (comma separated)",
  "  -b, --bases=INT               Resolution of bases per motif match  \n                                  (default=`5000')",
  "  -z, --size_minimum=INT        Minimum gene count for clusters of interest  \n                                  (default=`5')",
  "  -Z, --size_maximum=INT        Maximum motif count to consider a cluster \n                                  saturated  (default=`100')",
  "\nAdditional Data:",
  "  -n, --nucleosomes=filename    Nucleosome position file (ENCODE/PCL)",
  "\nMiscellaneous:",
  "  -e, --cache=filename          Cache file for sequence analysis",
  "\nOptional:",
  "  -t, --threads=INT             Maximum number of concurrent threads  \n                                  (default=`1')",
  "  -s, --skip=INT                Columns to skip in input PCL  (default=`2')",
  "  -r, --random=INT              Seed random generator  (default=`0')",
  "  -v, --verbosity=INT           Message verbosity  (default=`5')",
    0
};

typedef enum {ARG_NO
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


static char *
gengetopt_strdup (const char *s);

static
void clear_given (struct gengetopt_args_info *args_info)
{
  args_info->help_given = 0 ;
  args_info->version_given = 0 ;
  args_info->input_given = 0 ;
  args_info->fasta_given = 0 ;
  args_info->datasets_given = 0 ;
  args_info->output_given = 0 ;
  args_info->prob_gene_given = 0 ;
  args_info->pvalue_cond_given = 0 ;
  args_info->pvalue_motif_given = 0 ;
  args_info->k_given = 0 ;
  args_info->pvalue_merge_given = 0 ;
  args_info->cutoff_merge_given = 0 ;
  args_info->penalty_gap_given = 0 ;
  args_info->penalty_mismatch_given = 0 ;
  args_info->pvalue_correl_given = 0 ;
  args_info->number_correl_given = 0 ;
  args_info->sequences_given = 0 ;
  args_info->bases_given = 0 ;
  args_info->size_minimum_given = 0 ;
  args_info->size_maximum_given = 0 ;
  args_info->nucleosomes_given = 0 ;
  args_info->cache_given = 0 ;
  args_info->threads_given = 0 ;
  args_info->skip_given = 0 ;
  args_info->random_given = 0 ;
  args_info->verbosity_given = 0 ;
}

static
void clear_args (struct gengetopt_args_info *args_info)
{
  args_info->input_arg = NULL;
  args_info->input_orig = NULL;
  args_info->fasta_arg = NULL;
  args_info->fasta_orig = NULL;
  args_info->datasets_arg = NULL;
  args_info->datasets_orig = NULL;
  args_info->output_arg = NULL;
  args_info->output_orig = NULL;
  args_info->prob_gene_arg = 0.95;
  args_info->prob_gene_orig = NULL;
  args_info->pvalue_cond_arg = 0.05;
  args_info->pvalue_cond_orig = NULL;
  args_info->pvalue_motif_arg = 0.05;
  args_info->pvalue_motif_orig = NULL;
  args_info->k_arg = 7;
  args_info->k_orig = NULL;
  args_info->pvalue_merge_arg = 0.05;
  args_info->pvalue_merge_orig = NULL;
  args_info->cutoff_merge_arg = 2.5;
  args_info->cutoff_merge_orig = NULL;
  args_info->penalty_gap_arg = 1;
  args_info->penalty_gap_orig = NULL;
  args_info->penalty_mismatch_arg = 2.1;
  args_info->penalty_mismatch_orig = NULL;
  args_info->pvalue_correl_arg = 0.05;
  args_info->pvalue_correl_orig = NULL;
  args_info->number_correl_arg = 100000;
  args_info->number_correl_orig = NULL;
  args_info->sequences_arg = NULL;
  args_info->sequences_orig = NULL;
  args_info->bases_arg = 5000;
  args_info->bases_orig = NULL;
  args_info->size_minimum_arg = 5;
  args_info->size_minimum_orig = NULL;
  args_info->size_maximum_arg = 100;
  args_info->size_maximum_orig = NULL;
  args_info->nucleosomes_arg = NULL;
  args_info->nucleosomes_orig = NULL;
  args_info->cache_arg = NULL;
  args_info->cache_orig = NULL;
  args_info->threads_arg = 1;
  args_info->threads_orig = NULL;
  args_info->skip_arg = 2;
  args_info->skip_orig = NULL;
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
  args_info->input_help = gengetopt_args_info_help[3] ;
  args_info->fasta_help = gengetopt_args_info_help[4] ;
  args_info->datasets_help = gengetopt_args_info_help[5] ;
  args_info->output_help = gengetopt_args_info_help[6] ;
  args_info->prob_gene_help = gengetopt_args_info_help[8] ;
  args_info->pvalue_cond_help = gengetopt_args_info_help[9] ;
  args_info->pvalue_motif_help = gengetopt_args_info_help[10] ;
  args_info->k_help = gengetopt_args_info_help[12] ;
  args_info->pvalue_merge_help = gengetopt_args_info_help[13] ;
  args_info->cutoff_merge_help = gengetopt_args_info_help[14] ;
  args_info->penalty_gap_help = gengetopt_args_info_help[15] ;
  args_info->penalty_mismatch_help = gengetopt_args_info_help[16] ;
  args_info->pvalue_correl_help = gengetopt_args_info_help[18] ;
  args_info->number_correl_help = gengetopt_args_info_help[19] ;
  args_info->sequences_help = gengetopt_args_info_help[20] ;
  args_info->bases_help = gengetopt_args_info_help[21] ;
  args_info->size_minimum_help = gengetopt_args_info_help[22] ;
  args_info->size_maximum_help = gengetopt_args_info_help[23] ;
  args_info->nucleosomes_help = gengetopt_args_info_help[25] ;
  args_info->cache_help = gengetopt_args_info_help[27] ;
  args_info->threads_help = gengetopt_args_info_help[29] ;
  args_info->skip_help = gengetopt_args_info_help[30] ;
  args_info->random_help = gengetopt_args_info_help[31] ;
  args_info->verbosity_help = gengetopt_args_info_help[32] ;
  
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

  free_string_field (&(args_info->input_arg));
  free_string_field (&(args_info->input_orig));
  free_string_field (&(args_info->fasta_arg));
  free_string_field (&(args_info->fasta_orig));
  free_string_field (&(args_info->datasets_arg));
  free_string_field (&(args_info->datasets_orig));
  free_string_field (&(args_info->output_arg));
  free_string_field (&(args_info->output_orig));
  free_string_field (&(args_info->prob_gene_orig));
  free_string_field (&(args_info->pvalue_cond_orig));
  free_string_field (&(args_info->pvalue_motif_orig));
  free_string_field (&(args_info->k_orig));
  free_string_field (&(args_info->pvalue_merge_orig));
  free_string_field (&(args_info->cutoff_merge_orig));
  free_string_field (&(args_info->penalty_gap_orig));
  free_string_field (&(args_info->penalty_mismatch_orig));
  free_string_field (&(args_info->pvalue_correl_orig));
  free_string_field (&(args_info->number_correl_orig));
  free_string_field (&(args_info->sequences_arg));
  free_string_field (&(args_info->sequences_orig));
  free_string_field (&(args_info->bases_orig));
  free_string_field (&(args_info->size_minimum_orig));
  free_string_field (&(args_info->size_maximum_orig));
  free_string_field (&(args_info->nucleosomes_arg));
  free_string_field (&(args_info->nucleosomes_orig));
  free_string_field (&(args_info->cache_arg));
  free_string_field (&(args_info->cache_orig));
  free_string_field (&(args_info->threads_orig));
  free_string_field (&(args_info->skip_orig));
  free_string_field (&(args_info->random_orig));
  free_string_field (&(args_info->verbosity_orig));
  
  

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
  if (args_info->input_given)
    write_into_file(outfile, "input", args_info->input_orig, 0);
  if (args_info->fasta_given)
    write_into_file(outfile, "fasta", args_info->fasta_orig, 0);
  if (args_info->datasets_given)
    write_into_file(outfile, "datasets", args_info->datasets_orig, 0);
  if (args_info->output_given)
    write_into_file(outfile, "output", args_info->output_orig, 0);
  if (args_info->prob_gene_given)
    write_into_file(outfile, "prob_gene", args_info->prob_gene_orig, 0);
  if (args_info->pvalue_cond_given)
    write_into_file(outfile, "pvalue_cond", args_info->pvalue_cond_orig, 0);
  if (args_info->pvalue_motif_given)
    write_into_file(outfile, "pvalue_motif", args_info->pvalue_motif_orig, 0);
  if (args_info->k_given)
    write_into_file(outfile, "k", args_info->k_orig, 0);
  if (args_info->pvalue_merge_given)
    write_into_file(outfile, "pvalue_merge", args_info->pvalue_merge_orig, 0);
  if (args_info->cutoff_merge_given)
    write_into_file(outfile, "cutoff_merge", args_info->cutoff_merge_orig, 0);
  if (args_info->penalty_gap_given)
    write_into_file(outfile, "penalty_gap", args_info->penalty_gap_orig, 0);
  if (args_info->penalty_mismatch_given)
    write_into_file(outfile, "penalty_mismatch", args_info->penalty_mismatch_orig, 0);
  if (args_info->pvalue_correl_given)
    write_into_file(outfile, "pvalue_correl", args_info->pvalue_correl_orig, 0);
  if (args_info->number_correl_given)
    write_into_file(outfile, "number_correl", args_info->number_correl_orig, 0);
  if (args_info->sequences_given)
    write_into_file(outfile, "sequences", args_info->sequences_orig, 0);
  if (args_info->bases_given)
    write_into_file(outfile, "bases", args_info->bases_orig, 0);
  if (args_info->size_minimum_given)
    write_into_file(outfile, "size_minimum", args_info->size_minimum_orig, 0);
  if (args_info->size_maximum_given)
    write_into_file(outfile, "size_maximum", args_info->size_maximum_orig, 0);
  if (args_info->nucleosomes_given)
    write_into_file(outfile, "nucleosomes", args_info->nucleosomes_orig, 0);
  if (args_info->cache_given)
    write_into_file(outfile, "cache", args_info->cache_orig, 0);
  if (args_info->threads_given)
    write_into_file(outfile, "threads", args_info->threads_orig, 0);
  if (args_info->skip_given)
    write_into_file(outfile, "skip", args_info->skip_orig, 0);
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
        { "input",	1, NULL, 'i' },
        { "fasta",	1, NULL, 'f' },
        { "datasets",	1, NULL, 'd' },
        { "output",	1, NULL, 'o' },
        { "prob_gene",	1, NULL, 'p' },
        { "pvalue_cond",	1, NULL, 'P' },
        { "pvalue_motif",	1, NULL, 'm' },
        { "k",	1, NULL, 'k' },
        { "pvalue_merge",	1, NULL, 'g' },
        { "cutoff_merge",	1, NULL, 'G' },
        { "penalty_gap",	1, NULL, 'y' },
        { "penalty_mismatch",	1, NULL, 'Y' },
        { "pvalue_correl",	1, NULL, 'c' },
        { "number_correl",	1, NULL, 'C' },
        { "sequences",	1, NULL, 'q' },
        { "bases",	1, NULL, 'b' },
        { "size_minimum",	1, NULL, 'z' },
        { "size_maximum",	1, NULL, 'Z' },
        { "nucleosomes",	1, NULL, 'n' },
        { "cache",	1, NULL, 'e' },
        { "threads",	1, NULL, 't' },
        { "skip",	1, NULL, 's' },
        { "random",	1, NULL, 'r' },
        { "verbosity",	1, NULL, 'v' },
        { NULL,	0, NULL, 0 }
      };

      c = getopt_long (argc, argv, "hVi:f:d:o:p:P:m:k:g:G:y:Y:c:C:q:b:z:Z:n:e:t:s:r:v:", long_options, &option_index);

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
        case 'i':	/* Input PCL file.  */
        
        
          if (update_arg( (void *)&(args_info->input_arg), 
               &(args_info->input_orig), &(args_info->input_given),
              &(local_args_info.input_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "input", 'i',
              additional_error))
            goto failure;
        
          break;
        case 'f':	/* Input FASTA file.  */
        
        
          if (update_arg( (void *)&(args_info->fasta_arg), 
               &(args_info->fasta_orig), &(args_info->fasta_given),
              &(local_args_info.fasta_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "fasta", 'f',
              additional_error))
            goto failure;
        
          break;
        case 'd':	/* Condition groupings into dataset blocks.  */
        
        
          if (update_arg( (void *)&(args_info->datasets_arg), 
               &(args_info->datasets_orig), &(args_info->datasets_given),
              &(local_args_info.datasets_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "datasets", 'd',
              additional_error))
            goto failure;
        
          break;
        case 'o':	/* Directory for intermediate output files (PCLs).  */
        
        
          if (update_arg( (void *)&(args_info->output_arg), 
               &(args_info->output_orig), &(args_info->output_given),
              &(local_args_info.output_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "output", 'o',
              additional_error))
            goto failure;
        
          break;
        case 'p':	/* Probability threshhold for gene inclusion.  */
        
        
          if (update_arg( (void *)&(args_info->prob_gene_arg), 
               &(args_info->prob_gene_orig), &(args_info->prob_gene_given),
              &(local_args_info.prob_gene_given), optarg, 0, "0.95", ARG_DOUBLE,
              check_ambiguity, override, 0, 0,
              "prob_gene", 'p',
              additional_error))
            goto failure;
        
          break;
        case 'P':	/* P-value threshhold for condition inclusion.  */
        
        
          if (update_arg( (void *)&(args_info->pvalue_cond_arg), 
               &(args_info->pvalue_cond_orig), &(args_info->pvalue_cond_given),
              &(local_args_info.pvalue_cond_given), optarg, 0, "0.05", ARG_DOUBLE,
              check_ambiguity, override, 0, 0,
              "pvalue_cond", 'P',
              additional_error))
            goto failure;
        
          break;
        case 'm':	/* P-value threshhold for motif inclusion.  */
        
        
          if (update_arg( (void *)&(args_info->pvalue_motif_arg), 
               &(args_info->pvalue_motif_orig), &(args_info->pvalue_motif_given),
              &(local_args_info.pvalue_motif_given), optarg, 0, "0.05", ARG_DOUBLE,
              check_ambiguity, override, 0, 0,
              "pvalue_motif", 'm',
              additional_error))
            goto failure;
        
          break;
        case 'k':	/* Sequence kmer length.  */
        
        
          if (update_arg( (void *)&(args_info->k_arg), 
               &(args_info->k_orig), &(args_info->k_given),
              &(local_args_info.k_given), optarg, 0, "7", ARG_INT,
              check_ambiguity, override, 0, 0,
              "k", 'k',
              additional_error))
            goto failure;
        
          break;
        case 'g':	/* P-value threshhold for motif merging.  */
        
        
          if (update_arg( (void *)&(args_info->pvalue_merge_arg), 
               &(args_info->pvalue_merge_orig), &(args_info->pvalue_merge_given),
              &(local_args_info.pvalue_merge_given), optarg, 0, "0.05", ARG_DOUBLE,
              check_ambiguity, override, 0, 0,
              "pvalue_merge", 'g',
              additional_error))
            goto failure;
        
          break;
        case 'G':	/* Edit distance cutoff for motif merging.  */
        
        
          if (update_arg( (void *)&(args_info->cutoff_merge_arg), 
               &(args_info->cutoff_merge_orig), &(args_info->cutoff_merge_given),
              &(local_args_info.cutoff_merge_given), optarg, 0, "2.5", ARG_DOUBLE,
              check_ambiguity, override, 0, 0,
              "cutoff_merge", 'G',
              additional_error))
            goto failure;
        
          break;
        case 'y':	/* Edit distance penalty for gaps.  */
        
        
          if (update_arg( (void *)&(args_info->penalty_gap_arg), 
               &(args_info->penalty_gap_orig), &(args_info->penalty_gap_given),
              &(local_args_info.penalty_gap_given), optarg, 0, "1", ARG_DOUBLE,
              check_ambiguity, override, 0, 0,
              "penalty_gap", 'y',
              additional_error))
            goto failure;
        
          break;
        case 'Y':	/* Edit distance penalty for mismatches.  */
        
        
          if (update_arg( (void *)&(args_info->penalty_mismatch_arg), 
               &(args_info->penalty_mismatch_orig), &(args_info->penalty_mismatch_given),
              &(local_args_info.penalty_mismatch_given), optarg, 0, "2.1", ARG_DOUBLE,
              check_ambiguity, override, 0, 0,
              "penalty_mismatch", 'Y',
              additional_error))
            goto failure;
        
          break;
        case 'c':	/* P-value threshhold for significant correlation.  */
        
        
          if (update_arg( (void *)&(args_info->pvalue_correl_arg), 
               &(args_info->pvalue_correl_orig), &(args_info->pvalue_correl_given),
              &(local_args_info.pvalue_correl_given), optarg, 0, "0.05", ARG_DOUBLE,
              check_ambiguity, override, 0, 0,
              "pvalue_correl", 'c',
              additional_error))
            goto failure;
        
          break;
        case 'C':	/* Maximum number of pairs to sample for significant correlation.  */
        
        
          if (update_arg( (void *)&(args_info->number_correl_arg), 
               &(args_info->number_correl_orig), &(args_info->number_correl_given),
              &(local_args_info.number_correl_given), optarg, 0, "100000", ARG_INT,
              check_ambiguity, override, 0, 0,
              "number_correl", 'C',
              additional_error))
            goto failure;
        
          break;
        case 'q':	/* Sequence types to use (comma separated).  */
        
        
          if (update_arg( (void *)&(args_info->sequences_arg), 
               &(args_info->sequences_orig), &(args_info->sequences_given),
              &(local_args_info.sequences_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "sequences", 'q',
              additional_error))
            goto failure;
        
          break;
        case 'b':	/* Resolution of bases per motif match.  */
        
        
          if (update_arg( (void *)&(args_info->bases_arg), 
               &(args_info->bases_orig), &(args_info->bases_given),
              &(local_args_info.bases_given), optarg, 0, "5000", ARG_INT,
              check_ambiguity, override, 0, 0,
              "bases", 'b',
              additional_error))
            goto failure;
        
          break;
        case 'z':	/* Minimum gene count for clusters of interest.  */
        
        
          if (update_arg( (void *)&(args_info->size_minimum_arg), 
               &(args_info->size_minimum_orig), &(args_info->size_minimum_given),
              &(local_args_info.size_minimum_given), optarg, 0, "5", ARG_INT,
              check_ambiguity, override, 0, 0,
              "size_minimum", 'z',
              additional_error))
            goto failure;
        
          break;
        case 'Z':	/* Maximum motif count to consider a cluster saturated.  */
        
        
          if (update_arg( (void *)&(args_info->size_maximum_arg), 
               &(args_info->size_maximum_orig), &(args_info->size_maximum_given),
              &(local_args_info.size_maximum_given), optarg, 0, "100", ARG_INT,
              check_ambiguity, override, 0, 0,
              "size_maximum", 'Z',
              additional_error))
            goto failure;
        
          break;
        case 'n':	/* Nucleosome position file (ENCODE/PCL).  */
        
        
          if (update_arg( (void *)&(args_info->nucleosomes_arg), 
               &(args_info->nucleosomes_orig), &(args_info->nucleosomes_given),
              &(local_args_info.nucleosomes_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "nucleosomes", 'n',
              additional_error))
            goto failure;
        
          break;
        case 'e':	/* Cache file for sequence analysis.  */
        
        
          if (update_arg( (void *)&(args_info->cache_arg), 
               &(args_info->cache_orig), &(args_info->cache_given),
              &(local_args_info.cache_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "cache", 'e',
              additional_error))
            goto failure;
        
          break;
        case 't':	/* Maximum number of concurrent threads.  */
        
        
          if (update_arg( (void *)&(args_info->threads_arg), 
               &(args_info->threads_orig), &(args_info->threads_given),
              &(local_args_info.threads_given), optarg, 0, "1", ARG_INT,
              check_ambiguity, override, 0, 0,
              "threads", 't',
              additional_error))
            goto failure;
        
          break;
        case 's':	/* Columns to skip in input PCL.  */
        
        
          if (update_arg( (void *)&(args_info->skip_arg), 
               &(args_info->skip_orig), &(args_info->skip_given),
              &(local_args_info.skip_given), optarg, 0, "2", ARG_INT,
              check_ambiguity, override, 0, 0,
              "skip", 's',
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




  cmdline_parser_release (&local_args_info);

  if ( error )
    return (EXIT_FAILURE);

  return 0;

failure:
  
  cmdline_parser_release (&local_args_info);
  return (EXIT_FAILURE);
}
