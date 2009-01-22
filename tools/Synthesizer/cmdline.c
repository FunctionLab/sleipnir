/*
  File autogenerated by gengetopt version 2.22
  generated with the following command:
  /r01/tergeo/chuttenh/sleipnir/trunk/../extlib/gengetopt-2.22/bin/gengetopt -iSynthesizer.ggo --default-optional -u -N -e 

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

const char *gengetopt_args_info_purpose = "Creates configurable synthetic microarray and sequence data.";

const char *gengetopt_args_info_usage = "Usage: Synthesizer [OPTIONS]... [FILES]...";

const char *gengetopt_args_info_description = "";

const char *gengetopt_args_info_help[] = {
  "  -h, --help                   Print help and exit",
  "  -V, --version                Print version and exit",
  "\nMain:",
  "  -o, --output_pcl=filename    PCL expression output file",
  "  -O, --output_fasta=filename  FASTA sequence output file",
  "  -g, --genes=INT              Number of synthesized genes  (default=`5000')",
  "  -c, --conditions=INT         Number of synthesized conditions  \n                                 (default=`100')",
  "\nRegulators:",
  "  -n, --tfs=INT                Number of transcription factors  (default=`10')",
  "  -q, --tf_gene=DOUBLE         Probability of TF activity in a gene  \n                                 (default=`0.01')",
  "  -Q, --tf_condition=DOUBLE    Probability of TF activity in a condition  \n                                 (default=`0.1')",
  "  -t, --tf_min=INT             Minimum TFBS length  (default=`5')",
  "  -T, --tf_max=INT             Maximum TFBS length  (default=`12')",
  "\nExpression:",
  "  -m, --mean=DOUBLE            Expression mean  (default=`0')",
  "  -s, --stdev=DOUBLE           Expression standard deviation  (default=`1')",
  "  -M, --tf_mean=DOUBLE         Up/downregulation mean  (default=`2')",
  "  -S, --tf_stdev=DOUBLE        Up/downregulation standard deviation  \n                                 (default=`1')",
  "\nSequence:",
  "  -f, --fasta=filename         Input FASTA file",
  "  -d, --degree=INT             Degree of sequence model HMM  (default=`3')",
  "  -l, --seq_min=INT            Minimum sequence length  (default=`1000')",
  "  -L, --seq_max=INT            Maximum sequence length  (default=`3000')",
  "  -p, --tf_copm=INT            Minimum TFBS copies  (default=`1')",
  "  -P, --tf_copx=INT            Maximum TFBS copies  (default=`5')",
  "  -y, --tf_types=STRING        Sequence types containing TFBSs, comma separated",
  "\nOptional:",
  "  -w, --wrap=INT               Wrap width for FASTA output  (default=`60')",
  "  -r, --random=INT             Seed random generator  (default=`0')",
  "  -v, --verbosity=INT          Message verbosity  (default=`5')",
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
  args_info->output_pcl_given = 0 ;
  args_info->output_fasta_given = 0 ;
  args_info->genes_given = 0 ;
  args_info->conditions_given = 0 ;
  args_info->tfs_given = 0 ;
  args_info->tf_gene_given = 0 ;
  args_info->tf_condition_given = 0 ;
  args_info->tf_min_given = 0 ;
  args_info->tf_max_given = 0 ;
  args_info->mean_given = 0 ;
  args_info->stdev_given = 0 ;
  args_info->tf_mean_given = 0 ;
  args_info->tf_stdev_given = 0 ;
  args_info->fasta_given = 0 ;
  args_info->degree_given = 0 ;
  args_info->seq_min_given = 0 ;
  args_info->seq_max_given = 0 ;
  args_info->tf_copm_given = 0 ;
  args_info->tf_copx_given = 0 ;
  args_info->tf_types_given = 0 ;
  args_info->wrap_given = 0 ;
  args_info->random_given = 0 ;
  args_info->verbosity_given = 0 ;
}

static
void clear_args (struct gengetopt_args_info *args_info)
{
  args_info->output_pcl_arg = NULL;
  args_info->output_pcl_orig = NULL;
  args_info->output_fasta_arg = NULL;
  args_info->output_fasta_orig = NULL;
  args_info->genes_arg = 5000;
  args_info->genes_orig = NULL;
  args_info->conditions_arg = 100;
  args_info->conditions_orig = NULL;
  args_info->tfs_arg = 10;
  args_info->tfs_orig = NULL;
  args_info->tf_gene_arg = 0.01;
  args_info->tf_gene_orig = NULL;
  args_info->tf_condition_arg = 0.1;
  args_info->tf_condition_orig = NULL;
  args_info->tf_min_arg = 5;
  args_info->tf_min_orig = NULL;
  args_info->tf_max_arg = 12;
  args_info->tf_max_orig = NULL;
  args_info->mean_arg = 0;
  args_info->mean_orig = NULL;
  args_info->stdev_arg = 1;
  args_info->stdev_orig = NULL;
  args_info->tf_mean_arg = 2;
  args_info->tf_mean_orig = NULL;
  args_info->tf_stdev_arg = 1;
  args_info->tf_stdev_orig = NULL;
  args_info->fasta_arg = NULL;
  args_info->fasta_orig = NULL;
  args_info->degree_arg = 3;
  args_info->degree_orig = NULL;
  args_info->seq_min_arg = 1000;
  args_info->seq_min_orig = NULL;
  args_info->seq_max_arg = 3000;
  args_info->seq_max_orig = NULL;
  args_info->tf_copm_arg = 1;
  args_info->tf_copm_orig = NULL;
  args_info->tf_copx_arg = 5;
  args_info->tf_copx_orig = NULL;
  args_info->tf_types_arg = NULL;
  args_info->tf_types_orig = NULL;
  args_info->wrap_arg = 60;
  args_info->wrap_orig = NULL;
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
  args_info->output_pcl_help = gengetopt_args_info_help[3] ;
  args_info->output_fasta_help = gengetopt_args_info_help[4] ;
  args_info->genes_help = gengetopt_args_info_help[5] ;
  args_info->conditions_help = gengetopt_args_info_help[6] ;
  args_info->tfs_help = gengetopt_args_info_help[8] ;
  args_info->tf_gene_help = gengetopt_args_info_help[9] ;
  args_info->tf_condition_help = gengetopt_args_info_help[10] ;
  args_info->tf_min_help = gengetopt_args_info_help[11] ;
  args_info->tf_max_help = gengetopt_args_info_help[12] ;
  args_info->mean_help = gengetopt_args_info_help[14] ;
  args_info->stdev_help = gengetopt_args_info_help[15] ;
  args_info->tf_mean_help = gengetopt_args_info_help[16] ;
  args_info->tf_stdev_help = gengetopt_args_info_help[17] ;
  args_info->fasta_help = gengetopt_args_info_help[19] ;
  args_info->degree_help = gengetopt_args_info_help[20] ;
  args_info->seq_min_help = gengetopt_args_info_help[21] ;
  args_info->seq_max_help = gengetopt_args_info_help[22] ;
  args_info->tf_copm_help = gengetopt_args_info_help[23] ;
  args_info->tf_copx_help = gengetopt_args_info_help[24] ;
  args_info->tf_types_help = gengetopt_args_info_help[25] ;
  args_info->wrap_help = gengetopt_args_info_help[27] ;
  args_info->random_help = gengetopt_args_info_help[28] ;
  args_info->verbosity_help = gengetopt_args_info_help[29] ;
  
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
  free_string_field (&(args_info->output_pcl_arg));
  free_string_field (&(args_info->output_pcl_orig));
  free_string_field (&(args_info->output_fasta_arg));
  free_string_field (&(args_info->output_fasta_orig));
  free_string_field (&(args_info->genes_orig));
  free_string_field (&(args_info->conditions_orig));
  free_string_field (&(args_info->tfs_orig));
  free_string_field (&(args_info->tf_gene_orig));
  free_string_field (&(args_info->tf_condition_orig));
  free_string_field (&(args_info->tf_min_orig));
  free_string_field (&(args_info->tf_max_orig));
  free_string_field (&(args_info->mean_orig));
  free_string_field (&(args_info->stdev_orig));
  free_string_field (&(args_info->tf_mean_orig));
  free_string_field (&(args_info->tf_stdev_orig));
  free_string_field (&(args_info->fasta_arg));
  free_string_field (&(args_info->fasta_orig));
  free_string_field (&(args_info->degree_orig));
  free_string_field (&(args_info->seq_min_orig));
  free_string_field (&(args_info->seq_max_orig));
  free_string_field (&(args_info->tf_copm_orig));
  free_string_field (&(args_info->tf_copx_orig));
  free_string_field (&(args_info->tf_types_arg));
  free_string_field (&(args_info->tf_types_orig));
  free_string_field (&(args_info->wrap_orig));
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
  if (args_info->output_pcl_given)
    write_into_file(outfile, "output_pcl", args_info->output_pcl_orig, 0);
  if (args_info->output_fasta_given)
    write_into_file(outfile, "output_fasta", args_info->output_fasta_orig, 0);
  if (args_info->genes_given)
    write_into_file(outfile, "genes", args_info->genes_orig, 0);
  if (args_info->conditions_given)
    write_into_file(outfile, "conditions", args_info->conditions_orig, 0);
  if (args_info->tfs_given)
    write_into_file(outfile, "tfs", args_info->tfs_orig, 0);
  if (args_info->tf_gene_given)
    write_into_file(outfile, "tf_gene", args_info->tf_gene_orig, 0);
  if (args_info->tf_condition_given)
    write_into_file(outfile, "tf_condition", args_info->tf_condition_orig, 0);
  if (args_info->tf_min_given)
    write_into_file(outfile, "tf_min", args_info->tf_min_orig, 0);
  if (args_info->tf_max_given)
    write_into_file(outfile, "tf_max", args_info->tf_max_orig, 0);
  if (args_info->mean_given)
    write_into_file(outfile, "mean", args_info->mean_orig, 0);
  if (args_info->stdev_given)
    write_into_file(outfile, "stdev", args_info->stdev_orig, 0);
  if (args_info->tf_mean_given)
    write_into_file(outfile, "tf_mean", args_info->tf_mean_orig, 0);
  if (args_info->tf_stdev_given)
    write_into_file(outfile, "tf_stdev", args_info->tf_stdev_orig, 0);
  if (args_info->fasta_given)
    write_into_file(outfile, "fasta", args_info->fasta_orig, 0);
  if (args_info->degree_given)
    write_into_file(outfile, "degree", args_info->degree_orig, 0);
  if (args_info->seq_min_given)
    write_into_file(outfile, "seq_min", args_info->seq_min_orig, 0);
  if (args_info->seq_max_given)
    write_into_file(outfile, "seq_max", args_info->seq_max_orig, 0);
  if (args_info->tf_copm_given)
    write_into_file(outfile, "tf_copm", args_info->tf_copm_orig, 0);
  if (args_info->tf_copx_given)
    write_into_file(outfile, "tf_copx", args_info->tf_copx_orig, 0);
  if (args_info->tf_types_given)
    write_into_file(outfile, "tf_types", args_info->tf_types_orig, 0);
  if (args_info->wrap_given)
    write_into_file(outfile, "wrap", args_info->wrap_orig, 0);
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
        { "output_pcl",	1, NULL, 'o' },
        { "output_fasta",	1, NULL, 'O' },
        { "genes",	1, NULL, 'g' },
        { "conditions",	1, NULL, 'c' },
        { "tfs",	1, NULL, 'n' },
        { "tf_gene",	1, NULL, 'q' },
        { "tf_condition",	1, NULL, 'Q' },
        { "tf_min",	1, NULL, 't' },
        { "tf_max",	1, NULL, 'T' },
        { "mean",	1, NULL, 'm' },
        { "stdev",	1, NULL, 's' },
        { "tf_mean",	1, NULL, 'M' },
        { "tf_stdev",	1, NULL, 'S' },
        { "fasta",	1, NULL, 'f' },
        { "degree",	1, NULL, 'd' },
        { "seq_min",	1, NULL, 'l' },
        { "seq_max",	1, NULL, 'L' },
        { "tf_copm",	1, NULL, 'p' },
        { "tf_copx",	1, NULL, 'P' },
        { "tf_types",	1, NULL, 'y' },
        { "wrap",	1, NULL, 'w' },
        { "random",	1, NULL, 'r' },
        { "verbosity",	1, NULL, 'v' },
        { NULL,	0, NULL, 0 }
      };

      c = getopt_long (argc, argv, "hVo:O:g:c:n:q:Q:t:T:m:s:M:S:f:d:l:L:p:P:y:w:r:v:", long_options, &option_index);

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
        case 'o':	/* PCL expression output file.  */
        
        
          if (update_arg( (void *)&(args_info->output_pcl_arg), 
               &(args_info->output_pcl_orig), &(args_info->output_pcl_given),
              &(local_args_info.output_pcl_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "output_pcl", 'o',
              additional_error))
            goto failure;
        
          break;
        case 'O':	/* FASTA sequence output file.  */
        
        
          if (update_arg( (void *)&(args_info->output_fasta_arg), 
               &(args_info->output_fasta_orig), &(args_info->output_fasta_given),
              &(local_args_info.output_fasta_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "output_fasta", 'O',
              additional_error))
            goto failure;
        
          break;
        case 'g':	/* Number of synthesized genes.  */
        
        
          if (update_arg( (void *)&(args_info->genes_arg), 
               &(args_info->genes_orig), &(args_info->genes_given),
              &(local_args_info.genes_given), optarg, 0, "5000", ARG_INT,
              check_ambiguity, override, 0, 0,
              "genes", 'g',
              additional_error))
            goto failure;
        
          break;
        case 'c':	/* Number of synthesized conditions.  */
        
        
          if (update_arg( (void *)&(args_info->conditions_arg), 
               &(args_info->conditions_orig), &(args_info->conditions_given),
              &(local_args_info.conditions_given), optarg, 0, "100", ARG_INT,
              check_ambiguity, override, 0, 0,
              "conditions", 'c',
              additional_error))
            goto failure;
        
          break;
        case 'n':	/* Number of transcription factors.  */
        
        
          if (update_arg( (void *)&(args_info->tfs_arg), 
               &(args_info->tfs_orig), &(args_info->tfs_given),
              &(local_args_info.tfs_given), optarg, 0, "10", ARG_INT,
              check_ambiguity, override, 0, 0,
              "tfs", 'n',
              additional_error))
            goto failure;
        
          break;
        case 'q':	/* Probability of TF activity in a gene.  */
        
        
          if (update_arg( (void *)&(args_info->tf_gene_arg), 
               &(args_info->tf_gene_orig), &(args_info->tf_gene_given),
              &(local_args_info.tf_gene_given), optarg, 0, "0.01", ARG_DOUBLE,
              check_ambiguity, override, 0, 0,
              "tf_gene", 'q',
              additional_error))
            goto failure;
        
          break;
        case 'Q':	/* Probability of TF activity in a condition.  */
        
        
          if (update_arg( (void *)&(args_info->tf_condition_arg), 
               &(args_info->tf_condition_orig), &(args_info->tf_condition_given),
              &(local_args_info.tf_condition_given), optarg, 0, "0.1", ARG_DOUBLE,
              check_ambiguity, override, 0, 0,
              "tf_condition", 'Q',
              additional_error))
            goto failure;
        
          break;
        case 't':	/* Minimum TFBS length.  */
        
        
          if (update_arg( (void *)&(args_info->tf_min_arg), 
               &(args_info->tf_min_orig), &(args_info->tf_min_given),
              &(local_args_info.tf_min_given), optarg, 0, "5", ARG_INT,
              check_ambiguity, override, 0, 0,
              "tf_min", 't',
              additional_error))
            goto failure;
        
          break;
        case 'T':	/* Maximum TFBS length.  */
        
        
          if (update_arg( (void *)&(args_info->tf_max_arg), 
               &(args_info->tf_max_orig), &(args_info->tf_max_given),
              &(local_args_info.tf_max_given), optarg, 0, "12", ARG_INT,
              check_ambiguity, override, 0, 0,
              "tf_max", 'T',
              additional_error))
            goto failure;
        
          break;
        case 'm':	/* Expression mean.  */
        
        
          if (update_arg( (void *)&(args_info->mean_arg), 
               &(args_info->mean_orig), &(args_info->mean_given),
              &(local_args_info.mean_given), optarg, 0, "0", ARG_DOUBLE,
              check_ambiguity, override, 0, 0,
              "mean", 'm',
              additional_error))
            goto failure;
        
          break;
        case 's':	/* Expression standard deviation.  */
        
        
          if (update_arg( (void *)&(args_info->stdev_arg), 
               &(args_info->stdev_orig), &(args_info->stdev_given),
              &(local_args_info.stdev_given), optarg, 0, "1", ARG_DOUBLE,
              check_ambiguity, override, 0, 0,
              "stdev", 's',
              additional_error))
            goto failure;
        
          break;
        case 'M':	/* Up/downregulation mean.  */
        
        
          if (update_arg( (void *)&(args_info->tf_mean_arg), 
               &(args_info->tf_mean_orig), &(args_info->tf_mean_given),
              &(local_args_info.tf_mean_given), optarg, 0, "2", ARG_DOUBLE,
              check_ambiguity, override, 0, 0,
              "tf_mean", 'M',
              additional_error))
            goto failure;
        
          break;
        case 'S':	/* Up/downregulation standard deviation.  */
        
        
          if (update_arg( (void *)&(args_info->tf_stdev_arg), 
               &(args_info->tf_stdev_orig), &(args_info->tf_stdev_given),
              &(local_args_info.tf_stdev_given), optarg, 0, "1", ARG_DOUBLE,
              check_ambiguity, override, 0, 0,
              "tf_stdev", 'S',
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
        case 'd':	/* Degree of sequence model HMM.  */
        
        
          if (update_arg( (void *)&(args_info->degree_arg), 
               &(args_info->degree_orig), &(args_info->degree_given),
              &(local_args_info.degree_given), optarg, 0, "3", ARG_INT,
              check_ambiguity, override, 0, 0,
              "degree", 'd',
              additional_error))
            goto failure;
        
          break;
        case 'l':	/* Minimum sequence length.  */
        
        
          if (update_arg( (void *)&(args_info->seq_min_arg), 
               &(args_info->seq_min_orig), &(args_info->seq_min_given),
              &(local_args_info.seq_min_given), optarg, 0, "1000", ARG_INT,
              check_ambiguity, override, 0, 0,
              "seq_min", 'l',
              additional_error))
            goto failure;
        
          break;
        case 'L':	/* Maximum sequence length.  */
        
        
          if (update_arg( (void *)&(args_info->seq_max_arg), 
               &(args_info->seq_max_orig), &(args_info->seq_max_given),
              &(local_args_info.seq_max_given), optarg, 0, "3000", ARG_INT,
              check_ambiguity, override, 0, 0,
              "seq_max", 'L',
              additional_error))
            goto failure;
        
          break;
        case 'p':	/* Minimum TFBS copies.  */
        
        
          if (update_arg( (void *)&(args_info->tf_copm_arg), 
               &(args_info->tf_copm_orig), &(args_info->tf_copm_given),
              &(local_args_info.tf_copm_given), optarg, 0, "1", ARG_INT,
              check_ambiguity, override, 0, 0,
              "tf_copm", 'p',
              additional_error))
            goto failure;
        
          break;
        case 'P':	/* Maximum TFBS copies.  */
        
        
          if (update_arg( (void *)&(args_info->tf_copx_arg), 
               &(args_info->tf_copx_orig), &(args_info->tf_copx_given),
              &(local_args_info.tf_copx_given), optarg, 0, "5", ARG_INT,
              check_ambiguity, override, 0, 0,
              "tf_copx", 'P',
              additional_error))
            goto failure;
        
          break;
        case 'y':	/* Sequence types containing TFBSs, comma separated.  */
        
        
          if (update_arg( (void *)&(args_info->tf_types_arg), 
               &(args_info->tf_types_orig), &(args_info->tf_types_given),
              &(local_args_info.tf_types_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "tf_types", 'y',
              additional_error))
            goto failure;
        
          break;
        case 'w':	/* Wrap width for FASTA output.  */
        
        
          if (update_arg( (void *)&(args_info->wrap_arg), 
               &(args_info->wrap_orig), &(args_info->wrap_given),
              &(local_args_info.wrap_given), optarg, 0, "60", ARG_INT,
              check_ambiguity, override, 0, 0,
              "wrap", 'w',
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
