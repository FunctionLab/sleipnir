/*
  File autogenerated by gengetopt version 2.22
  generated with the following command:
  ..\..\..\..\extlib\proj\vs2005\release\gengetopt.exe --default-optional -N -e --output-dir=c:\Users\eblis\Documents\Visual Studio 2005\Projects\Sleipnir\trunk\tools\Dat2Dab\ -i ..\..\..\tools\Dat2Dab\Dat2Dab.ggo 

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

const char *gengetopt_args_info_purpose = "Text/binary data file interconversion";

const char *gengetopt_args_info_usage = "Usage: Dat2Dab [OPTIONS]...";

const char *gengetopt_args_info_description = "";

const char *gengetopt_args_info_help[] = {
  "  -h, --help              Print help and exit",
  "  -V, --version           Print version and exit",
  "\nMain:",
  "  -i, --input=filename    Input DAT/DAB file",
  "  -o, --output=filename   Output DAT/DAB file",
  "\nPreprocessing:",
  "  -f, --flip              Calculate one minus values  (default=off)",
  "  -n, --normalize         Normalize to the range [0,1]  (default=off)",
  "  -z, --zscore            Convert values to z-scores  (default=off)",
  "  -r, --rank              Rank transform data  (default=off)",
  "\nFiltering:",
  "  -g, --genes=filename    Process only genes from the given set",
  "  -G, --genex=filename    Exclude all genes from the given set",
  "  -c, --cutoff=DOUBLE     Exclude edges below cutoff",
  "  -e, --zero              Zero missing values  (default=off)",
  "  -d, --duplicates        Allow dissimilar duplicate values  (default=off)",
  "  -u, --subsample=FLOAT   Fraction of output to randomly subsample  \n                            (default=`1')",
  "\nLookups:",
  "  -l, --lookup1=STRING    First lookup gene",
  "  -L, --lookup2=STRING    Second lookup gene",
  "  -t, --lookups=filename  Lookup gene set",
  "  -T, --genelist          Only list genes  (default=off)",
  "  -P, --paircount         Only count pairs above cutoff  (default=off)",
  "\nOptional:",
  "  -p, --remap=filename    Gene name remapping file",
  "  -b, --table             Produce table formatted output  (default=off)",
  "  -s, --skip=INT          Columns to skip in input PCL  (default=`2')",
  "  -m, --memmap            Memory map input/output  (default=off)",
  "  -v, --verbosity=INT     Message verbosity  (default=`5')",
    0
};

typedef enum {ARG_NO
  , ARG_FLAG
  , ARG_STRING
  , ARG_INT
  , ARG_FLOAT
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
  args_info->output_given = 0 ;
  args_info->flip_given = 0 ;
  args_info->normalize_given = 0 ;
  args_info->zscore_given = 0 ;
  args_info->rank_given = 0 ;
  args_info->genes_given = 0 ;
  args_info->genex_given = 0 ;
  args_info->cutoff_given = 0 ;
  args_info->zero_given = 0 ;
  args_info->duplicates_given = 0 ;
  args_info->subsample_given = 0 ;
  args_info->lookup1_given = 0 ;
  args_info->lookup2_given = 0 ;
  args_info->lookups_given = 0 ;
  args_info->genelist_given = 0 ;
  args_info->paircount_given = 0 ;
  args_info->remap_given = 0 ;
  args_info->table_given = 0 ;
  args_info->skip_given = 0 ;
  args_info->memmap_given = 0 ;
  args_info->verbosity_given = 0 ;
}

static
void clear_args (struct gengetopt_args_info *args_info)
{
  args_info->input_arg = NULL;
  args_info->input_orig = NULL;
  args_info->output_arg = NULL;
  args_info->output_orig = NULL;
  args_info->flip_flag = 0;
  args_info->normalize_flag = 0;
  args_info->zscore_flag = 0;
  args_info->rank_flag = 0;
  args_info->genes_arg = NULL;
  args_info->genes_orig = NULL;
  args_info->genex_arg = NULL;
  args_info->genex_orig = NULL;
  args_info->cutoff_orig = NULL;
  args_info->zero_flag = 0;
  args_info->duplicates_flag = 0;
  args_info->subsample_arg = 1;
  args_info->subsample_orig = NULL;
  args_info->lookup1_arg = NULL;
  args_info->lookup1_orig = NULL;
  args_info->lookup2_arg = NULL;
  args_info->lookup2_orig = NULL;
  args_info->lookups_arg = NULL;
  args_info->lookups_orig = NULL;
  args_info->genelist_flag = 0;
  args_info->paircount_flag = 0;
  args_info->remap_arg = NULL;
  args_info->remap_orig = NULL;
  args_info->table_flag = 0;
  args_info->skip_arg = 2;
  args_info->skip_orig = NULL;
  args_info->memmap_flag = 0;
  args_info->verbosity_arg = 5;
  args_info->verbosity_orig = NULL;
  
}

static
void init_args_info(struct gengetopt_args_info *args_info)
{


  args_info->help_help = gengetopt_args_info_help[0] ;
  args_info->version_help = gengetopt_args_info_help[1] ;
  args_info->input_help = gengetopt_args_info_help[3] ;
  args_info->output_help = gengetopt_args_info_help[4] ;
  args_info->flip_help = gengetopt_args_info_help[6] ;
  args_info->normalize_help = gengetopt_args_info_help[7] ;
  args_info->zscore_help = gengetopt_args_info_help[8] ;
  args_info->rank_help = gengetopt_args_info_help[9] ;
  args_info->genes_help = gengetopt_args_info_help[11] ;
  args_info->genex_help = gengetopt_args_info_help[12] ;
  args_info->cutoff_help = gengetopt_args_info_help[13] ;
  args_info->zero_help = gengetopt_args_info_help[14] ;
  args_info->duplicates_help = gengetopt_args_info_help[15] ;
  args_info->subsample_help = gengetopt_args_info_help[16] ;
  args_info->lookup1_help = gengetopt_args_info_help[18] ;
  args_info->lookup2_help = gengetopt_args_info_help[19] ;
  args_info->lookups_help = gengetopt_args_info_help[20] ;
  args_info->genelist_help = gengetopt_args_info_help[21] ;
  args_info->paircount_help = gengetopt_args_info_help[22] ;
  args_info->remap_help = gengetopt_args_info_help[24] ;
  args_info->table_help = gengetopt_args_info_help[25] ;
  args_info->skip_help = gengetopt_args_info_help[26] ;
  args_info->memmap_help = gengetopt_args_info_help[27] ;
  args_info->verbosity_help = gengetopt_args_info_help[28] ;
  
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
  free_string_field (&(args_info->output_arg));
  free_string_field (&(args_info->output_orig));
  free_string_field (&(args_info->genes_arg));
  free_string_field (&(args_info->genes_orig));
  free_string_field (&(args_info->genex_arg));
  free_string_field (&(args_info->genex_orig));
  free_string_field (&(args_info->cutoff_orig));
  free_string_field (&(args_info->subsample_orig));
  free_string_field (&(args_info->lookup1_arg));
  free_string_field (&(args_info->lookup1_orig));
  free_string_field (&(args_info->lookup2_arg));
  free_string_field (&(args_info->lookup2_orig));
  free_string_field (&(args_info->lookups_arg));
  free_string_field (&(args_info->lookups_orig));
  free_string_field (&(args_info->remap_arg));
  free_string_field (&(args_info->remap_orig));
  free_string_field (&(args_info->skip_orig));
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
  if (args_info->output_given)
    write_into_file(outfile, "output", args_info->output_orig, 0);
  if (args_info->flip_given)
    write_into_file(outfile, "flip", 0, 0 );
  if (args_info->normalize_given)
    write_into_file(outfile, "normalize", 0, 0 );
  if (args_info->zscore_given)
    write_into_file(outfile, "zscore", 0, 0 );
  if (args_info->rank_given)
    write_into_file(outfile, "rank", 0, 0 );
  if (args_info->genes_given)
    write_into_file(outfile, "genes", args_info->genes_orig, 0);
  if (args_info->genex_given)
    write_into_file(outfile, "genex", args_info->genex_orig, 0);
  if (args_info->cutoff_given)
    write_into_file(outfile, "cutoff", args_info->cutoff_orig, 0);
  if (args_info->zero_given)
    write_into_file(outfile, "zero", 0, 0 );
  if (args_info->duplicates_given)
    write_into_file(outfile, "duplicates", 0, 0 );
  if (args_info->subsample_given)
    write_into_file(outfile, "subsample", args_info->subsample_orig, 0);
  if (args_info->lookup1_given)
    write_into_file(outfile, "lookup1", args_info->lookup1_orig, 0);
  if (args_info->lookup2_given)
    write_into_file(outfile, "lookup2", args_info->lookup2_orig, 0);
  if (args_info->lookups_given)
    write_into_file(outfile, "lookups", args_info->lookups_orig, 0);
  if (args_info->genelist_given)
    write_into_file(outfile, "genelist", 0, 0 );
  if (args_info->paircount_given)
    write_into_file(outfile, "paircount", 0, 0 );
  if (args_info->remap_given)
    write_into_file(outfile, "remap", args_info->remap_orig, 0);
  if (args_info->table_given)
    write_into_file(outfile, "table", 0, 0 );
  if (args_info->skip_given)
    write_into_file(outfile, "skip", args_info->skip_orig, 0);
  if (args_info->memmap_given)
    write_into_file(outfile, "memmap", 0, 0 );
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
  case ARG_FLAG:
    *((int *)field) = !*((int *)field);
    break;
  case ARG_INT:
    if (val) *((int *)field) = strtol (val, &stop_char, 0);
    break;
  case ARG_FLOAT:
    if (val) *((float *)field) = (float)strtod (val, &stop_char);
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
  case ARG_FLOAT:
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
        { "input",	1, NULL, 'i' },
        { "output",	1, NULL, 'o' },
        { "flip",	0, NULL, 'f' },
        { "normalize",	0, NULL, 'n' },
        { "zscore",	0, NULL, 'z' },
        { "rank",	0, NULL, 'r' },
        { "genes",	1, NULL, 'g' },
        { "genex",	1, NULL, 'G' },
        { "cutoff",	1, NULL, 'c' },
        { "zero",	0, NULL, 'e' },
        { "duplicates",	0, NULL, 'd' },
        { "subsample",	1, NULL, 'u' },
        { "lookup1",	1, NULL, 'l' },
        { "lookup2",	1, NULL, 'L' },
        { "lookups",	1, NULL, 't' },
        { "genelist",	0, NULL, 'T' },
        { "paircount",	0, NULL, 'P' },
        { "remap",	1, NULL, 'p' },
        { "table",	0, NULL, 'b' },
        { "skip",	1, NULL, 's' },
        { "memmap",	0, NULL, 'm' },
        { "verbosity",	1, NULL, 'v' },
        { NULL,	0, NULL, 0 }
      };

      c = getopt_long (argc, argv, "hVi:o:fnzrg:G:c:edu:l:L:t:TPp:bs:mv:", long_options, &option_index);

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
        case 'i':	/* Input DAT/DAB file.  */
        
        
          if (update_arg( (void *)&(args_info->input_arg), 
               &(args_info->input_orig), &(args_info->input_given),
              &(local_args_info.input_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "input", 'i',
              additional_error))
            goto failure;
        
          break;
        case 'o':	/* Output DAT/DAB file.  */
        
        
          if (update_arg( (void *)&(args_info->output_arg), 
               &(args_info->output_orig), &(args_info->output_given),
              &(local_args_info.output_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "output", 'o',
              additional_error))
            goto failure;
        
          break;
        case 'f':	/* Calculate one minus values.  */
        
        
          if (update_arg((void *)&(args_info->flip_flag), 0, &(args_info->flip_given),
              &(local_args_info.flip_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "flip", 'f',
              additional_error))
            goto failure;
        
          break;
        case 'n':	/* Normalize to the range [0,1].  */
        
        
          if (update_arg((void *)&(args_info->normalize_flag), 0, &(args_info->normalize_given),
              &(local_args_info.normalize_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "normalize", 'n',
              additional_error))
            goto failure;
        
          break;
        case 'z':	/* Convert values to z-scores.  */
        
        
          if (update_arg((void *)&(args_info->zscore_flag), 0, &(args_info->zscore_given),
              &(local_args_info.zscore_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "zscore", 'z',
              additional_error))
            goto failure;
        
          break;
        case 'r':	/* Rank transform data.  */
        
        
          if (update_arg((void *)&(args_info->rank_flag), 0, &(args_info->rank_given),
              &(local_args_info.rank_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "rank", 'r',
              additional_error))
            goto failure;
        
          break;
        case 'g':	/* Process only genes from the given set.  */
        
        
          if (update_arg( (void *)&(args_info->genes_arg), 
               &(args_info->genes_orig), &(args_info->genes_given),
              &(local_args_info.genes_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "genes", 'g',
              additional_error))
            goto failure;
        
          break;
        case 'G':	/* Exclude all genes from the given set.  */
        
        
          if (update_arg( (void *)&(args_info->genex_arg), 
               &(args_info->genex_orig), &(args_info->genex_given),
              &(local_args_info.genex_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "genex", 'G',
              additional_error))
            goto failure;
        
          break;
        case 'c':	/* Exclude edges below cutoff.  */
        
        
          if (update_arg( (void *)&(args_info->cutoff_arg), 
               &(args_info->cutoff_orig), &(args_info->cutoff_given),
              &(local_args_info.cutoff_given), optarg, 0, 0, ARG_DOUBLE,
              check_ambiguity, override, 0, 0,
              "cutoff", 'c',
              additional_error))
            goto failure;
        
          break;
        case 'e':	/* Zero missing values.  */
        
        
          if (update_arg((void *)&(args_info->zero_flag), 0, &(args_info->zero_given),
              &(local_args_info.zero_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "zero", 'e',
              additional_error))
            goto failure;
        
          break;
        case 'd':	/* Allow dissimilar duplicate values.  */
        
        
          if (update_arg((void *)&(args_info->duplicates_flag), 0, &(args_info->duplicates_given),
              &(local_args_info.duplicates_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "duplicates", 'd',
              additional_error))
            goto failure;
        
          break;
        case 'u':	/* Fraction of output to randomly subsample.  */
        
        
          if (update_arg( (void *)&(args_info->subsample_arg), 
               &(args_info->subsample_orig), &(args_info->subsample_given),
              &(local_args_info.subsample_given), optarg, 0, "1", ARG_FLOAT,
              check_ambiguity, override, 0, 0,
              "subsample", 'u',
              additional_error))
            goto failure;
        
          break;
        case 'l':	/* First lookup gene.  */
        
        
          if (update_arg( (void *)&(args_info->lookup1_arg), 
               &(args_info->lookup1_orig), &(args_info->lookup1_given),
              &(local_args_info.lookup1_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "lookup1", 'l',
              additional_error))
            goto failure;
        
          break;
        case 'L':	/* Second lookup gene.  */
        
        
          if (update_arg( (void *)&(args_info->lookup2_arg), 
               &(args_info->lookup2_orig), &(args_info->lookup2_given),
              &(local_args_info.lookup2_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "lookup2", 'L',
              additional_error))
            goto failure;
        
          break;
        case 't':	/* Lookup gene set.  */
        
        
          if (update_arg( (void *)&(args_info->lookups_arg), 
               &(args_info->lookups_orig), &(args_info->lookups_given),
              &(local_args_info.lookups_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "lookups", 't',
              additional_error))
            goto failure;
        
          break;
        case 'T':	/* Only list genes.  */
        
        
          if (update_arg((void *)&(args_info->genelist_flag), 0, &(args_info->genelist_given),
              &(local_args_info.genelist_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "genelist", 'T',
              additional_error))
            goto failure;
        
          break;
        case 'P':	/* Only count pairs above cutoff.  */
        
        
          if (update_arg((void *)&(args_info->paircount_flag), 0, &(args_info->paircount_given),
              &(local_args_info.paircount_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "paircount", 'P',
              additional_error))
            goto failure;
        
          break;
        case 'p':	/* Gene name remapping file.  */
        
        
          if (update_arg( (void *)&(args_info->remap_arg), 
               &(args_info->remap_orig), &(args_info->remap_given),
              &(local_args_info.remap_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "remap", 'p',
              additional_error))
            goto failure;
        
          break;
        case 'b':	/* Produce table formatted output.  */
        
        
          if (update_arg((void *)&(args_info->table_flag), 0, &(args_info->table_given),
              &(local_args_info.table_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "table", 'b',
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
        case 'm':	/* Memory map input/output.  */
        
        
          if (update_arg((void *)&(args_info->memmap_flag), 0, &(args_info->memmap_given),
              &(local_args_info.memmap_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "memmap", 'm',
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
