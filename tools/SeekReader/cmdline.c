/*
  File autogenerated by gengetopt version 2.22.5
  generated with the following command:
  /usr/bin/gengetopt -iSeekReader.ggo --default-optional -u -N -e 

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

const char *gengetopt_args_info_purpose = "Reads db files";

const char *gengetopt_args_info_usage = "Usage: SeekReader [OPTIONS]... [FILES]...";

const char *gengetopt_args_info_description = "";

const char *gengetopt_args_info_help[] = {
  "      --help                    Print help and exit",
  "  -V, --version                 Print version and exit",
  "\nDiagnosis:",
  "  -D, --databaselet             Display values from databaselet(s)  \n                                  (default=off)",
  "  -A, --dataset                 Check which datasets contain query of interest, \n                                  based on .gpres file  (default=off)",
  "  -W, --weight                  Test dataset weights  (default=off)",
  "  -U, --weight2                 Test dataset weights 2  (default=off)",
  "  -C, --comp_ranking            Compare two rankings (*.gscore files)  \n                                  (default=off)",
  "\nWeight:",
  "  -E, --dweight_dir=directory   Dataset weight directory  (default=`NA')",
  "  -n, --dweight_num=INT         Number of .dweight files  (default=`1000')",
  "  -M, --dweight_map=filename    Dataset mapping file  (default=`NA')",
  "  -F, --dweight_test_dir=directory\n                                Test dataset weight directory  (default=`NA')",
  "  -G, --dweight_test_num=INT    Test number of .dweight files  (default=`1000')",
  "\nCompare Rankings:",
  "  -H, --gscore_dir1=directory   Gene score directory 1  (default=`NA')",
  "  -h, --gscore_dir2=directory   Gene score directory 2  (default=`NA')",
  "  -I, --gscore_num1=INT         Number of .gscore files  (default=`1000')",
  "\nMain:",
  "  -O, --order_stat_single_gene_query\n                                Order statistics mode (single-gene query)  \n                                  (default=off)",
  "  -x, --db=filename             Input dataset-platform definition",
  "  -X, --dset_list=filename      Input a set of datasets",
  "  -i, --input=filename          Input gene mapping",
  "  -q, --single_query=filename   Query gene list",
  "  -d, --dir_in=directory        Database directory",
  "  -p, --dir_prep_in=directory   Prep directory (containing .gavg, .gpres files)",
  "  -r, --dir_gvar_in=directory   Prep directory (containing .gexpvar files)  \n                                  (default=`NA')",
  "  -s, --dir_sinfo_in=directory  Sinfo directory (containing .sinfo files)  \n                                  (default=`NA')",
  "  -N, --is_nibble               Whether the input DB is nibble type  \n                                  (default=off)",
  "  -P, --platform_dir=directory  Platform directory",
  "  -v, --gvar_cutoff=FLOAT       Query gene's variance in the dataset cutoff  \n                                  (default=`-1')",
  "  -Q, --multi_query=filename    File containing multiple queries  \n                                  (default=`NA')",
  "  -o, --output_file=filename    Output file  (default=`NA')",
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


static char *
gengetopt_strdup (const char *s);

static
void clear_given (struct gengetopt_args_info *args_info)
{
  args_info->help_given = 0 ;
  args_info->version_given = 0 ;
  args_info->databaselet_given = 0 ;
  args_info->dataset_given = 0 ;
  args_info->weight_given = 0 ;
  args_info->weight2_given = 0 ;
  args_info->comp_ranking_given = 0 ;
  args_info->dweight_dir_given = 0 ;
  args_info->dweight_num_given = 0 ;
  args_info->dweight_map_given = 0 ;
  args_info->dweight_test_dir_given = 0 ;
  args_info->dweight_test_num_given = 0 ;
  args_info->gscore_dir1_given = 0 ;
  args_info->gscore_dir2_given = 0 ;
  args_info->gscore_num1_given = 0 ;
  args_info->order_stat_single_gene_query_given = 0 ;
  args_info->db_given = 0 ;
  args_info->dset_list_given = 0 ;
  args_info->input_given = 0 ;
  args_info->single_query_given = 0 ;
  args_info->dir_in_given = 0 ;
  args_info->dir_prep_in_given = 0 ;
  args_info->dir_gvar_in_given = 0 ;
  args_info->dir_sinfo_in_given = 0 ;
  args_info->is_nibble_given = 0 ;
  args_info->platform_dir_given = 0 ;
  args_info->gvar_cutoff_given = 0 ;
  args_info->multi_query_given = 0 ;
  args_info->output_file_given = 0 ;
}

static
void clear_args (struct gengetopt_args_info *args_info)
{
  FIX_UNUSED (args_info);
  args_info->databaselet_flag = 0;
  args_info->dataset_flag = 0;
  args_info->weight_flag = 0;
  args_info->weight2_flag = 0;
  args_info->comp_ranking_flag = 0;
  args_info->dweight_dir_arg = gengetopt_strdup ("NA");
  args_info->dweight_dir_orig = NULL;
  args_info->dweight_num_arg = 1000;
  args_info->dweight_num_orig = NULL;
  args_info->dweight_map_arg = gengetopt_strdup ("NA");
  args_info->dweight_map_orig = NULL;
  args_info->dweight_test_dir_arg = gengetopt_strdup ("NA");
  args_info->dweight_test_dir_orig = NULL;
  args_info->dweight_test_num_arg = 1000;
  args_info->dweight_test_num_orig = NULL;
  args_info->gscore_dir1_arg = gengetopt_strdup ("NA");
  args_info->gscore_dir1_orig = NULL;
  args_info->gscore_dir2_arg = gengetopt_strdup ("NA");
  args_info->gscore_dir2_orig = NULL;
  args_info->gscore_num1_arg = 1000;
  args_info->gscore_num1_orig = NULL;
  args_info->order_stat_single_gene_query_flag = 0;
  args_info->db_arg = NULL;
  args_info->db_orig = NULL;
  args_info->dset_list_arg = NULL;
  args_info->dset_list_orig = NULL;
  args_info->input_arg = NULL;
  args_info->input_orig = NULL;
  args_info->single_query_arg = NULL;
  args_info->single_query_orig = NULL;
  args_info->dir_in_arg = NULL;
  args_info->dir_in_orig = NULL;
  args_info->dir_prep_in_arg = NULL;
  args_info->dir_prep_in_orig = NULL;
  args_info->dir_gvar_in_arg = gengetopt_strdup ("NA");
  args_info->dir_gvar_in_orig = NULL;
  args_info->dir_sinfo_in_arg = gengetopt_strdup ("NA");
  args_info->dir_sinfo_in_orig = NULL;
  args_info->is_nibble_flag = 0;
  args_info->platform_dir_arg = NULL;
  args_info->platform_dir_orig = NULL;
  args_info->gvar_cutoff_arg = -1;
  args_info->gvar_cutoff_orig = NULL;
  args_info->multi_query_arg = gengetopt_strdup ("NA");
  args_info->multi_query_orig = NULL;
  args_info->output_file_arg = gengetopt_strdup ("NA");
  args_info->output_file_orig = NULL;
  
}

static
void init_args_info(struct gengetopt_args_info *args_info)
{


  args_info->help_help = gengetopt_args_info_help[0] ;
  args_info->version_help = gengetopt_args_info_help[1] ;
  args_info->databaselet_help = gengetopt_args_info_help[3] ;
  args_info->dataset_help = gengetopt_args_info_help[4] ;
  args_info->weight_help = gengetopt_args_info_help[5] ;
  args_info->weight2_help = gengetopt_args_info_help[6] ;
  args_info->comp_ranking_help = gengetopt_args_info_help[7] ;
  args_info->dweight_dir_help = gengetopt_args_info_help[9] ;
  args_info->dweight_num_help = gengetopt_args_info_help[10] ;
  args_info->dweight_map_help = gengetopt_args_info_help[11] ;
  args_info->dweight_test_dir_help = gengetopt_args_info_help[12] ;
  args_info->dweight_test_num_help = gengetopt_args_info_help[13] ;
  args_info->gscore_dir1_help = gengetopt_args_info_help[15] ;
  args_info->gscore_dir2_help = gengetopt_args_info_help[16] ;
  args_info->gscore_num1_help = gengetopt_args_info_help[17] ;
  args_info->order_stat_single_gene_query_help = gengetopt_args_info_help[19] ;
  args_info->db_help = gengetopt_args_info_help[20] ;
  args_info->dset_list_help = gengetopt_args_info_help[21] ;
  args_info->input_help = gengetopt_args_info_help[22] ;
  args_info->single_query_help = gengetopt_args_info_help[23] ;
  args_info->dir_in_help = gengetopt_args_info_help[24] ;
  args_info->dir_prep_in_help = gengetopt_args_info_help[25] ;
  args_info->dir_gvar_in_help = gengetopt_args_info_help[26] ;
  args_info->dir_sinfo_in_help = gengetopt_args_info_help[27] ;
  args_info->is_nibble_help = gengetopt_args_info_help[28] ;
  args_info->platform_dir_help = gengetopt_args_info_help[29] ;
  args_info->gvar_cutoff_help = gengetopt_args_info_help[30] ;
  args_info->multi_query_help = gengetopt_args_info_help[31] ;
  args_info->output_file_help = gengetopt_args_info_help[32] ;
  
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
  free_string_field (&(args_info->dweight_dir_arg));
  free_string_field (&(args_info->dweight_dir_orig));
  free_string_field (&(args_info->dweight_num_orig));
  free_string_field (&(args_info->dweight_map_arg));
  free_string_field (&(args_info->dweight_map_orig));
  free_string_field (&(args_info->dweight_test_dir_arg));
  free_string_field (&(args_info->dweight_test_dir_orig));
  free_string_field (&(args_info->dweight_test_num_orig));
  free_string_field (&(args_info->gscore_dir1_arg));
  free_string_field (&(args_info->gscore_dir1_orig));
  free_string_field (&(args_info->gscore_dir2_arg));
  free_string_field (&(args_info->gscore_dir2_orig));
  free_string_field (&(args_info->gscore_num1_orig));
  free_string_field (&(args_info->db_arg));
  free_string_field (&(args_info->db_orig));
  free_string_field (&(args_info->dset_list_arg));
  free_string_field (&(args_info->dset_list_orig));
  free_string_field (&(args_info->input_arg));
  free_string_field (&(args_info->input_orig));
  free_string_field (&(args_info->single_query_arg));
  free_string_field (&(args_info->single_query_orig));
  free_string_field (&(args_info->dir_in_arg));
  free_string_field (&(args_info->dir_in_orig));
  free_string_field (&(args_info->dir_prep_in_arg));
  free_string_field (&(args_info->dir_prep_in_orig));
  free_string_field (&(args_info->dir_gvar_in_arg));
  free_string_field (&(args_info->dir_gvar_in_orig));
  free_string_field (&(args_info->dir_sinfo_in_arg));
  free_string_field (&(args_info->dir_sinfo_in_orig));
  free_string_field (&(args_info->platform_dir_arg));
  free_string_field (&(args_info->platform_dir_orig));
  free_string_field (&(args_info->gvar_cutoff_orig));
  free_string_field (&(args_info->multi_query_arg));
  free_string_field (&(args_info->multi_query_orig));
  free_string_field (&(args_info->output_file_arg));
  free_string_field (&(args_info->output_file_orig));
  
  
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
  if (args_info->databaselet_given)
    write_into_file(outfile, "databaselet", 0, 0 );
  if (args_info->dataset_given)
    write_into_file(outfile, "dataset", 0, 0 );
  if (args_info->weight_given)
    write_into_file(outfile, "weight", 0, 0 );
  if (args_info->weight2_given)
    write_into_file(outfile, "weight2", 0, 0 );
  if (args_info->comp_ranking_given)
    write_into_file(outfile, "comp_ranking", 0, 0 );
  if (args_info->dweight_dir_given)
    write_into_file(outfile, "dweight_dir", args_info->dweight_dir_orig, 0);
  if (args_info->dweight_num_given)
    write_into_file(outfile, "dweight_num", args_info->dweight_num_orig, 0);
  if (args_info->dweight_map_given)
    write_into_file(outfile, "dweight_map", args_info->dweight_map_orig, 0);
  if (args_info->dweight_test_dir_given)
    write_into_file(outfile, "dweight_test_dir", args_info->dweight_test_dir_orig, 0);
  if (args_info->dweight_test_num_given)
    write_into_file(outfile, "dweight_test_num", args_info->dweight_test_num_orig, 0);
  if (args_info->gscore_dir1_given)
    write_into_file(outfile, "gscore_dir1", args_info->gscore_dir1_orig, 0);
  if (args_info->gscore_dir2_given)
    write_into_file(outfile, "gscore_dir2", args_info->gscore_dir2_orig, 0);
  if (args_info->gscore_num1_given)
    write_into_file(outfile, "gscore_num1", args_info->gscore_num1_orig, 0);
  if (args_info->order_stat_single_gene_query_given)
    write_into_file(outfile, "order_stat_single_gene_query", 0, 0 );
  if (args_info->db_given)
    write_into_file(outfile, "db", args_info->db_orig, 0);
  if (args_info->dset_list_given)
    write_into_file(outfile, "dset_list", args_info->dset_list_orig, 0);
  if (args_info->input_given)
    write_into_file(outfile, "input", args_info->input_orig, 0);
  if (args_info->single_query_given)
    write_into_file(outfile, "single_query", args_info->single_query_orig, 0);
  if (args_info->dir_in_given)
    write_into_file(outfile, "dir_in", args_info->dir_in_orig, 0);
  if (args_info->dir_prep_in_given)
    write_into_file(outfile, "dir_prep_in", args_info->dir_prep_in_orig, 0);
  if (args_info->dir_gvar_in_given)
    write_into_file(outfile, "dir_gvar_in", args_info->dir_gvar_in_orig, 0);
  if (args_info->dir_sinfo_in_given)
    write_into_file(outfile, "dir_sinfo_in", args_info->dir_sinfo_in_orig, 0);
  if (args_info->is_nibble_given)
    write_into_file(outfile, "is_nibble", 0, 0 );
  if (args_info->platform_dir_given)
    write_into_file(outfile, "platform_dir", args_info->platform_dir_orig, 0);
  if (args_info->gvar_cutoff_given)
    write_into_file(outfile, "gvar_cutoff", args_info->gvar_cutoff_orig, 0);
  if (args_info->multi_query_given)
    write_into_file(outfile, "multi_query", args_info->multi_query_orig, 0);
  if (args_info->output_file_given)
    write_into_file(outfile, "output_file", args_info->output_file_orig, 0);
  

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
        { "help",	0, NULL, 0 },
        { "version",	0, NULL, 'V' },
        { "databaselet",	0, NULL, 'D' },
        { "dataset",	0, NULL, 'A' },
        { "weight",	0, NULL, 'W' },
        { "weight2",	0, NULL, 'U' },
        { "comp_ranking",	0, NULL, 'C' },
        { "dweight_dir",	1, NULL, 'E' },
        { "dweight_num",	1, NULL, 'n' },
        { "dweight_map",	1, NULL, 'M' },
        { "dweight_test_dir",	1, NULL, 'F' },
        { "dweight_test_num",	1, NULL, 'G' },
        { "gscore_dir1",	1, NULL, 'H' },
        { "gscore_dir2",	1, NULL, 'h' },
        { "gscore_num1",	1, NULL, 'I' },
        { "order_stat_single_gene_query",	0, NULL, 'O' },
        { "db",	1, NULL, 'x' },
        { "dset_list",	1, NULL, 'X' },
        { "input",	1, NULL, 'i' },
        { "single_query",	1, NULL, 'q' },
        { "dir_in",	1, NULL, 'd' },
        { "dir_prep_in",	1, NULL, 'p' },
        { "dir_gvar_in",	1, NULL, 'r' },
        { "dir_sinfo_in",	1, NULL, 's' },
        { "is_nibble",	0, NULL, 'N' },
        { "platform_dir",	1, NULL, 'P' },
        { "gvar_cutoff",	1, NULL, 'v' },
        { "multi_query",	1, NULL, 'Q' },
        { "output_file",	1, NULL, 'o' },
        { 0,  0, 0, 0 }
      };

      c = getopt_long (argc, argv, "VDAWUCE:n:M:F:G:H:h:I:Ox:X:i:q:d:p:r:s:NP:v:Q:o:", long_options, &option_index);

      if (c == -1) break;	/* Exit from `while (1)' loop.  */

      switch (c)
        {
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
        case 'D':	/* Display values from databaselet(s).  */
        
        
          if (update_arg((void *)&(args_info->databaselet_flag), 0, &(args_info->databaselet_given),
              &(local_args_info.databaselet_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "databaselet", 'D',
              additional_error))
            goto failure;
        
          break;
        case 'A':	/* Check which datasets contain query of interest, based on .gpres file.  */
        
        
          if (update_arg((void *)&(args_info->dataset_flag), 0, &(args_info->dataset_given),
              &(local_args_info.dataset_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "dataset", 'A',
              additional_error))
            goto failure;
        
          break;
        case 'W':	/* Test dataset weights.  */
        
        
          if (update_arg((void *)&(args_info->weight_flag), 0, &(args_info->weight_given),
              &(local_args_info.weight_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "weight", 'W',
              additional_error))
            goto failure;
        
          break;
        case 'U':	/* Test dataset weights 2.  */
        
        
          if (update_arg((void *)&(args_info->weight2_flag), 0, &(args_info->weight2_given),
              &(local_args_info.weight2_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "weight2", 'U',
              additional_error))
            goto failure;
        
          break;
        case 'C':	/* Compare two rankings (*.gscore files).  */
        
        
          if (update_arg((void *)&(args_info->comp_ranking_flag), 0, &(args_info->comp_ranking_given),
              &(local_args_info.comp_ranking_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "comp_ranking", 'C',
              additional_error))
            goto failure;
        
          break;
        case 'E':	/* Dataset weight directory.  */
        
        
          if (update_arg( (void *)&(args_info->dweight_dir_arg), 
               &(args_info->dweight_dir_orig), &(args_info->dweight_dir_given),
              &(local_args_info.dweight_dir_given), optarg, 0, "NA", ARG_STRING,
              check_ambiguity, override, 0, 0,
              "dweight_dir", 'E',
              additional_error))
            goto failure;
        
          break;
        case 'n':	/* Number of .dweight files.  */
        
        
          if (update_arg( (void *)&(args_info->dweight_num_arg), 
               &(args_info->dweight_num_orig), &(args_info->dweight_num_given),
              &(local_args_info.dweight_num_given), optarg, 0, "1000", ARG_INT,
              check_ambiguity, override, 0, 0,
              "dweight_num", 'n',
              additional_error))
            goto failure;
        
          break;
        case 'M':	/* Dataset mapping file.  */
        
        
          if (update_arg( (void *)&(args_info->dweight_map_arg), 
               &(args_info->dweight_map_orig), &(args_info->dweight_map_given),
              &(local_args_info.dweight_map_given), optarg, 0, "NA", ARG_STRING,
              check_ambiguity, override, 0, 0,
              "dweight_map", 'M',
              additional_error))
            goto failure;
        
          break;
        case 'F':	/* Test dataset weight directory.  */
        
        
          if (update_arg( (void *)&(args_info->dweight_test_dir_arg), 
               &(args_info->dweight_test_dir_orig), &(args_info->dweight_test_dir_given),
              &(local_args_info.dweight_test_dir_given), optarg, 0, "NA", ARG_STRING,
              check_ambiguity, override, 0, 0,
              "dweight_test_dir", 'F',
              additional_error))
            goto failure;
        
          break;
        case 'G':	/* Test number of .dweight files.  */
        
        
          if (update_arg( (void *)&(args_info->dweight_test_num_arg), 
               &(args_info->dweight_test_num_orig), &(args_info->dweight_test_num_given),
              &(local_args_info.dweight_test_num_given), optarg, 0, "1000", ARG_INT,
              check_ambiguity, override, 0, 0,
              "dweight_test_num", 'G',
              additional_error))
            goto failure;
        
          break;
        case 'H':	/* Gene score directory 1.  */
        
        
          if (update_arg( (void *)&(args_info->gscore_dir1_arg), 
               &(args_info->gscore_dir1_orig), &(args_info->gscore_dir1_given),
              &(local_args_info.gscore_dir1_given), optarg, 0, "NA", ARG_STRING,
              check_ambiguity, override, 0, 0,
              "gscore_dir1", 'H',
              additional_error))
            goto failure;
        
          break;
        case 'h':	/* Gene score directory 2.  */
        
        
          if (update_arg( (void *)&(args_info->gscore_dir2_arg), 
               &(args_info->gscore_dir2_orig), &(args_info->gscore_dir2_given),
              &(local_args_info.gscore_dir2_given), optarg, 0, "NA", ARG_STRING,
              check_ambiguity, override, 0, 0,
              "gscore_dir2", 'h',
              additional_error))
            goto failure;
        
          break;
        case 'I':	/* Number of .gscore files.  */
        
        
          if (update_arg( (void *)&(args_info->gscore_num1_arg), 
               &(args_info->gscore_num1_orig), &(args_info->gscore_num1_given),
              &(local_args_info.gscore_num1_given), optarg, 0, "1000", ARG_INT,
              check_ambiguity, override, 0, 0,
              "gscore_num1", 'I',
              additional_error))
            goto failure;
        
          break;
        case 'O':	/* Order statistics mode (single-gene query).  */
        
        
          if (update_arg((void *)&(args_info->order_stat_single_gene_query_flag), 0, &(args_info->order_stat_single_gene_query_given),
              &(local_args_info.order_stat_single_gene_query_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "order_stat_single_gene_query", 'O',
              additional_error))
            goto failure;
        
          break;
        case 'x':	/* Input dataset-platform definition.  */
        
        
          if (update_arg( (void *)&(args_info->db_arg), 
               &(args_info->db_orig), &(args_info->db_given),
              &(local_args_info.db_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "db", 'x',
              additional_error))
            goto failure;
        
          break;
        case 'X':	/* Input a set of datasets.  */
        
        
          if (update_arg( (void *)&(args_info->dset_list_arg), 
               &(args_info->dset_list_orig), &(args_info->dset_list_given),
              &(local_args_info.dset_list_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "dset_list", 'X',
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
        case 'q':	/* Query gene list.  */
        
        
          if (update_arg( (void *)&(args_info->single_query_arg), 
               &(args_info->single_query_orig), &(args_info->single_query_given),
              &(local_args_info.single_query_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "single_query", 'q',
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
        case 'r':	/* Prep directory (containing .gexpvar files).  */
        
        
          if (update_arg( (void *)&(args_info->dir_gvar_in_arg), 
               &(args_info->dir_gvar_in_orig), &(args_info->dir_gvar_in_given),
              &(local_args_info.dir_gvar_in_given), optarg, 0, "NA", ARG_STRING,
              check_ambiguity, override, 0, 0,
              "dir_gvar_in", 'r',
              additional_error))
            goto failure;
        
          break;
        case 's':	/* Sinfo directory (containing .sinfo files).  */
        
        
          if (update_arg( (void *)&(args_info->dir_sinfo_in_arg), 
               &(args_info->dir_sinfo_in_orig), &(args_info->dir_sinfo_in_given),
              &(local_args_info.dir_sinfo_in_given), optarg, 0, "NA", ARG_STRING,
              check_ambiguity, override, 0, 0,
              "dir_sinfo_in", 's',
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
        case 'P':	/* Platform directory.  */
        
        
          if (update_arg( (void *)&(args_info->platform_dir_arg), 
               &(args_info->platform_dir_orig), &(args_info->platform_dir_given),
              &(local_args_info.platform_dir_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "platform_dir", 'P',
              additional_error))
            goto failure;
        
          break;
        case 'v':	/* Query gene's variance in the dataset cutoff.  */
        
        
          if (update_arg( (void *)&(args_info->gvar_cutoff_arg), 
               &(args_info->gvar_cutoff_orig), &(args_info->gvar_cutoff_given),
              &(local_args_info.gvar_cutoff_given), optarg, 0, "-1", ARG_FLOAT,
              check_ambiguity, override, 0, 0,
              "gvar_cutoff", 'v',
              additional_error))
            goto failure;
        
          break;
        case 'Q':	/* File containing multiple queries.  */
        
        
          if (update_arg( (void *)&(args_info->multi_query_arg), 
               &(args_info->multi_query_orig), &(args_info->multi_query_given),
              &(local_args_info.multi_query_given), optarg, 0, "NA", ARG_STRING,
              check_ambiguity, override, 0, 0,
              "multi_query", 'Q',
              additional_error))
            goto failure;
        
          break;
        case 'o':	/* Output file.  */
        
        
          if (update_arg( (void *)&(args_info->output_file_arg), 
               &(args_info->output_file_orig), &(args_info->output_file_given),
              &(local_args_info.output_file_given), optarg, 0, "NA", ARG_STRING,
              check_ambiguity, override, 0, 0,
              "output_file", 'o',
              additional_error))
            goto failure;
        
          break;

        case 0:	/* Long option with no short option */
          if (strcmp (long_options[option_index].name, "help") == 0) {
            cmdline_parser_print_help ();
            cmdline_parser_free (&local_args_info);
            exit (EXIT_SUCCESS);
          }

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
