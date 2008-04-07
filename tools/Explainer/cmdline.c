/*
  File autogenerated by gengetopt version 2.22
  generated with the following command:
  /Genomics/Users/mark/sleipnir/extlib/gengetopt-2.22/src/gengetopt -iExplainer.ggo --default-optional -C -N -e 

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

const char *gengetopt_args_info_purpose = "Allows evaluation of genes contributing to a good prediction score";

const char *gengetopt_args_info_usage = "Usage: Explainer [OPTIONS]...";

const char *gengetopt_args_info_description = "";

const char *gengetopt_args_info_help[] = {
  "  -h, --help               Print help and exit",
  "  -V, --version            Print version and exit",
  "\nMain:",
  "  -i, --input=filename     Similarity DAT/DAB file",
  "  -w, --answers=filename   Answer DAT/DAB file",
  "\nMiscellaneous:",
  "  -k, --count=INT          Number of pairs to display  (default=`-1')",
  "  -p, --positives          Include only positive pairs  (default=off)",
  "  -e, --everything         Include pairs without answers  (default=off)",
  "  -u, --unknowns=STRING    Treatment of unknown genes  (possible \n                             values=\"exclude\", \"include\", \"only\" \n                             default=`exclude')",
  "  -x, --fraction=DOUBLE    Random fraction of results to calculate  \n                             (default=`1')",
  "\nLearning/Evaluation:",
  "  -g, --genes=filename     Gene inclusion file",
  "  -G, --genex=filename     Gene exclusion file",
  "  -c, --genet=filename     Term inclusion file",
  "\nPreprocessing:",
  "  -n, --normalize          Normalize to the range [0,1]  (default=off)",
  "  -t, --invert             Invert correlations to distances  (default=off)",
  "  -r, --reverse            Reverse sort order  (default=off)",
  "\nFunction Catalogs:",
  "  -o, --go_onto=filename   GO ontology",
  "  -a, --go_anno=filename   GO annotations",
  "  -f, --features=filename  SGD gene features",
  "\nOptional:",
  "  -m, --memmap             Memory map input files  (default=off)",
  "  -C, --config=filename    Command line config file  (default=`Explainer.ini')",
  "  -v, --verbosity=INT      Message verbosity  (default=`5')",
    0
};

typedef enum {ARG_NO
  , ARG_FLAG
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

static int
cmdline_parser_required2 (struct gengetopt_args_info *args_info, const char *prog_name, const char *additional_error);
struct line_list
{
  char * string_arg;
  struct line_list * next;
};

static struct line_list *cmd_line_list = 0;
static struct line_list *cmd_line_list_tmp = 0;

static void
free_cmd_list(void)
{
  /* free the list of a previous call */
  if (cmd_line_list)
    {
      while (cmd_line_list) {
        cmd_line_list_tmp = cmd_line_list;
        cmd_line_list = cmd_line_list->next;
        free (cmd_line_list_tmp->string_arg);
        free (cmd_line_list_tmp);
      }
    }
}


char *cmdline_parser_unknowns_values[] = {"exclude", "include", "only", 0} ;	/* Possible values for unknowns.  */

static char *
gengetopt_strdup (const char *s);

static
void clear_given (struct gengetopt_args_info *args_info)
{
  args_info->help_given = 0 ;
  args_info->version_given = 0 ;
  args_info->input_given = 0 ;
  args_info->answers_given = 0 ;
  args_info->count_given = 0 ;
  args_info->positives_given = 0 ;
  args_info->everything_given = 0 ;
  args_info->unknowns_given = 0 ;
  args_info->fraction_given = 0 ;
  args_info->genes_given = 0 ;
  args_info->genex_given = 0 ;
  args_info->genet_given = 0 ;
  args_info->normalize_given = 0 ;
  args_info->invert_given = 0 ;
  args_info->reverse_given = 0 ;
  args_info->go_onto_given = 0 ;
  args_info->go_anno_given = 0 ;
  args_info->features_given = 0 ;
  args_info->memmap_given = 0 ;
  args_info->config_given = 0 ;
  args_info->verbosity_given = 0 ;
}

static
void clear_args (struct gengetopt_args_info *args_info)
{
  args_info->input_arg = NULL;
  args_info->input_orig = NULL;
  args_info->answers_arg = NULL;
  args_info->answers_orig = NULL;
  args_info->count_arg = -1;
  args_info->count_orig = NULL;
  args_info->positives_flag = 0;
  args_info->everything_flag = 0;
  args_info->unknowns_arg = gengetopt_strdup ("exclude");
  args_info->unknowns_orig = NULL;
  args_info->fraction_arg = 1;
  args_info->fraction_orig = NULL;
  args_info->genes_arg = NULL;
  args_info->genes_orig = NULL;
  args_info->genex_arg = NULL;
  args_info->genex_orig = NULL;
  args_info->genet_arg = NULL;
  args_info->genet_orig = NULL;
  args_info->normalize_flag = 0;
  args_info->invert_flag = 0;
  args_info->reverse_flag = 0;
  args_info->go_onto_arg = NULL;
  args_info->go_onto_orig = NULL;
  args_info->go_anno_arg = NULL;
  args_info->go_anno_orig = NULL;
  args_info->features_arg = NULL;
  args_info->features_orig = NULL;
  args_info->memmap_flag = 0;
  args_info->config_arg = gengetopt_strdup ("Explainer.ini");
  args_info->config_orig = NULL;
  args_info->verbosity_arg = 5;
  args_info->verbosity_orig = NULL;
  
}

static
void init_args_info(struct gengetopt_args_info *args_info)
{


  args_info->help_help = gengetopt_args_info_help[0] ;
  args_info->version_help = gengetopt_args_info_help[1] ;
  args_info->input_help = gengetopt_args_info_help[3] ;
  args_info->answers_help = gengetopt_args_info_help[4] ;
  args_info->count_help = gengetopt_args_info_help[6] ;
  args_info->positives_help = gengetopt_args_info_help[7] ;
  args_info->everything_help = gengetopt_args_info_help[8] ;
  args_info->unknowns_help = gengetopt_args_info_help[9] ;
  args_info->fraction_help = gengetopt_args_info_help[10] ;
  args_info->genes_help = gengetopt_args_info_help[12] ;
  args_info->genex_help = gengetopt_args_info_help[13] ;
  args_info->genet_help = gengetopt_args_info_help[14] ;
  args_info->normalize_help = gengetopt_args_info_help[16] ;
  args_info->invert_help = gengetopt_args_info_help[17] ;
  args_info->reverse_help = gengetopt_args_info_help[18] ;
  args_info->go_onto_help = gengetopt_args_info_help[20] ;
  args_info->go_anno_help = gengetopt_args_info_help[21] ;
  args_info->features_help = gengetopt_args_info_help[22] ;
  args_info->memmap_help = gengetopt_args_info_help[24] ;
  args_info->config_help = gengetopt_args_info_help[25] ;
  args_info->verbosity_help = gengetopt_args_info_help[26] ;
  
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
  free_string_field (&(args_info->answers_arg));
  free_string_field (&(args_info->answers_orig));
  free_string_field (&(args_info->count_orig));
  free_string_field (&(args_info->unknowns_arg));
  free_string_field (&(args_info->unknowns_orig));
  free_string_field (&(args_info->fraction_orig));
  free_string_field (&(args_info->genes_arg));
  free_string_field (&(args_info->genes_orig));
  free_string_field (&(args_info->genex_arg));
  free_string_field (&(args_info->genex_orig));
  free_string_field (&(args_info->genet_arg));
  free_string_field (&(args_info->genet_orig));
  free_string_field (&(args_info->go_onto_arg));
  free_string_field (&(args_info->go_onto_orig));
  free_string_field (&(args_info->go_anno_arg));
  free_string_field (&(args_info->go_anno_orig));
  free_string_field (&(args_info->features_arg));
  free_string_field (&(args_info->features_orig));
  free_string_field (&(args_info->config_arg));
  free_string_field (&(args_info->config_orig));
  free_string_field (&(args_info->verbosity_orig));
  
  

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
check_possible_values(const char *val, char *values[])
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
write_into_file(FILE *outfile, const char *opt, const char *arg, char *values[])
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
  if (args_info->input_given)
    write_into_file(outfile, "input", args_info->input_orig, 0);
  if (args_info->answers_given)
    write_into_file(outfile, "answers", args_info->answers_orig, 0);
  if (args_info->count_given)
    write_into_file(outfile, "count", args_info->count_orig, 0);
  if (args_info->positives_given)
    write_into_file(outfile, "positives", 0, 0 );
  if (args_info->everything_given)
    write_into_file(outfile, "everything", 0, 0 );
  if (args_info->unknowns_given)
    write_into_file(outfile, "unknowns", args_info->unknowns_orig, cmdline_parser_unknowns_values);
  if (args_info->fraction_given)
    write_into_file(outfile, "fraction", args_info->fraction_orig, 0);
  if (args_info->genes_given)
    write_into_file(outfile, "genes", args_info->genes_orig, 0);
  if (args_info->genex_given)
    write_into_file(outfile, "genex", args_info->genex_orig, 0);
  if (args_info->genet_given)
    write_into_file(outfile, "genet", args_info->genet_orig, 0);
  if (args_info->normalize_given)
    write_into_file(outfile, "normalize", 0, 0 );
  if (args_info->invert_given)
    write_into_file(outfile, "invert", 0, 0 );
  if (args_info->reverse_given)
    write_into_file(outfile, "reverse", 0, 0 );
  if (args_info->go_onto_given)
    write_into_file(outfile, "go_onto", args_info->go_onto_orig, 0);
  if (args_info->go_anno_given)
    write_into_file(outfile, "go_anno", args_info->go_anno_orig, 0);
  if (args_info->features_given)
    write_into_file(outfile, "features", args_info->features_orig, 0);
  if (args_info->memmap_given)
    write_into_file(outfile, "memmap", 0, 0 );
  if (args_info->config_given)
    write_into_file(outfile, "config", args_info->config_orig, 0);
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
  case ARG_FLAG:
    *((int *)field) = !*((int *)field);
    break;
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
        { "answers",	1, NULL, 'w' },
        { "count",	1, NULL, 'k' },
        { "positives",	0, NULL, 'p' },
        { "everything",	0, NULL, 'e' },
        { "unknowns",	1, NULL, 'u' },
        { "fraction",	1, NULL, 'x' },
        { "genes",	1, NULL, 'g' },
        { "genex",	1, NULL, 'G' },
        { "genet",	1, NULL, 'c' },
        { "normalize",	0, NULL, 'n' },
        { "invert",	0, NULL, 't' },
        { "reverse",	0, NULL, 'r' },
        { "go_onto",	1, NULL, 'o' },
        { "go_anno",	1, NULL, 'a' },
        { "features",	1, NULL, 'f' },
        { "memmap",	0, NULL, 'm' },
        { "config",	1, NULL, 'C' },
        { "verbosity",	1, NULL, 'v' },
        { NULL,	0, NULL, 0 }
      };

      c = getopt_long (argc, argv, "hVi:w:k:peu:x:g:G:c:ntro:a:f:mC:v:", long_options, &option_index);

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
        case 'i':	/* Similarity DAT/DAB file.  */
        
        
          if (update_arg( (void *)&(args_info->input_arg), 
               &(args_info->input_orig), &(args_info->input_given),
              &(local_args_info.input_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "input", 'i',
              additional_error))
            goto failure;
        
          break;
        case 'w':	/* Answer DAT/DAB file.  */
        
        
          if (update_arg( (void *)&(args_info->answers_arg), 
               &(args_info->answers_orig), &(args_info->answers_given),
              &(local_args_info.answers_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "answers", 'w',
              additional_error))
            goto failure;
        
          break;
        case 'k':	/* Number of pairs to display.  */
        
        
          if (update_arg( (void *)&(args_info->count_arg), 
               &(args_info->count_orig), &(args_info->count_given),
              &(local_args_info.count_given), optarg, 0, "-1", ARG_INT,
              check_ambiguity, override, 0, 0,
              "count", 'k',
              additional_error))
            goto failure;
        
          break;
        case 'p':	/* Include only positive pairs.  */
        
        
          if (update_arg((void *)&(args_info->positives_flag), 0, &(args_info->positives_given),
              &(local_args_info.positives_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "positives", 'p',
              additional_error))
            goto failure;
        
          break;
        case 'e':	/* Include pairs without answers.  */
        
        
          if (update_arg((void *)&(args_info->everything_flag), 0, &(args_info->everything_given),
              &(local_args_info.everything_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "everything", 'e',
              additional_error))
            goto failure;
        
          break;
        case 'u':	/* Treatment of unknown genes.  */
        
        
          if (update_arg( (void *)&(args_info->unknowns_arg), 
               &(args_info->unknowns_orig), &(args_info->unknowns_given),
              &(local_args_info.unknowns_given), optarg, cmdline_parser_unknowns_values, "exclude", ARG_STRING,
              check_ambiguity, override, 0, 0,
              "unknowns", 'u',
              additional_error))
            goto failure;
        
          break;
        case 'x':	/* Random fraction of results to calculate.  */
        
        
          if (update_arg( (void *)&(args_info->fraction_arg), 
               &(args_info->fraction_orig), &(args_info->fraction_given),
              &(local_args_info.fraction_given), optarg, 0, "1", ARG_DOUBLE,
              check_ambiguity, override, 0, 0,
              "fraction", 'x',
              additional_error))
            goto failure;
        
          break;
        case 'g':	/* Gene inclusion file.  */
        
        
          if (update_arg( (void *)&(args_info->genes_arg), 
               &(args_info->genes_orig), &(args_info->genes_given),
              &(local_args_info.genes_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "genes", 'g',
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
        case 'c':	/* Term inclusion file.  */
        
        
          if (update_arg( (void *)&(args_info->genet_arg), 
               &(args_info->genet_orig), &(args_info->genet_given),
              &(local_args_info.genet_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "genet", 'c',
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
        case 't':	/* Invert correlations to distances.  */
        
        
          if (update_arg((void *)&(args_info->invert_flag), 0, &(args_info->invert_given),
              &(local_args_info.invert_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "invert", 't',
              additional_error))
            goto failure;
        
          break;
        case 'r':	/* Reverse sort order.  */
        
        
          if (update_arg((void *)&(args_info->reverse_flag), 0, &(args_info->reverse_given),
              &(local_args_info.reverse_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "reverse", 'r',
              additional_error))
            goto failure;
        
          break;
        case 'o':	/* GO ontology.  */
        
        
          if (update_arg( (void *)&(args_info->go_onto_arg), 
               &(args_info->go_onto_orig), &(args_info->go_onto_given),
              &(local_args_info.go_onto_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "go_onto", 'o',
              additional_error))
            goto failure;
        
          break;
        case 'a':	/* GO annotations.  */
        
        
          if (update_arg( (void *)&(args_info->go_anno_arg), 
               &(args_info->go_anno_orig), &(args_info->go_anno_given),
              &(local_args_info.go_anno_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "go_anno", 'a',
              additional_error))
            goto failure;
        
          break;
        case 'f':	/* SGD gene features.  */
        
        
          if (update_arg( (void *)&(args_info->features_arg), 
               &(args_info->features_orig), &(args_info->features_given),
              &(local_args_info.features_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "features", 'f',
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
        case 'C':	/* Command line config file.  */
        
        
          if (update_arg( (void *)&(args_info->config_arg), 
               &(args_info->config_orig), &(args_info->config_given),
              &(local_args_info.config_given), optarg, 0, "Explainer.ini", ARG_STRING,
              check_ambiguity, override, 0, 0,
              "config", 'C',
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

  return 0;

failure:
  
  cmdline_parser_release (&local_args_info);
  return (EXIT_FAILURE);
}

#ifndef CONFIG_FILE_LINE_SIZE
#define CONFIG_FILE_LINE_SIZE 2048
#endif
#define ADDITIONAL_ERROR " in configuration file "

#define CONFIG_FILE_LINE_BUFFER_SIZE (CONFIG_FILE_LINE_SIZE+3)
/* 3 is for "--" and "=" */

static int
_cmdline_parser_configfile (char * const filename, int *my_argc)
{
  FILE* file;
  char my_argv[CONFIG_FILE_LINE_BUFFER_SIZE+1];
  char linebuf[CONFIG_FILE_LINE_SIZE];
  int line_num = 0;
  int result = 0, equal;
  char *fopt, *farg;
  char *str_index;
  size_t len, next_token;
  char delimiter;

  if ((file = fopen(filename, "r")) == NULL)
    {
      fprintf (stderr, "%s: Error opening configuration file '%s'\n",
               CMDLINE_PARSER_PACKAGE, filename);
      return EXIT_FAILURE;
    }

  while ((fgets(linebuf, CONFIG_FILE_LINE_SIZE, file)) != NULL)
    {
      ++line_num;
      my_argv[0] = '\0';
      len = strlen(linebuf);
      if (len > (CONFIG_FILE_LINE_BUFFER_SIZE-1))
        {
          fprintf (stderr, "%s:%s:%d: Line too long in configuration file\n",
                   CMDLINE_PARSER_PACKAGE, filename, line_num);
          result = EXIT_FAILURE;
          break;
        }

      /* find first non-whitespace character in the line */
      next_token = strspn (linebuf, " \t\r\n");
      str_index  = linebuf + next_token;

      if ( str_index[0] == '\0' || str_index[0] == '#')
        continue; /* empty line or comment line is skipped */

      fopt = str_index;

      /* truncate fopt at the end of the first non-valid character */
      next_token = strcspn (fopt, " \t\r\n=");

      if (fopt[next_token] == '\0') /* the line is over */
        {
          farg  = NULL;
          equal = 0;
          goto noarg;
        }

      /* remember if equal sign is present */
      equal = (fopt[next_token] == '=');
      fopt[next_token++] = '\0';

      /* advance pointers to the next token after the end of fopt */
      next_token += strspn (fopt + next_token, " \t\r\n");

      /* check for the presence of equal sign, and if so, skip it */
      if ( !equal )
        if ((equal = (fopt[next_token] == '=')))
          {
            next_token++;
            next_token += strspn (fopt + next_token, " \t\r\n");
          }
      str_index  += next_token;

      /* find argument */
      farg = str_index;
      if ( farg[0] == '\"' || farg[0] == '\'' )
        { /* quoted argument */
          str_index = strchr (++farg, str_index[0] ); /* skip opening quote */
          if (! str_index)
            {
              fprintf
                (stderr,
                 "%s:%s:%d: unterminated string in configuration file\n",
                 CMDLINE_PARSER_PACKAGE, filename, line_num);
              result = EXIT_FAILURE;
              break;
            }
        }
      else
        { /* read up the remaining part up to a delimiter */
          next_token = strcspn (farg, " \t\r\n#\'\"");
          str_index += next_token;
        }

      /* truncate farg at the delimiter and store it for further check */
      delimiter = *str_index, *str_index++ = '\0';

      /* everything but comment is illegal at the end of line */
      if (delimiter != '\0' && delimiter != '#')
        {
          str_index += strspn(str_index, " \t\r\n");
          if (*str_index != '\0' && *str_index != '#')
            {
              fprintf
                (stderr,
                 "%s:%s:%d: malformed string in configuration file\n",
                 CMDLINE_PARSER_PACKAGE, filename, line_num);
              result = EXIT_FAILURE;
              break;
            }
        }

    noarg:
      if (!strcmp(fopt,"include")) {
        if (farg && *farg) {
          result = _cmdline_parser_configfile(farg, my_argc);
        } else {
          fprintf(stderr, "%s:%s:%d: include requires a filename argument.\n",
                  CMDLINE_PARSER_PACKAGE, filename, line_num);
        }
        continue;
      }
      len = strlen(fopt);
      strcat (my_argv, len > 1 ? "--" : "-");
      strcat (my_argv, fopt);
      if (len > 1 && ((farg && *farg) || equal))
        strcat (my_argv, "=");
      if (farg && *farg)
        strcat (my_argv, farg);
      ++(*my_argc);

      cmd_line_list_tmp = (struct line_list *) malloc (sizeof (struct line_list));
      cmd_line_list_tmp->next = cmd_line_list;
      cmd_line_list = cmd_line_list_tmp;
      cmd_line_list->string_arg = gengetopt_strdup(my_argv);
    } /* while */

  if (file)
    fclose(file);
  return result;
}

int
cmdline_parser_configfile (char * const filename,
                           struct gengetopt_args_info *args_info,
                           int override, int initialize, int check_required)
{
  struct cmdline_parser_params params;

  params.override = override;
  params.initialize = initialize;
  params.check_required = check_required;
  params.check_ambiguity = 0;
  params.print_errors = 1;
  
  return cmdline_parser_config_file (filename, args_info, &params);
}

int
cmdline_parser_config_file (char * const filename,
                           struct gengetopt_args_info *args_info,
                           struct cmdline_parser_params *params)
{
  int i, result;
  int my_argc = 1;
  char **my_argv_arg;
  char *additional_error;

  /* store the program name */
  cmd_line_list_tmp = (struct line_list *) malloc (sizeof (struct line_list));
  cmd_line_list_tmp->next = cmd_line_list;
  cmd_line_list = cmd_line_list_tmp;
  cmd_line_list->string_arg = gengetopt_strdup (CMDLINE_PARSER_PACKAGE);

  result = _cmdline_parser_configfile(filename, &my_argc);

  if (result != EXIT_FAILURE) {
    my_argv_arg = (char **) malloc((my_argc+1) * sizeof(char *));
    cmd_line_list_tmp = cmd_line_list;

    for (i = my_argc - 1; i >= 0; --i) {
      my_argv_arg[i] = cmd_line_list_tmp->string_arg;
      cmd_line_list_tmp = cmd_line_list_tmp->next;
    }

    my_argv_arg[my_argc] = 0;

    additional_error = (char *)malloc(strlen(filename) + strlen(ADDITIONAL_ERROR) + 1);
    strcpy (additional_error, ADDITIONAL_ERROR);
    strcat (additional_error, filename);
    result =
      cmdline_parser_internal (my_argc, my_argv_arg, args_info,
                              params,
                              additional_error);

    free (additional_error);
    free (my_argv_arg);
  }

  free_cmd_list();
  return result;
}
