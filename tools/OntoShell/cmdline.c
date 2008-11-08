/*
  File autogenerated by gengetopt version 2.22
  generated with the following command:
  /r01/tergeo/chuttenh/sleipnir/trunk/../extlib/gengetopt-2.22/bin/gengetopt -iOntoShell.ggo --default-optional -C -N -e 

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

const char *gengetopt_args_info_purpose = "Ontology exploration utility";

const char *gengetopt_args_info_usage = "Usage: OntoShell [OPTIONS]...";

const char *gengetopt_args_info_description = "";

const char *gengetopt_args_info_help[] = {
  "  -h, --help                 Print help and exit",
  "  -V, --version              Print version and exit",
  "\nMain:",
  "  -x, --exec=STRING          Command to execute",
  "  -l, --altids               Force alternate (human) IDs  (default=off)",
  "  -b, --dbids                Include GO database IDs  (default=off)",
  "\nFunction Catalogs:",
  "  -o, --go_onto=filename     GO ontology",
  "  -g, --go_anno=filename     GO annotations",
  "  -k, --kegg=filename        KEGG ontology",
  "  -K, --kegg_org=STRING      KEGG organism  (default=`SCE')",
  "  -m, --mips_onto=filename   MIPS ontology",
  "  -a, --mips_anno=filename   MIPS annotations",
  "  -M, --mipsp_onto=filename  MIPS phenotypes ontology",
  "  -A, --mipsp_anno=filename  MIPS phenotypes annotations",
  "  -f, --features=filename    SGD gene features",
  "\nOptional:",
  "  -z, --zeroes               Tab-complete zero entries  (default=off)",
  "  -c, --config=filename      Command line config file  \n                               (default=`OntoShell.ini')",
  "  -v, --verbosity=INT        Message verbosity  (default=`5')",
    0
};

typedef enum {ARG_NO
  , ARG_FLAG
  , ARG_STRING
  , ARG_INT
} cmdline_parser_arg_type;

static
void clear_given (struct gengetopt_args_info *args_info);
static
void clear_args (struct gengetopt_args_info *args_info);

static int
cmdline_parser_internal (int argc, char * const *argv, struct gengetopt_args_info *args_info,
                        struct cmdline_parser_params *params, const char *additional_error);

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


static char *
gengetopt_strdup (const char *s);

static
void clear_given (struct gengetopt_args_info *args_info)
{
  args_info->help_given = 0 ;
  args_info->version_given = 0 ;
  args_info->exec_given = 0 ;
  args_info->altids_given = 0 ;
  args_info->dbids_given = 0 ;
  args_info->go_onto_given = 0 ;
  args_info->go_anno_given = 0 ;
  args_info->kegg_given = 0 ;
  args_info->kegg_org_given = 0 ;
  args_info->mips_onto_given = 0 ;
  args_info->mips_anno_given = 0 ;
  args_info->mipsp_onto_given = 0 ;
  args_info->mipsp_anno_given = 0 ;
  args_info->features_given = 0 ;
  args_info->zeroes_given = 0 ;
  args_info->config_given = 0 ;
  args_info->verbosity_given = 0 ;
}

static
void clear_args (struct gengetopt_args_info *args_info)
{
  args_info->exec_arg = NULL;
  args_info->exec_orig = NULL;
  args_info->altids_flag = 0;
  args_info->dbids_flag = 0;
  args_info->go_onto_arg = NULL;
  args_info->go_onto_orig = NULL;
  args_info->go_anno_arg = NULL;
  args_info->go_anno_orig = NULL;
  args_info->kegg_arg = NULL;
  args_info->kegg_orig = NULL;
  args_info->kegg_org_arg = gengetopt_strdup ("SCE");
  args_info->kegg_org_orig = NULL;
  args_info->mips_onto_arg = NULL;
  args_info->mips_onto_orig = NULL;
  args_info->mips_anno_arg = NULL;
  args_info->mips_anno_orig = NULL;
  args_info->mipsp_onto_arg = NULL;
  args_info->mipsp_onto_orig = NULL;
  args_info->mipsp_anno_arg = NULL;
  args_info->mipsp_anno_orig = NULL;
  args_info->features_arg = NULL;
  args_info->features_orig = NULL;
  args_info->zeroes_flag = 0;
  args_info->config_arg = gengetopt_strdup ("OntoShell.ini");
  args_info->config_orig = NULL;
  args_info->verbosity_arg = 5;
  args_info->verbosity_orig = NULL;
  
}

static
void init_args_info(struct gengetopt_args_info *args_info)
{


  args_info->help_help = gengetopt_args_info_help[0] ;
  args_info->version_help = gengetopt_args_info_help[1] ;
  args_info->exec_help = gengetopt_args_info_help[3] ;
  args_info->altids_help = gengetopt_args_info_help[4] ;
  args_info->dbids_help = gengetopt_args_info_help[5] ;
  args_info->go_onto_help = gengetopt_args_info_help[7] ;
  args_info->go_anno_help = gengetopt_args_info_help[8] ;
  args_info->kegg_help = gengetopt_args_info_help[9] ;
  args_info->kegg_org_help = gengetopt_args_info_help[10] ;
  args_info->mips_onto_help = gengetopt_args_info_help[11] ;
  args_info->mips_anno_help = gengetopt_args_info_help[12] ;
  args_info->mipsp_onto_help = gengetopt_args_info_help[13] ;
  args_info->mipsp_anno_help = gengetopt_args_info_help[14] ;
  args_info->features_help = gengetopt_args_info_help[15] ;
  args_info->zeroes_help = gengetopt_args_info_help[17] ;
  args_info->config_help = gengetopt_args_info_help[18] ;
  args_info->verbosity_help = gengetopt_args_info_help[19] ;
  
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

  free_string_field (&(args_info->exec_arg));
  free_string_field (&(args_info->exec_orig));
  free_string_field (&(args_info->go_onto_arg));
  free_string_field (&(args_info->go_onto_orig));
  free_string_field (&(args_info->go_anno_arg));
  free_string_field (&(args_info->go_anno_orig));
  free_string_field (&(args_info->kegg_arg));
  free_string_field (&(args_info->kegg_orig));
  free_string_field (&(args_info->kegg_org_arg));
  free_string_field (&(args_info->kegg_org_orig));
  free_string_field (&(args_info->mips_onto_arg));
  free_string_field (&(args_info->mips_onto_orig));
  free_string_field (&(args_info->mips_anno_arg));
  free_string_field (&(args_info->mips_anno_orig));
  free_string_field (&(args_info->mipsp_onto_arg));
  free_string_field (&(args_info->mipsp_onto_orig));
  free_string_field (&(args_info->mipsp_anno_arg));
  free_string_field (&(args_info->mipsp_anno_orig));
  free_string_field (&(args_info->features_arg));
  free_string_field (&(args_info->features_orig));
  free_string_field (&(args_info->config_arg));
  free_string_field (&(args_info->config_orig));
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
  if (args_info->exec_given)
    write_into_file(outfile, "exec", args_info->exec_orig, 0);
  if (args_info->altids_given)
    write_into_file(outfile, "altids", 0, 0 );
  if (args_info->dbids_given)
    write_into_file(outfile, "dbids", 0, 0 );
  if (args_info->go_onto_given)
    write_into_file(outfile, "go_onto", args_info->go_onto_orig, 0);
  if (args_info->go_anno_given)
    write_into_file(outfile, "go_anno", args_info->go_anno_orig, 0);
  if (args_info->kegg_given)
    write_into_file(outfile, "kegg", args_info->kegg_orig, 0);
  if (args_info->kegg_org_given)
    write_into_file(outfile, "kegg_org", args_info->kegg_org_orig, 0);
  if (args_info->mips_onto_given)
    write_into_file(outfile, "mips_onto", args_info->mips_onto_orig, 0);
  if (args_info->mips_anno_given)
    write_into_file(outfile, "mips_anno", args_info->mips_anno_orig, 0);
  if (args_info->mipsp_onto_given)
    write_into_file(outfile, "mipsp_onto", args_info->mipsp_onto_orig, 0);
  if (args_info->mipsp_anno_given)
    write_into_file(outfile, "mipsp_anno", args_info->mipsp_anno_orig, 0);
  if (args_info->features_given)
    write_into_file(outfile, "features", args_info->features_orig, 0);
  if (args_info->zeroes_given)
    write_into_file(outfile, "zeroes", 0, 0 );
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
        { "exec",	1, NULL, 'x' },
        { "altids",	0, NULL, 'l' },
        { "dbids",	0, NULL, 'b' },
        { "go_onto",	1, NULL, 'o' },
        { "go_anno",	1, NULL, 'g' },
        { "kegg",	1, NULL, 'k' },
        { "kegg_org",	1, NULL, 'K' },
        { "mips_onto",	1, NULL, 'm' },
        { "mips_anno",	1, NULL, 'a' },
        { "mipsp_onto",	1, NULL, 'M' },
        { "mipsp_anno",	1, NULL, 'A' },
        { "features",	1, NULL, 'f' },
        { "zeroes",	0, NULL, 'z' },
        { "config",	1, NULL, 'c' },
        { "verbosity",	1, NULL, 'v' },
        { NULL,	0, NULL, 0 }
      };

      c = getopt_long (argc, argv, "hVx:lbo:g:k:K:m:a:M:A:f:zc:v:", long_options, &option_index);

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
        case 'x':	/* Command to execute.  */
        
        
          if (update_arg( (void *)&(args_info->exec_arg), 
               &(args_info->exec_orig), &(args_info->exec_given),
              &(local_args_info.exec_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "exec", 'x',
              additional_error))
            goto failure;
        
          break;
        case 'l':	/* Force alternate (human) IDs.  */
        
        
          if (update_arg((void *)&(args_info->altids_flag), 0, &(args_info->altids_given),
              &(local_args_info.altids_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "altids", 'l',
              additional_error))
            goto failure;
        
          break;
        case 'b':	/* Include GO database IDs.  */
        
        
          if (update_arg((void *)&(args_info->dbids_flag), 0, &(args_info->dbids_given),
              &(local_args_info.dbids_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "dbids", 'b',
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
        case 'g':	/* GO annotations.  */
        
        
          if (update_arg( (void *)&(args_info->go_anno_arg), 
               &(args_info->go_anno_orig), &(args_info->go_anno_given),
              &(local_args_info.go_anno_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "go_anno", 'g',
              additional_error))
            goto failure;
        
          break;
        case 'k':	/* KEGG ontology.  */
        
        
          if (update_arg( (void *)&(args_info->kegg_arg), 
               &(args_info->kegg_orig), &(args_info->kegg_given),
              &(local_args_info.kegg_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "kegg", 'k',
              additional_error))
            goto failure;
        
          break;
        case 'K':	/* KEGG organism.  */
        
        
          if (update_arg( (void *)&(args_info->kegg_org_arg), 
               &(args_info->kegg_org_orig), &(args_info->kegg_org_given),
              &(local_args_info.kegg_org_given), optarg, 0, "SCE", ARG_STRING,
              check_ambiguity, override, 0, 0,
              "kegg_org", 'K',
              additional_error))
            goto failure;
        
          break;
        case 'm':	/* MIPS ontology.  */
        
        
          if (update_arg( (void *)&(args_info->mips_onto_arg), 
               &(args_info->mips_onto_orig), &(args_info->mips_onto_given),
              &(local_args_info.mips_onto_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "mips_onto", 'm',
              additional_error))
            goto failure;
        
          break;
        case 'a':	/* MIPS annotations.  */
        
        
          if (update_arg( (void *)&(args_info->mips_anno_arg), 
               &(args_info->mips_anno_orig), &(args_info->mips_anno_given),
              &(local_args_info.mips_anno_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "mips_anno", 'a',
              additional_error))
            goto failure;
        
          break;
        case 'M':	/* MIPS phenotypes ontology.  */
        
        
          if (update_arg( (void *)&(args_info->mipsp_onto_arg), 
               &(args_info->mipsp_onto_orig), &(args_info->mipsp_onto_given),
              &(local_args_info.mipsp_onto_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "mipsp_onto", 'M',
              additional_error))
            goto failure;
        
          break;
        case 'A':	/* MIPS phenotypes annotations.  */
        
        
          if (update_arg( (void *)&(args_info->mipsp_anno_arg), 
               &(args_info->mipsp_anno_orig), &(args_info->mipsp_anno_given),
              &(local_args_info.mipsp_anno_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "mipsp_anno", 'A',
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
        case 'z':	/* Tab-complete zero entries.  */
        
        
          if (update_arg((void *)&(args_info->zeroes_flag), 0, &(args_info->zeroes_given),
              &(local_args_info.zeroes_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "zeroes", 'z',
              additional_error))
            goto failure;
        
          break;
        case 'c':	/* Command line config file.  */
        
        
          if (update_arg( (void *)&(args_info->config_arg), 
               &(args_info->config_orig), &(args_info->config_given),
              &(local_args_info.config_given), optarg, 0, "OntoShell.ini", ARG_STRING,
              check_ambiguity, override, 0, 0,
              "config", 'c',
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
