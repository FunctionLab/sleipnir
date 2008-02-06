/* cmdline.h */

/* File autogenerated by gengetopt version 2.13.1  */

#ifndef CMDLINE_H
#define CMDLINE_H

/* If we use autoconf.  */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#ifndef CMDLINE_PARSER_PACKAGE
#define CMDLINE_PARSER_PACKAGE "Data2Sql"
#endif

#ifndef CMDLINE_PARSER_VERSION
#define CMDLINE_PARSER_VERSION "1.0"
#endif

struct gengetopt_args_info
{
  char * input_arg;	/* Input gene mapping.  */
  int verbosity_arg;	/* Message verbosity (default='5').  */
  int memmap_flag;	/* Memory map input/output (default=off).  */
  int datasets_flag;	/* Output datasets table (default=off).  */
  int block_arg;	/* Block size for SQL chunking (default='1000').  */
  char * host_arg;	/* Database hostname (default='localhost').  */
  char * port_arg;	/* Database port (default='3306').  */
  char * database_arg;	/* Database name.  */
  char * table_arg;	/* Database table name (default='datapairs').  */
  char * username_arg;	/* Database username (default='root').  */
  char * password_arg;	/* Database password.  */
  
  int help_given ;	/* Whether help was given.  */
  int version_given ;	/* Whether version was given.  */
  int input_given ;	/* Whether input was given.  */
  int verbosity_given ;	/* Whether verbosity was given.  */
  int memmap_given ;	/* Whether memmap was given.  */
  int datasets_given ;	/* Whether datasets was given.  */
  int block_given ;	/* Whether block was given.  */
  int host_given ;	/* Whether host was given.  */
  int port_given ;	/* Whether port was given.  */
  int database_given ;	/* Whether database was given.  */
  int table_given ;	/* Whether table was given.  */
  int username_given ;	/* Whether username was given.  */
  int password_given ;	/* Whether password was given.  */

  char **inputs ; /* unamed options */
  unsigned inputs_num ; /* unamed options number */
} ;

int cmdline_parser (int argc, char * const *argv, struct gengetopt_args_info *args_info);
int cmdline_parser2 (int argc, char * const *argv, struct gengetopt_args_info *args_info, int override, int initialize, int check_required);

void cmdline_parser_print_help(void);
void cmdline_parser_print_version(void);

void cmdline_parser_init (struct gengetopt_args_info *args_info);
void cmdline_parser_free (struct gengetopt_args_info *args_info);

int cmdline_parser_required (struct gengetopt_args_info *args_info, const char *prog_name);


#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /* CMDLINE_H */
