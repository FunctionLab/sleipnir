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
#define CMDLINE_PARSER_PACKAGE "BNServer"
#endif

#ifndef CMDLINE_PARSER_VERSION
#define CMDLINE_PARSER_VERSION "1.0"
#endif

struct gengetopt_args_info
{
  char * database_arg;	/* Database directory (default='.').  */
  char * input_arg;	/* Context IDs and names.  */
  char * contexts_arg;	/* Context/gene mapping.  */
  char * backgrounds_arg;	/* Background connectivities for all genes.  */
  char * files_arg;	/* File directory (default='.').  */
  char * graphviz_arg;	/* Graphviz executable path (default='fdp').  */
  char * config_arg;	/* Command line config file (default='BNServer.ini').  */
  int verbosity_arg;	/* Message verbosity (default='5').  */
  char * networks_arg;	/* Bayes net directory (default='.').  */
  char * default_arg;	/* Bayes net for no context.  */
  int xdsl_flag;	/* Use XDSL files instead of DSL (default=on).  */
  int minimal_in_flag;	/* Read stored contexts and minimal Bayes nets (default=off).  */
  char * minimal_out_arg;	/* Store contexts and minimal Bayes nets.  */
  char * kegg_arg;	/* KEGG ontology.  */
  char * kegg_org_arg;	/* KEGG organism (default='HSA').  */
  char * go_onto_arg;	/* GO ontology.  */
  char * go_anno_arg;	/* GO annotations.  */
  int port_arg;	/* Server port (default='1234').  */
  int timeout_arg;	/* Server timeout (default='100').  */
  
  int help_given ;	/* Whether help was given.  */
  int version_given ;	/* Whether version was given.  */
  int database_given ;	/* Whether database was given.  */
  int input_given ;	/* Whether input was given.  */
  int contexts_given ;	/* Whether contexts was given.  */
  int backgrounds_given ;	/* Whether backgrounds was given.  */
  int files_given ;	/* Whether files was given.  */
  int graphviz_given ;	/* Whether graphviz was given.  */
  int config_given ;	/* Whether config was given.  */
  int verbosity_given ;	/* Whether verbosity was given.  */
  int networks_given ;	/* Whether networks was given.  */
  int default_given ;	/* Whether default was given.  */
  int xdsl_given ;	/* Whether xdsl was given.  */
  int minimal_in_given ;	/* Whether minimal_in was given.  */
  int minimal_out_given ;	/* Whether minimal_out was given.  */
  int kegg_given ;	/* Whether kegg was given.  */
  int kegg_org_given ;	/* Whether kegg_org was given.  */
  int go_onto_given ;	/* Whether go_onto was given.  */
  int go_anno_given ;	/* Whether go_anno was given.  */
  int port_given ;	/* Whether port was given.  */
  int timeout_given ;	/* Whether timeout was given.  */

} ;

int cmdline_parser (int argc, char * const *argv, struct gengetopt_args_info *args_info);
int cmdline_parser2 (int argc, char * const *argv, struct gengetopt_args_info *args_info, int override, int initialize, int check_required);

void cmdline_parser_print_help(void);
void cmdline_parser_print_version(void);

void cmdline_parser_init (struct gengetopt_args_info *args_info);
void cmdline_parser_free (struct gengetopt_args_info *args_info);

int cmdline_parser_configfile (char * const filename, struct gengetopt_args_info *args_info, int override, int initialize, int check_required);

int cmdline_parser_required (struct gengetopt_args_info *args_info, const char *prog_name);


#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /* CMDLINE_H */
