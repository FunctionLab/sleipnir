/*
  File autogenerated by gengetopt version 2.22
  generated with the following command:
  /memex/qzhu/usr/bin/gengetopt -iDBCombiner.ggo --default-optional -u -N -e 

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

const char *gengetopt_args_info_purpose = "Combines a list of DB files with the same gene content";

const char *gengetopt_args_info_usage = "Usage: DBCombiner [OPTIONS]... [FILES]...";

const char *gengetopt_args_info_description = "";

const char *gengetopt_args_info_help[] = {
        "  -h, --help                   Print help and exit",
        "  -V, --version                Print version and exit",
        "\nMode:",
        "  -C, --combine                Combine a set of DB's, each coming from a \n                                 different dataset subset  (default=off)",
        "  -R, --reorganize             Reorganize a set of DB's, such as from 21000 DB \n                                 files to 1000 DB files, ie expanding/shrinking \n                                 the number of genes a DB contains  \n                                 (default=off)",
        "\nMain:",
        "  -i, --input=filename         Input gene mapping",
        "\nCombine Mode:",
        "  -x, --db=filename            Input a set of databaselet filenames (including \n                                 path)",
        "  -D, --dir_out=directory      Output database directory  (default=`.')",
        "  -N, --is_nibble              Whether the input DB is nibble type  \n                                 (default=off)",
        "  -s, --split                  Split to one-gene per file  (default=off)",
        "\nReorganize Mode:",
        "  -A, --dataset=filename       Dataset-platform mapping file",
        "  -d, --db_dir=directory       Source DB collection directory",
        "  -n, --src_db_num=INT         Source DB number of files",
        "  -b, --dest_db_num=INT        Destination DB number of files",
        "  -B, --dest_db_dir=directory  Destination DB directory",
        0
};

typedef enum {
    ARG_NO, ARG_FLAG, ARG_STRING, ARG_INT
} cmdline_parser_arg_type;

static
void clear_given(struct gengetopt_args_info *args_info);

static
void clear_args(struct gengetopt_args_info *args_info);

static int
cmdline_parser_internal(int argc, char *const *argv, struct gengetopt_args_info *args_info,
                        struct cmdline_parser_params *params, const char *additional_error);

static int
cmdline_parser_required2(struct gengetopt_args_info *args_info, const char *prog_name, const char *additional_error);

static char *
gengetopt_strdup(const char *s);

static
void clear_given(struct gengetopt_args_info *args_info) {
    args_info->help_given = 0;
    args_info->version_given = 0;
    args_info->combine_given = 0;
    args_info->reorganize_given = 0;
    args_info->input_given = 0;
    args_info->db_given = 0;
    args_info->dir_out_given = 0;
    args_info->is_nibble_given = 0;
    args_info->split_given = 0;
    args_info->dataset_given = 0;
    args_info->db_dir_given = 0;
    args_info->src_db_num_given = 0;
    args_info->dest_db_num_given = 0;
    args_info->dest_db_dir_given = 0;
}

static
void clear_args(struct gengetopt_args_info *args_info) {
    args_info->combine_flag = 0;
    args_info->reorganize_flag = 0;
    args_info->input_arg = NULL;
    args_info->input_orig = NULL;
    args_info->db_arg = NULL;
    args_info->db_orig = NULL;
    args_info->dir_out_arg = gengetopt_strdup(".");
    args_info->dir_out_orig = NULL;
    args_info->is_nibble_flag = 0;
    args_info->split_flag = 0;
    args_info->dataset_arg = NULL;
    args_info->dataset_orig = NULL;
    args_info->db_dir_arg = NULL;
    args_info->db_dir_orig = NULL;
    args_info->src_db_num_orig = NULL;
    args_info->dest_db_num_orig = NULL;
    args_info->dest_db_dir_arg = NULL;
    args_info->dest_db_dir_orig = NULL;

}

static
void init_args_info(struct gengetopt_args_info *args_info) {


    args_info->help_help = gengetopt_args_info_help[0];
    args_info->version_help = gengetopt_args_info_help[1];
    args_info->combine_help = gengetopt_args_info_help[3];
    args_info->reorganize_help = gengetopt_args_info_help[4];
    args_info->input_help = gengetopt_args_info_help[6];
    args_info->db_help = gengetopt_args_info_help[8];
    args_info->dir_out_help = gengetopt_args_info_help[9];
    args_info->is_nibble_help = gengetopt_args_info_help[10];
    args_info->split_help = gengetopt_args_info_help[11];
    args_info->dataset_help = gengetopt_args_info_help[13];
    args_info->db_dir_help = gengetopt_args_info_help[14];
    args_info->src_db_num_help = gengetopt_args_info_help[15];
    args_info->dest_db_num_help = gengetopt_args_info_help[16];
    args_info->dest_db_dir_help = gengetopt_args_info_help[17];

}

void
cmdline_parser_print_version(void) {
    printf("%s %s\n", CMDLINE_PARSER_PACKAGE, CMDLINE_PARSER_VERSION);
}

static void print_help_common(void) {
    cmdline_parser_print_version();

    if (strlen(gengetopt_args_info_purpose) > 0)
        printf("\n%s\n", gengetopt_args_info_purpose);

    if (strlen(gengetopt_args_info_usage) > 0)
        printf("\n%s\n", gengetopt_args_info_usage);

    printf("\n");

    if (strlen(gengetopt_args_info_description) > 0)
        printf("%s\n", gengetopt_args_info_description);
}

void
cmdline_parser_print_help(void) {
    int i = 0;
    print_help_common();
    while (gengetopt_args_info_help[i])
        printf("%s\n", gengetopt_args_info_help[i++]);
}

void
cmdline_parser_init(struct gengetopt_args_info *args_info) {
    clear_given(args_info);
    clear_args(args_info);
    init_args_info(args_info);

    args_info->inputs = NULL;
    args_info->inputs_num = 0;
}

void
cmdline_parser_params_init(struct cmdline_parser_params *params) {
    if (params) {
        params->override = 0;
        params->initialize = 1;
        params->check_required = 1;
        params->check_ambiguity = 0;
        params->print_errors = 1;
    }
}

struct cmdline_parser_params *
cmdline_parser_params_create(void) {
    struct cmdline_parser_params *params =
            (struct cmdline_parser_params *) malloc(sizeof(struct cmdline_parser_params));
    cmdline_parser_params_init(params);
    return params;
}

static void
free_string_field(char **s) {
    if (*s) {
        free(*s);
        *s = 0;
    }
}


static void
cmdline_parser_release(struct gengetopt_args_info *args_info) {
    unsigned int i;
    free_string_field(&(args_info->input_arg));
    free_string_field(&(args_info->input_orig));
    free_string_field(&(args_info->db_arg));
    free_string_field(&(args_info->db_orig));
    free_string_field(&(args_info->dir_out_arg));
    free_string_field(&(args_info->dir_out_orig));
    free_string_field(&(args_info->dataset_arg));
    free_string_field(&(args_info->dataset_orig));
    free_string_field(&(args_info->db_dir_arg));
    free_string_field(&(args_info->db_dir_orig));
    free_string_field(&(args_info->src_db_num_orig));
    free_string_field(&(args_info->dest_db_num_orig));
    free_string_field(&(args_info->dest_db_dir_arg));
    free_string_field(&(args_info->dest_db_dir_orig));


    for (i = 0; i < args_info->inputs_num; ++i)
        free(args_info->inputs[i]);

    if (args_info->inputs_num)
        free(args_info->inputs);

    clear_given(args_info);
}


static void
write_into_file(FILE *outfile, const char *opt, const char *arg, char *values[]) {
    if (arg) {
        fprintf(outfile, "%s=\"%s\"\n", opt, arg);
    } else {
        fprintf(outfile, "%s\n", opt);
    }
}


int
cmdline_parser_dump(FILE *outfile, struct gengetopt_args_info *args_info) {
    int i = 0;

    if (!outfile) {
        fprintf(stderr, "%s: cannot dump options to stream\n", CMDLINE_PARSER_PACKAGE);
        return EXIT_FAILURE;
    }

    if (args_info->help_given)
        write_into_file(outfile, "help", 0, 0);
    if (args_info->version_given)
        write_into_file(outfile, "version", 0, 0);
    if (args_info->combine_given)
        write_into_file(outfile, "combine", 0, 0);
    if (args_info->reorganize_given)
        write_into_file(outfile, "reorganize", 0, 0);
    if (args_info->input_given)
        write_into_file(outfile, "input", args_info->input_orig, 0);
    if (args_info->db_given)
        write_into_file(outfile, "db", args_info->db_orig, 0);
    if (args_info->dir_out_given)
        write_into_file(outfile, "dir_out", args_info->dir_out_orig, 0);
    if (args_info->is_nibble_given)
        write_into_file(outfile, "is_nibble", 0, 0);
    if (args_info->split_given)
        write_into_file(outfile, "split", 0, 0);
    if (args_info->dataset_given)
        write_into_file(outfile, "dataset", args_info->dataset_orig, 0);
    if (args_info->db_dir_given)
        write_into_file(outfile, "db_dir", args_info->db_dir_orig, 0);
    if (args_info->src_db_num_given)
        write_into_file(outfile, "src_db_num", args_info->src_db_num_orig, 0);
    if (args_info->dest_db_num_given)
        write_into_file(outfile, "dest_db_num", args_info->dest_db_num_orig, 0);
    if (args_info->dest_db_dir_given)
        write_into_file(outfile, "dest_db_dir", args_info->dest_db_dir_orig, 0);


    i = EXIT_SUCCESS;
    return i;
}

int
cmdline_parser_file_save(const char *filename, struct gengetopt_args_info *args_info) {
    FILE *outfile;
    int i = 0;

    outfile = fopen(filename, "w");

    if (!outfile) {
        fprintf(stderr, "%s: cannot open file for writing: %s\n", CMDLINE_PARSER_PACKAGE, filename);
        return EXIT_FAILURE;
    }

    i = cmdline_parser_dump(outfile, args_info);
    fclose(outfile);

    return i;
}

void
cmdline_parser_free(struct gengetopt_args_info *args_info) {
    cmdline_parser_release(args_info);
}

/** @brief replacement of strdup, which is not standard */
char *
gengetopt_strdup(const char *s) {
    char *result = NULL;
    if (!s)
        return result;

    result = (char *) malloc(strlen(s) + 1);
    if (result == (char *) 0)
        return (char *) 0;
    strcpy(result, s);
    return result;
}

int
cmdline_parser(int argc, char *const *argv, struct gengetopt_args_info *args_info) {
    return cmdline_parser2(argc, argv, args_info, 0, 1, 1);
}

int
cmdline_parser_ext(int argc, char *const *argv, struct gengetopt_args_info *args_info,
                   struct cmdline_parser_params *params) {
    int result;
    result = cmdline_parser_internal(argc, argv, args_info, params, NULL);

    return result;
}

int
cmdline_parser2(int argc, char *const *argv, struct gengetopt_args_info *args_info, int override, int initialize,
                int check_required) {
    int result;
    struct cmdline_parser_params params;

    params.override = override;
    params.initialize = initialize;
    params.check_required = check_required;
    params.check_ambiguity = 0;
    params.print_errors = 1;

    result = cmdline_parser_internal(argc, argv, args_info, &params, NULL);

    return result;
}

int
cmdline_parser_required(struct gengetopt_args_info *args_info, const char *prog_name) {
    int result = EXIT_SUCCESS;

    if (cmdline_parser_required2(args_info, prog_name, NULL) > 0)
        result = EXIT_FAILURE;

    return result;
}

int
cmdline_parser_required2(struct gengetopt_args_info *args_info, const char *prog_name, const char *additional_error) {
    int error = 0;

    /* checks for required options */
    if (!args_info->input_given) {
        fprintf(stderr, "%s: '--input' ('-i') option required%s\n", prog_name,
                (additional_error ? additional_error : ""));
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
               const char *additional_error) {
    char *stop_char = 0;
    const char *val = value;
    int found;
    char **string_field;

    stop_char = 0;
    found = 0;

    if (!multiple_option && prev_given && (*prev_given || (check_ambiguity && *field_given))) {
        if (short_opt != '-')
            fprintf(stderr, "%s: `--%s' (`-%c') option given more than once%s\n",
                    package_name, long_opt, short_opt,
                    (additional_error ? additional_error : ""));
        else
            fprintf(stderr, "%s: `--%s' option given more than once%s\n",
                    package_name, long_opt,
                    (additional_error ? additional_error : ""));
        return 1; /* failure */
    }


    if (field_given && *field_given && !override)
        return 0;
    if (prev_given)
        (*prev_given)++;
    if (field_given)
        (*field_given)++;
    if (possible_values)
        val = possible_values[found];

    switch (arg_type) {
        case ARG_FLAG:
            *((int *) field) = !*((int *) field);
            break;
        case ARG_INT:
            if (val) *((int *) field) = strtol(val, &stop_char, 0);
            break;
        case ARG_STRING:
            if (val) {
                string_field = (char **) field;
                if (!no_free && *string_field)
                    free(*string_field); /* free previous string */
                *string_field = gengetopt_strdup(val);
            }
            break;
        default:
            break;
    };

    /* check numeric conversion */
    switch (arg_type) {
        case ARG_INT:
            if (val && !(stop_char && *stop_char == '\0')) {
                fprintf(stderr, "%s: invalid numeric value: %s\n", package_name, val);
                return 1; /* failure */
            }
            break;
        default:;
    };

    /* store the original value */
    switch (arg_type) {
        case ARG_NO:
        case ARG_FLAG:
            break;
        default:
            if (value && orig_field) {
                if (no_free) {
                    *orig_field = value;
                } else {
                    if (*orig_field)
                        free(*orig_field); /* free previous string */
                    *orig_field = gengetopt_strdup(value);
                }
            }
    };

    return 0; /* OK */
}


int
cmdline_parser_internal(int argc, char *const *argv, struct gengetopt_args_info *args_info,
                        struct cmdline_parser_params *params, const char *additional_error) {
    int c;    /* Character of the parsed option.  */

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
        cmdline_parser_init(args_info);

    cmdline_parser_init(&local_args_info);

    optarg = 0;
    optind = 0;
    opterr = params->print_errors;
    optopt = '?';

    while (1) {
        int option_index = 0;

        static struct option long_options[] = {
                {"help",        0, NULL, 'h'},
                {"version",     0, NULL, 'V'},
                {"combine",     0, NULL, 'C'},
                {"reorganize",  0, NULL, 'R'},
                {"input",       1, NULL, 'i'},
                {"db",          1, NULL, 'x'},
                {"dir_out",     1, NULL, 'D'},
                {"is_nibble",   0, NULL, 'N'},
                {"split",       0, NULL, 's'},
                {"dataset",     1, NULL, 'A'},
                {"db_dir",      1, NULL, 'd'},
                {"src_db_num",  1, NULL, 'n'},
                {"dest_db_num", 1, NULL, 'b'},
                {"dest_db_dir", 1, NULL, 'B'},
                {NULL,          0, NULL, 0}
        };

        c = getopt_long(argc, argv, "hVCRi:x:D:NsA:d:n:b:B:", long_options, &option_index);

        if (c == -1) break;    /* Exit from `while (1)' loop.  */

        switch (c) {
            case 'h':    /* Print help and exit.  */
                cmdline_parser_print_help();
                cmdline_parser_free(&local_args_info);
                exit(EXIT_SUCCESS);

            case 'V':    /* Print version and exit.  */


                if (update_arg(0,
                               0, &(args_info->version_given),
                               &(local_args_info.version_given), optarg, 0, 0, ARG_NO,
                               check_ambiguity, override, 0, 0,
                               "version", 'V',
                               additional_error))
                    goto failure;
                cmdline_parser_free(&local_args_info);
                return 0;

                break;
            case 'C':    /* Combine a set of DB's, each coming from a different dataset subset.  */


                if (update_arg((void *) &(args_info->combine_flag), 0, &(args_info->combine_given),
                               &(local_args_info.combine_given), optarg, 0, 0, ARG_FLAG,
                               check_ambiguity, override, 1, 0, "combine", 'C',
                               additional_error))
                    goto failure;

                break;
            case 'R':    /* Reorganize a set of DB's, such as from 21000 DB files to 1000 DB files, ie expanding/shrinking the number of genes a DB contains.  */


                if (update_arg((void *) &(args_info->reorganize_flag), 0, &(args_info->reorganize_given),
                               &(local_args_info.reorganize_given), optarg, 0, 0, ARG_FLAG,
                               check_ambiguity, override, 1, 0, "reorganize", 'R',
                               additional_error))
                    goto failure;

                break;
            case 'i':    /* Input gene mapping.  */


                if (update_arg((void *) &(args_info->input_arg),
                               &(args_info->input_orig), &(args_info->input_given),
                               &(local_args_info.input_given), optarg, 0, 0, ARG_STRING,
                               check_ambiguity, override, 0, 0,
                               "input", 'i',
                               additional_error))
                    goto failure;

                break;
            case 'x':    /* Input a set of databaselet filenames (including path).  */


                if (update_arg((void *) &(args_info->db_arg),
                               &(args_info->db_orig), &(args_info->db_given),
                               &(local_args_info.db_given), optarg, 0, 0, ARG_STRING,
                               check_ambiguity, override, 0, 0,
                               "db", 'x',
                               additional_error))
                    goto failure;

                break;
            case 'D':    /* Output database directory.  */


                if (update_arg((void *) &(args_info->dir_out_arg),
                               &(args_info->dir_out_orig), &(args_info->dir_out_given),
                               &(local_args_info.dir_out_given), optarg, 0, ".", ARG_STRING,
                               check_ambiguity, override, 0, 0,
                               "dir_out", 'D',
                               additional_error))
                    goto failure;

                break;
            case 'N':    /* Whether the input DB is nibble type.  */


                if (update_arg((void *) &(args_info->is_nibble_flag), 0, &(args_info->is_nibble_given),
                               &(local_args_info.is_nibble_given), optarg, 0, 0, ARG_FLAG,
                               check_ambiguity, override, 1, 0, "is_nibble", 'N',
                               additional_error))
                    goto failure;

                break;
            case 's':    /* Split to one-gene per file.  */


                if (update_arg((void *) &(args_info->split_flag), 0, &(args_info->split_given),
                               &(local_args_info.split_given), optarg, 0, 0, ARG_FLAG,
                               check_ambiguity, override, 1, 0, "split", 's',
                               additional_error))
                    goto failure;

                break;
            case 'A':    /* Dataset-platform mapping file.  */


                if (update_arg((void *) &(args_info->dataset_arg),
                               &(args_info->dataset_orig), &(args_info->dataset_given),
                               &(local_args_info.dataset_given), optarg, 0, 0, ARG_STRING,
                               check_ambiguity, override, 0, 0,
                               "dataset", 'A',
                               additional_error))
                    goto failure;

                break;
            case 'd':    /* Source DB collection directory.  */


                if (update_arg((void *) &(args_info->db_dir_arg),
                               &(args_info->db_dir_orig), &(args_info->db_dir_given),
                               &(local_args_info.db_dir_given), optarg, 0, 0, ARG_STRING,
                               check_ambiguity, override, 0, 0,
                               "db_dir", 'd',
                               additional_error))
                    goto failure;

                break;
            case 'n':    /* Source DB number of files.  */


                if (update_arg((void *) &(args_info->src_db_num_arg),
                               &(args_info->src_db_num_orig), &(args_info->src_db_num_given),
                               &(local_args_info.src_db_num_given), optarg, 0, 0, ARG_INT,
                               check_ambiguity, override, 0, 0,
                               "src_db_num", 'n',
                               additional_error))
                    goto failure;

                break;
            case 'b':    /* Destination DB number of files.  */


                if (update_arg((void *) &(args_info->dest_db_num_arg),
                               &(args_info->dest_db_num_orig), &(args_info->dest_db_num_given),
                               &(local_args_info.dest_db_num_given), optarg, 0, 0, ARG_INT,
                               check_ambiguity, override, 0, 0,
                               "dest_db_num", 'b',
                               additional_error))
                    goto failure;

                break;
            case 'B':    /* Destination DB directory.  */


                if (update_arg((void *) &(args_info->dest_db_dir_arg),
                               &(args_info->dest_db_dir_orig), &(args_info->dest_db_dir_given),
                               &(local_args_info.dest_db_dir_given), optarg, 0, 0, ARG_STRING,
                               check_ambiguity, override, 0, 0,
                               "dest_db_dir", 'B',
                               additional_error))
                    goto failure;

                break;

            case 0:    /* Long option with no short option */
            case '?':    /* Invalid option.  */
                /* `getopt_long' already printed an error message.  */
                goto failure;

            default:    /* bug: option not considered.  */
                fprintf(stderr, "%s: option unknown: %c%s\n", CMDLINE_PARSER_PACKAGE, c,
                        (additional_error ? additional_error : ""));
                abort();
        } /* switch */
    } /* while */



    if (check_required) {
        error += cmdline_parser_required2(args_info, argv[0], additional_error);
    }

    cmdline_parser_release(&local_args_info);

    if (error)
        return (EXIT_FAILURE);

    if (optind < argc) {
        int i = 0;
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
                (char **) (malloc((args_info->inputs_num) * sizeof(char *)));
        while (optind < argc)
            if (argv[optind++] != argv[0])
                args_info->inputs[i++] = gengetopt_strdup(argv[optind - 1]);
    }

    return 0;

    failure:

    cmdline_parser_release(&local_args_info);
    return (EXIT_FAILURE);
}
