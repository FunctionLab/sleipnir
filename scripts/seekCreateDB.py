# This script is an adaptation of the prepare_seek.py script that 
# comes with the breast cancer example dataset on the Seek website

import os
import re
import sys
import resource
import argparse

def read_dataset_list(n):
    f = open(n)
    c = 0
    m = []
    for l in f:
        c += 1
        l = l.rstrip("\n").split("\t")
        if len(l)!=3:
            sys.stdout.write("Error: does not contain 3 columns!\n");
            sys.stdout.write("Error line number: %d\n" % c);
            break
        else:
            m.append(tuple(l))
    f.close()
    return m

def write_dataset_list(dset_list, n):
    fw = open(n, "w")
    for (file_name, name, platform) in dset_list:
        m = re.match("(.*).pcl", file_name)
        datasetName = m.group(1)
        fw.write("%s\n" % datasetName)
    fw.close()

def write_dataset_platform_map(dset_list, n):
    fw = open(n, "w")
    for (file_name, name, platform) in dset_list:
        # m = re.match("(.*).pcl", file_name)
        # datasetName = m.group(1)
        fw.write("%s\t%s\n" % (name, platform))
    fw.close()


def write_db_list(n, numDBFiles):
    fw = open(n, "w")
    for i in range(0, numDBFiles):
        a_str = '{:08d}'.format(i)
        fw.write("../db/%s.db\n" % a_str)
    fw.close()


def createSeekDB(sleipnirBinDir, inputDatasetFile, pclDir, refDir, output_dir, numDBFiles=1000 ):
    pclDir = os.path.abspath(pclDir)
    output_dir = os.path.abspath(output_dir)
    sleipnirBinDir = os.path.abspath(sleipnirBinDir)
    refDir = os.path.abspath(refDir)

    if not os.path.exists(os.path.join(sleipnirBinDir, "Distancer")):
        print("Error: Distancer binary not found")
        return False
    if not os.path.exists(os.path.join(sleipnirBinDir, "SeekPrep")):
        print("Error: SeekPrep binary not found")
        return False
    if not os.path.exists(os.path.join(sleipnirBinDir, "Data2DB")):
        print("Error: Data2DB binary not found")
        return False
    if not os.path.exists(os.path.join(sleipnirBinDir, "PCL2Bin")):
        print("Error: PCL2Bin binary not found")
        return False
    
    # check max open files setting is sufficient, i.e. ulimit -n
    softFileLimit, _ = resource.getrlimit(resource.RLIMIT_NOFILE)
    if softFileLimit < numDBFiles:
        print("Max open files limit is insufficient: {}".format(softFileLimit))
        return False

    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)

    os.system("cp %s/quant2 %s/." % (refDir, output_dir))
    os.system("cp %s/gene_map.txt %s/." % (refDir, output_dir))
    if not os.path.exists(f'{output_dir}/quant2'):
        print("Missing quant file")
        return False
    if not os.path.exists(f'{output_dir}/gene_map.txt'):
        print("Missing gene_map file")
        return False

    datasets = read_dataset_list(inputDatasetFile)
    if len(datasets)==0:
        sys.stderr.write("Empty dataset list: %s!\n" % inputDatasetFile)
        sys.exit(1)

    if not os.path.isdir(pclDir):
        sys.stderr.write("PCL directory not found %s!\n" % pclDir)
        sys.exit(1)

    os.chdir(pclDir)
    for (dataset, name, platform) in datasets:
        if not os.path.exists(dataset):
            sys.stderr.write("Error, dataset %s does not exist!\n" % dataset)
            sys.exit(1)
        if not dataset.endswith(".pcl"):
            sys.stderr.write("Error, file name (%s) does not end in .pcl\n" % dataset)
            sys.exit(1)

    sys.stdout.write("Calculating correlation matrix...\n")
    if not os.path.exists(output_dir+"/dab"):
        os.mkdir(output_dir + "/dab")
    os.chdir(output_dir + "/dab")
    for (dataset, name, platform) in datasets:
        m = re.match("(.*).pcl", dataset)
        datasetName = m.group(1)
        cmd = "%s/Distancer -i %s/%s.pcl -o %s.dab -s 0 -t 1" % (sleipnirBinDir, pclDir, datasetName, datasetName)
        sys.stdout.write("Distancer: Processing %s...\n" % dataset)
        print(cmd)
        os.system(cmd)
        if os.path.exists("%s.quant" % datasetName):
            os.system("rm -rf %s.quant" % datasetName)
        cmd = "ln -s %s/quant2 %s.quant" % (output_dir, datasetName)
        print(cmd)
        os.system(cmd)
    sys.stdout.write("Done!\n")

    sys.stdout.write("Calculating gene z-score average...\n")
    if not os.path.exists(output_dir+"/prep"):
        os.mkdir(output_dir + "/prep")
    os.chdir(output_dir + "/prep")    
    for (dataset, name, platform) in datasets:
        m = re.match("(.*).pcl", dataset)
        datasetName = m.group(1)
        sys.stdout.write("SeekPrep Processing %s...\n" % datasetName)
        cmd = "%s/SeekPrep -d -a -B ../dab/%s.dab -i %s/gene_map.txt -D ." % (sleipnirBinDir, datasetName, output_dir)
        print(cmd)
        os.system(cmd)
        cmd = "%s/SeekPrep -d -p -B ../dab/%s.dab -i %s/gene_map.txt -D ." % (sleipnirBinDir, datasetName, output_dir)
        print(cmd)
        os.system(cmd)
    sys.stdout.write("Done!\n")

    sys.stdout.write("Joining correlation matrices followed by splitting by genes...\n")
    write_dataset_list(datasets, output_dir + "/dataset.list")
    if not os.path.exists(output_dir+"/db"):
        os.mkdir(output_dir + "/db")
    os.chdir(output_dir + "/db")
    cmd = "%s/Data2DB -i %s/gene_map.txt -D . -u -B 50  -f %d -x %s -d ../dab" % (sleipnirBinDir, output_dir, numDBFiles, output_dir + "/dataset.list")
    print(cmd)
    os.system(cmd)
    sys.stdout.write("Done!\n")

    sys.stdout.write("Calculating platform-wide gene score average...\n")
    if not os.path.exists(output_dir + "/plat"):
        os.mkdir(output_dir + "/plat")
    os.chdir(output_dir + "/plat")
    write_dataset_platform_map(datasets, output_dir + "/dataset_platform")
    write_db_list(output_dir + "/db_list", numDBFiles)
    cmd = "%s/SeekPrep -i %s/gene_map.txt -D . -f -P -b ../db_list -I ../prep -A ../dataset_platform -Q ../quant2" % (sleipnirBinDir, output_dir)
    print(cmd)
    os.system(cmd)
    sys.stdout.write("Done!\n")

    sys.stdout.write("Calculating sinfo files...\n")
    if not os.path.exists(output_dir + "/pclbin"):
        os.mkdir(output_dir + "/pclbin")
    if not os.path.exists(output_dir + "/sinfo"):
        os.mkdir(output_dir + "/sinfo")
    os.chdir(output_dir + "/sinfo")
    for (dataset, name, platform) in datasets:
        m = re.match("(.*).pcl", dataset)
        datasetName = m.group(1)
        sys.stdout.write("PCL2Bin Processing %s...\n" % dataset)
        cmd = "%s/PCL2Bin -i %s/%s.pcl -o ../pclbin/%s.pcl.bin" % (sleipnirBinDir, pclDir, datasetName, datasetName)
        print(cmd)
        os.system(cmd)
        cmd = "%s/SeekPrep -i %s/gene_map.txt -D . -e -V ../pclbin/%s.pcl.bin -s" % (sleipnirBinDir, output_dir, datasetName)
        print(cmd)
        os.system(cmd)
    sys.stdout.write("Done!\n")

    # TODO - create dataset size file
    sys.stdout.write("Creating dataset size...\n")
    os.chdir(output_dir)
    os.system("echo '' > ./dataset_size")
    for (dataset, name, platform) in datasets:
        cmd = "%s/SeekPrep -e -S -i /tmp -D ./ -V ./pclbin/%s.pcl.bin >> ./dataset_size" % (sleipnirBinDir, datasetName)
    sys.stdout.write("Done!\n")
    #output directory
    #directories to be created: dab, db, prep, sinfo
    return True

#set up script for SEEK
if __name__=="__main__":
    #input directory (containing the PCL's)
    #input file: dataset list (3-column format: file, name, platform)
    #input directory (containing the setting file: quant2, gene_map.txt)
    #path to sleipnir build
    argParser = argparse.ArgumentParser()
    argParser.add_argument('--pclDir', '-p', type=str, required=True,
                           help='Directory containing the PCL files for the new datasets')
    argParser.add_argument('--datasetFile', '-d', type=str, required=True,
                           help='Text file listing the datasets, one dataset per line, three columns (pcl_file, name, platform)')
    argParser.add_argument('--sleipnirBinDir', '-s', type=str, required=True,
                           help='Directory where the Sleipnir binaries are installed')
    argParser.add_argument('--refDir', '-r', type=str, required=True,
                           help='Directory of existing reference files (i.e. for reference gene_list and quant files)')
    argParser.add_argument('--outDir', '-o', type=str, required=True,
                           help='Output directory to write new database into')
    argParser.add_argument('--numDBFiles', '-n', type=int, required=False, default=1000,
                           help='Number of output DB files to spread gene data across (should match refDB number)')
    args = argParser.parse_args()

    createSeekDB(sleipnirBinDir = args.sleipnirBinDir,
                 inputDatasetFile = args.datasetFile,
                 pclDir = args.pclDir,
                 refDir = args.refDir,
                 output_dir = args.outDir,
                 numDBFiles = args.numDBFiles)

