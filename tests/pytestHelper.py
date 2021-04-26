import os
import glob
import subprocess

testDir = os.path.abspath(os.path.dirname(__file__))
testInputsDir = os.path.join(testDir, 'test_inputs')
testOutputsDir = os.path.join(testDir, 'test_outputs')
sampleBcDir = os.path.join(testOutputsDir, 'sampleBC')
sleipnirDir = os.path.dirname(testDir)
sleipnirBin = os.path.join(sleipnirDir, 'Debug')
seekScriptsDir = os.path.join(sleipnirDir, 'scripts', 'seek')

def createSampleDatabase():
    dbDir = os.path.join(sampleBcDir, 'db')
    dbFiles = glob.glob1(dbDir, '*.db')
    if len(dbFiles) < 100:
        os.makedirs(sampleBcDir, exist_ok=True)
        # Copy the reference files needed to create the breast cancer sample
        #  project from the test_inputs directory to the destination in test_outputs
        breastCancerRefDir = os.path.join(testInputsDir, 'breast-cancer-sample')
        os.system(f'cp -a {breastCancerRefDir}/* {sampleBcDir}/')
        # Create the sample DB
        cmd = f'python {seekScriptsDir}/seekCreateDB.py --all -g bc_gene_map.txt ' \
            f'-n 100 -m 4 -i {sampleBcDir} -o {sampleBcDir} -b {sleipnirBin} --dab-use-gene-set'
        ret = subprocess.run(cmd, shell=True)
        assert ret.returncode == 0
    return sampleBcDir