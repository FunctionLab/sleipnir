"""
This utility parses refinebio metadata to extract the platform
type for each dataset. The chip_plat_map is a standard mapping
from chip type to platform type.
Example download command from refinebio:
    https://api.refine.bio/v1/experiments/GSE24528/
The chip_plat_map is derived from:
    "GEO IMPORT"
    Masters Thesis in Bioinformatics, Faisal Ibne Rezwan
    Department of Mathematical Sciences
    Chalmers University of Technology Sweden, January 2007
"""
import os
import re
import json
import argparse

chip_plat_map = {
    'HG_U95Av2': 'GPL91',
    'ATH1-121501': 'GPL198',
    'MG_U74A': 'GPL32',
    'MG_U74B': 'GPL33',
    'MG_U74C': 'GPL34',
    'MG_U74Av2': 'GPL81',
    'MG_U74Bv2': 'GPL82',
    'MG_U74Cv2': 'GPL83',
    'HG_U95A': 'GPL91',
    'HG_U95B': 'GPL92',
    'HG_U95C': 'GPL93',
    'HG_U95D': 'GPL94',
    'HG_U95E': 'GPL95',
    'Mu11KsubA': 'GPL175',
    'Mu11KsubB': 'GPL176',
    'AG1': 'GPL71',
    'DrosGenome1': 'GPL72',
    'RG_U34A': 'GPL85',
    'RG_U34B': 'GPL86',
    'RG_U34C': 'GPL87',
    'RT_U34': 'GPL89',
    'RN_U34': 'GPL88',
    'MOE430A': 'GPL339',
    'MOE430B': 'GPL340',
    'RAE230A': 'GPL341',
    'RAE230B': 'GPL342',
    'YG_S98': 'GPL90',
    'Ecoli': 'GPL73',
    'Ecoli_ASv2': 'GPL199',
    'Pae_G1a': 'GPL84',
    'Barley1': 'GPL1340',
    'HuGeneFL': 'GPL80',
    'HG-U133A': 'GPL96',
    'HG-U133B': 'GPL97',
    'Drosophila_2': 'GPL1322',
    'Mouse430A_2': 'GPL339',
    'HG-U133A_2': 'GPL571',
    'Zebrafish': 'GPL1319',
    'HG-Focus': 'GPL201',
    'Rat230_2': 'GPL1355',
    'HG-U133_Plus_2': 'GPL570',
    'Mouse430_2': 'GPL1261',
    'Yeast_2': 'GPL2529',
    'U133_X3P': 'GPL1352',
    'E_coli_2': 'GPL3154',
    'wheat': 'GPL3802',
    'Celegans': 'GPL200',
    'Xenopus_laevis': 'GPL1318',
    'Mapping10K_Xba131': 'GPL1266',
    'Mapping10K_Xba142': 'GPL2641',
    'HT_HG-U133A': 'GPL3921',
    'HT_HG-U133_Plus_PM': 'GPL13158',
    'Illumina_HumanHT-12_V3.0': ' GPL6947',
    'U133AAofAv2':  'GPL4685',
    'ZebGene-1_1-st': 'GPL18378',
    'Illumina_MouseWG-6_v2.0': 'GPL6887',
    'Illumina_Mouse-6_v1.1': 'GPL6105',
    'Illumina_MouseRef-8_V1.1': 'GPL6103',
    'MoGene-1_0-st': 'GPL6246',
    'MoGene-2_0-st': 'GPL16570',

}
# TODO - these are potential chipsets to add to map above
"""
Clariom_S_Mouse {GPL23038, GPL26644}
DroGene-1_0-st {GPL18094, GPL20842}
DroGene-1_1-st {GPL20767, GPL17506}
EleGene-1_0-st {GPL19230}
Hs_PromPR {GPL5082}
HT_MG-430A {GPL8759}
HT_MG-430_PM
HuEx-1_0-st
HuGene-1_0-st
HuGene-1_1-st
HuGene-2_0-st
miRNA-3
Mm_PromPR
MoEx-1_0-st
MoGene-1_0-st {many chose most samples}
MoGene-1_1-st
MoGene-2_0-st {many chose most samples}
MoGene-2_1-st
MOUSEDIVm520650
PrimeView
ZebGene-1_0-st
# Up to here the also have GPL Ids in refine_bio.json,
# below they do not
Illumina_MouseWG-6_v2.0 {GPL6887, GPL17543}
Illumina_Mouse-6_v1.1 {GPL6105, GPL6481}
Illumina_MouseRef-8_V1.1 {GPL6103}
"""

def getGsePlatform(metadata):
    platform = metadata['annotations'][0]['data']['platform_id'][0]
    return platform

def getE_Platform(platformName):
    # platformName = metadata['samples'][0]['platform_name']
    result = re.search('\[(.*)\]', platformName)
    if result is not None:
        chip = result.group(1)
    else:
        chip = platformName
    platform = chip_plat_map.get(chip)
    if platform is None:
        print(f'No platformId map for \'{chip}\'')
    return platform


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--dir', '-d', default=None,
                        help="directory of refine.bio metadata files to parse")
    args = parser.parse_args()

    if args.dir is None:
        raise Exception("Must specify -d <directory> to parse")

    for file in os.listdir(args.dir):
        fullpath = os.path.join(args.dir, file)
        # print(fullpath)
        if os.path.isfile(fullpath):
            with open(fullpath) as fp:
                platform = 'Default'
                try:
                    metadata = json.load(fp)
                except Exception as err:
                    print(f"Error: Json load: {file}")
                    continue
                try:
                    if file.startswith("GSE"):
                        platform = getGsePlatform(metadata)
                    elif file.startswith("E-"):
                        platName = metadata['samples'][0]['platform_name']
                        platform = getE_Platform(platName)
                    else:
                        print(f'Error: Filename {file} doesnt match expected pattern')
                except Exception as err:
                    print(f'Error {file}: {err}')
                    continue
                print(f'{file}    {platform}')


