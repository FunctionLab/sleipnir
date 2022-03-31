import os
import re
import json
import argparse
from collections import Counter
from mapDsetToPlatform import getE_Platform

def perSpeciesDsetList(data, outdir):
    multiSpeciesList = []
    dsetMap = {}
    for key, val in data.items():
        # print(f'{key} -> {val}\n')
        speciesList = val['organism_names']
        if len(speciesList) > 1:
            multiSpeciesList.append(key)
        for species in speciesList:
            arr = dsetMap.get(species)
            if arr is None:
                arr = []
                dsetMap[species] = arr
            arr.append(key)

    print(f'Multi-species count: {len(multiSpeciesList)}')
    with open(os.path.join(outdir, 'multiSpeciesDsets.txt'), 'w') as fp:
        fp.write("\n".join(multiSpeciesList))

    for key, val in dsetMap.items():
        with open(os.path.join(outdir, f'{key}.txt'), 'w') as fp:
            fp.write("\n".join(val))
            fp.write("\n")


def getPlatformFromSamples(dsetName, item):
    samples = item.get("samples")
    if samples is not None and len(samples) > 0:
        platformNameList = []
        platformIdList = []
        if type(samples) is list:
            for sample in samples:
                platName = sample.get('platform_name')
                platId = sample.get('platform_id')
                if platName is not None:
                    platformNameList.append(platName)
                if platId is not None:
                    platformIdList.append(platId)
        elif type(samples) is dict:
            for key, val in samples.items():
                platName = val.get('platform_name')
                platId = val.get('platform_id')
                if platName is not None:
                    platformNameList.append(platName)
                if platId is not None:
                    platformIdList.append(platId)
        else:
            print(f"Samples section not of supported type, {type(samples)}")
            return None
        numPlatformNames = len(set(platformNameList)) # consolidate identical items
        numPlatformIds = len(set(platformIdList)) # consolidate identical items
        if numPlatformIds == 1:
            # check there are no spaces in the platformId
            if bool(re.search(r"\s", platformIdList[0])) is False:
                return platformIdList[0]
            else:
                print(f"Space found in platformId {platformIdList[0]}")
        if numPlatformNames == 1:
            return getE_Platform(platformNameList[0])
        if numPlatformNames > 1 or numPlatformIds > 1:
            # print(f'Multiple platform names or Ids in samples for dataset {dsetName}')
            # fall through and return None
            pass
        elif numPlatformNames == 0 or numPlatformIds == 0:
            # print(f'No platform names or Ids in samples for dataset {dsetName}')
            # fall through and return None
            pass
    return None


def dsetPlatformMap_prev(data, outdir, matchSpecies=None):
    """
    Previous version of dsetPlatformMap which parses the platform_id for each
    dataset. This previous version used platform_id whereas the new version
    uses platform_accession_code as the platform identifier
    """
    speciesCounter = Counter()
    if matchSpecies:
        matchSpecies = matchSpecies.upper()
    with open(os.path.join(outdir, 'dset_map.txt'), 'w') as fp:
        for key, val in data.items():
            organisms = None
            if val['organism_names'] is not None:
                organisms = [species.upper() for species in val['organism_names']]
            speciesCounter.update(organisms)
            if matchSpecies and matchSpecies not in organisms:
                continue
            platformId = None
            # First try the platform_id list
            platformIdList = val.get("platform_ids")
            if platformIdList is not None and len(platformIdList) == 1:
                platformId = platformIdList[0]
                # check if platformId has spaces in it
                if bool(re.search(r"\s", platformId)):
                    # not a valid platformId, it has spaces in the name
                    platformId = None
            if platformId is None:
                # Next try the platform ids included with the samples
                platformId = getPlatformFromSamples(key, val)
            if platformId is None:
                # Finally try the platform_names list
                platformNames = val.get("platform_names")
                if platformNames is not None and len(platformNames) > 0:
                    numPlatforms = len(set(platformNames)) # consolidate identical items
                    if numPlatforms == 1:
                        platformId = getE_Platform(platformNames[0])
                    else:
                        assert numPlatforms > 1
                        # print(f'Multiple platform names for dataset {key}')
            if platformId is None or platformId == 'Missing':
                print(f'No or non-unique platform_id for dataset {key}')
                continue
            line = f"{key}\t{key}.{platformId}.pcl\t{key}.{platformId}\t{platformId}"
            fp.write(line + os.linesep)
    print(f"Species counts in file: {dict(speciesCounter)}")
    totalCount = sum(speciesCounter.values())
    print(f"Total dataset count: {totalCount}")


# Using platform_accession_code
def dsetPlatformMap(data, platNameMap, outdir, matchSpecies=None):
    """
    Function to get the platform identifier associated with each
    dataset, (in this case the platform_accession_code is used).
    The output is written to dset_map.txt file with on dataset per
    line and it's associated platform.
    """
    assert platNameMap is not None
    speciesCounter = Counter()
    countMatch = 0
    countMissing = 0
    if matchSpecies:
        matchSpecies = matchSpecies.upper()
    with open(os.path.join(outdir, 'dset_map.txt'), 'w') as fp:
        for key, val in data.items():
            organisms = None
            if val['organism_names'] is not None:
                organisms = [species.upper() for species in val['organism_names']]
            speciesCounter.update(organisms)
            if matchSpecies and matchSpecies not in organisms:
                continue
            platformId = None
            # First try the platform_accession_codes list
            platformAccList = val.get("platform_accession_codes")
            if platformAccList is not None and len(set(platformAccList)) == 1:
                platformId = platformAccList[0]
                # check if platformId has spaces in it
                if bool(re.search(r"\s", platformId)):
                    # not a valid platformId, it has spaces in the name
                    print(f"Space found in platformId: {platformId}")
                    platformId = None
            if platformId is None:
                # Try the platform_names list
                platformNames = val.get("platform_names")
                if platformNames is not None and len(platformNames) > 0:
                    numPlatforms = len(set(platformNames)) # consolidate identical items
                    if numPlatforms == 1:
                        platformId = platNameMap.get(platformNames[0])
                    else:
                        assert numPlatforms > 1
                        # print(f'Multiple platform names for dataset {key}')
            # if platformId is None:
            #     # Next try the platform ids included with the samples
            #     platformId = getPlatformFromSamples(key, val)
            if platformId is None or platformId == 'Missing':
                countMissing += 1
                print(f'No or non-unique platform_id for dataset {key}')
                continue
            countMatch += 1
            line = f"{key}\t{key}.{platformId}.pcl\t{key}.{platformId}\t{platformId}"
            fp.write(line + os.linesep)
    # print(f"Species counts in file: {dict(speciesCounter)}")
    print(f"CountMatch: {countMatch}, countMissing: {countMissing}")
    totalCount = sum(speciesCounter.values())
    print(f"Total dataset count: {totalCount}")


def createPlatNameMap(refPlatData):
    platNameMap = dict()
    for item in refPlatData:
        name = item.get('platform_name', None)
        accCode = item.get('platform_accession_code', None)
        if None not in [name, accCode]:
            if " " in accCode:
                accCode = accCode.replace(" ", "_")
            prevAcc = platNameMap.get(name, None)
            if prevAcc is not None and accCode != prevAcc:
                print(f"Multiple versions of platform {prevAcc} {accCode}")
                if "_" not in prevAcc:
                    # keep the previous one
                    continue
            platNameMap[name] = accCode
    return platNameMap


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--file', '-f', default='refine_bio_nonhuman_meta.json',
                        help="refine_bio json file to parse")
    parser.add_argument('--refPlatFile', '-p', default=None,
                        help="json file mapping platform name and platform_accession_code, "
                             "dowload from: https://api.refine.bio/v1/platforms/")
    parser.add_argument('--outdir', '-o', default='./',
                        help="output directory to write results")
    parser.add_argument('--species', '-s', default=None,
                        help="species to limit results to")
    parser.add_argument('--make-dset-map', '-m', default=False, action='store_true',
                           help='make dataset platform map')
    parser.add_argument('--make-species-dset-list', '-l', default=False, action='store_true',
                           help='make lists of datasets per species')
    args = parser.parse_args()

    with open(args.file) as fp:
        data = json.load(fp)

    platNameMap = None
    if args.refPlatFile is not None:
        with open(args.refPlatFile) as fp:
            refPlatData = json.load(fp)
        platNameMap = createPlatNameMap(refPlatData)

    if args.make_species_dset_list:
        perSpeciesDsetList(data, args.outdir)

    if args.make_dset_map:
        # dsetPlatformMap_prev(data, args.outdir, args.species.upper())
        assert platNameMap is not None
        dsetPlatformMap(data, platNameMap, args.outdir, args.species.upper())
