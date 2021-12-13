import os
import json
import argparse
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
        platformList = []
        for key, val in samples:
            plat = val.get('platform_name')
            if plat is not None:
                platformList.push(plat)
        numPlatforms = len(set(platformList)) # consolidate identical items
        if numPlatforms == 1:
            return platformList[0]
        elif numPlatforms > 1:
            print(f'Multiple platform names in samples for dataset {dsetName}')
            # fall through and return None
        elif numPlatforms == 0:
            print(f'No platform names in samples for dataset {dsetName}')
            # fall through and return None
    return None


def dsetPlatformMap(data, outdir, matchSpecies=None):
    with open(os.path.join(outdir, 'dset_map.txt'), 'w') as fp:
        for key, val in data.items():
            if matchSpecies and matchSpecies not in val['organism_names']:
                continue
            platformIdList = val.get("platform_ids")
            platformId = None
            if platformIdList is not None and len(platformIdList) > 0:
                if len(platformIdList) > 1:
                    print(f'More than one entry in platformIdList: {key}')
                    continue
                platformId = platformIdList[0]
            else:
                # fall back to getting from the platform name
                platformNames = val.get("platform_names")
                if platformNames is not None and len(platformNames) > 0:
                    numPlatforms = len(set(platformNames)) # consolidate identical items
                    if numPlatforms == 1:
                        platformId = getE_Platform(platformNames[0])
                    else:
                        assert numPlatforms > 1
                        print(f'Multiple platform names for dataset {key}')
                else:
                    # fall back to looking at the samples
                    platformId = getPlatformFromSamples(key, val)
                if platformId is None or platformId == 'Missing':
                    print(f'No platform_id for dataset {key}')
                    continue
            line = f"{key}\t{key}.{platformId}.pcl\t{key}.{platformId}\t{platformId}"
            fp.write(line + os.linesep)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--file', '-f', default='refine_bio_nonhuman_meta.json',
                        help="refine_bio json file to parse")
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

    if args.make_species_dset_list:
        perSpeciesDsetList(data, args.outdir)

    if args.make_dset_map:
        dsetPlatformMap(data, args.outdir, args.species.upper())


    