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


def dsetPlatformMap(data, outdir, matchSpecies=None):
    with open(os.path.join(outdir, 'dset_map.txt'), 'w') as fp:
        for key, val in data.items():
            if matchSpecies and matchSpecies not in val['organism_names']:
                continue
            platformIdList = val.get("platform_ids")
            if platformIdList is not None and len(platformIdList) > 0:
                if len(platformIdList) > 1:
                    print(f'More than one entry in platformIdList: {key}')
                    continue
                platformId = platformIdList[0]
            else:
                # fall back to getting from the platform name
                platformNames = val.get("platform_names")
                platformId = getE_Platform(platformNames[0])
                if platformId == 'Missing':
                    print(f'No platform_id for dataset {key}')
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
        dsetPlatformMap(data, args.outdir, args.species)


    