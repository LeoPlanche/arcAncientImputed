#!/usr/bin/env python3

import argparse
import os
from utils import createMap, computeSimilarities, computeDistances

    parser = argparse.ArgumentParser()
    subparser = parser.add_subparsers(dest = 'mode')

    # Create introgression map from decode folder
    createmap_subparser = subparser.add_parser('create_map', help='Create introgression map from decode folder')
    createmap_subparser.add_argument("-decode",help="path to the decode folder",type=str,required = True )
    createmap_subparser.add_argument("-out",help="outputfile (defaults to stdout)",type=str,required = True )
 
    # Add similarity to introgression map
    similarity_subparser = subparser.add_parser('similarity', help='Add similarity to introgression map')
    similarity_subparser.add_argument("-decode",
    help="path to the decode folder if -folder option is called, else path to a single decode file",
    type=str,required = True )

    similarity_subparser.add_argument('-p1', '--mapped_individual', required=True,
    type=str, action='store',
    help='VCF for the individual in the introgression map.')

    similarity_subparser.add_argument(
    '-p2', '--reference_individual', nargs='+',required=True,
    type=str, action='store',
    help='VCF for the individual to be compared to.',
    )

    similarity_subparser.add_argument(
    '-n2', '--name_reference', required=True,
    type=str, action='store',
    help='Name of the reference individual.',
    )   

    similarity_subparser.add_argument("-folder", 
    help="Apply the similarity to all files in the given folder, else to a single given file", 
    action='store_true')
 


    # Add distance to introgression map
    distance_subparser = subparser.add_parser('distance', help='Add distance to introgression map')
    distance_subparser.add_argument("-decode",
    help="path to the decode folder if -folder option is called, else path to a single decode file",
    type=str,required = True )

    distance_subparser.add_argument('-p1', '--mapped_individual', required=True,
    type=str, action='store',
    help='VCF for the individual in the introgression map.')

    distance_subparser.add_argument(
    '-p2', '--reference_individual', nargs='+',required=True,
    type=str, action='store',
    help='VCF for the individual to be compared to.',
    )

    distance_subparser.add_argument(
    '-n2', '--name_reference', required=True,
    type=str, action='store',
    help='Name of the reference individual.',
    )   

    distance_subparser.add_argument(
    '-p3', '--outgroup_individuals', nargs='+',required=True,
    type=str, action='store',
    help='VCF for the individual to be compared as an outgroup.',
    )

    distance_subparser.add_argument(
    '-n3', '--names_outgroup', required=True,
    type=str, action='store',
    help='Names of the reference individual, as a json file.',
    )

    distance_subparser.add_argument("-folder", 
    help="Apply the distance to all files in the given folder, else to a single given file", 
    action='store_true')

    args = parser.parse_args()

    if args.mode == 'create_map':
        createMap(args.decode, args.out)
    elif args.mode == 'similarity':
        if args.folder:
            for filename in os.listdir(args.decode):
                print(filename)
                file_path = os.path.join(args.decode, filename)
                computeSimilarities(file_path, args.mapped_individual, args.reference_individual, args.name_reference)
        else:
            computeSimilarities(args.decode, args.mapped_individual, args.reference_individual, args.name_reference)
    elif args.mode == 'distance':
        if args.folder:
            for filename in os.listdir(args.decode):
                print(filename)
                file_path = os.path.join(args.decode, filename)
                computeDistances(file_path, args.mapped_individual, args.reference_individual, args.name_reference, args.outgroup_individuals, args.names_outgroup)
        else:
            computeDistances(args.decode, args.mapped_individual, args.reference_individual, args.name_reference, args.outgroup_individuals, args.names_outgroup)
  
if __name__ == "__main__":
    main()