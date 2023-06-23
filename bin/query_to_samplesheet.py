#!/usr/bin/env python

from pymongo import MongoClient
import pandas
import csv
import os
import argparse
import sys
import logging

def get_sample(filter,max_number):
    # Database connection
    '''
    TODO: Use pipeline or environment variable(s) to fill authentication data.
    '''
    hostname = "10.9.63.24"
    port = 27017
    username = "user01"
    password = "mtb_project2023"
    database_name = "mtb_database"
    collection = "metadata_db"

    mongo_uri = f"mongodb://{username}:{password}@{hostname}:{port}"

    mongo_client = MongoClient(mongo_uri)

    mongo_db = mongo_client[database_name]

    mongo_coll = mongo_db[collection]

    # Query

    query_fields = {
        "Sample_Name" : 1,
        "files_name" : 1
    }
    if max_number == None:
        return mongo_coll.find(filter, query_fields)
    else:
        return mongo_coll.find(filter, query_fields).limit(max_number)

def make_samplesheet(mongo_obj):
    # Convert to nextflow pipeline input format as csv
    '''
    Pipeline input format:
    sample,fastq_1,fastq_2
    SRR000001,SRR000001_1.fastq.gz,SRR000001_2.fastq.gz
    SRR000002,SRR000002.fastq.gz,
    '''

    fieldnames = [ "sample", "fastq_1", "fastq_2" ]
    '''
    TODO: Decide on share path to fastq later. Ideally handed by database.
    '''
    fastq_path = "/home/thavin/nas.folder/test"
    with open("samplesheet.csv", "w", newline="") as csvfile:
        w = csv.DictWriter(csvfile, fieldnames=fieldnames, extrasaction="raise")
        w.writeheader()
        for record in mongo_obj:
            sample = record['Sample_Name']
            fastq = record['files_name'].strip('][').replace("\'","").split(', ')

            if len(fastq) > 2:
                logger.error(f"[{sample}]: More than 2 fastq file was found!!")
                sys.exit(2)

            fastq_1 = os.path.join(fastq_path,fastq[0])
            if len(fastq) == 1:
                fastq_2 = ""
            elif len(fastq) > 1:
                fastq_2 = os.path.join(fastq_path,fastq[1])

            w.writerow({"sample":sample, "fastq_1":fastq_1, "fastq_2":fastq_2})


def parse_args(argv=None):
    """
    Define and parse command line arguments.
    TODO: Parse positional/kw argument for filter criteria in query term to MTB databases (prolly from string)
          (e.g. Lineage, Cluster assignment, Location) into dictionary type.
    TODO: Add criteria to query only unanalysed MTB records. (Currently this query all from database).
    """
    parser = argparse.ArgumentParser(
        description="Query new fastq of MTB from mongodb to create the samplesheet for the clustering pipeline",
        epilog="Example: python check_samplesheet.py",
    )
    parser.add_argument(
        "-f",
        "--filter",
        metavar="filter",
        default={},
        help="Filter use in the query of MTB",
    )

    parser.add_argument(
        "-m",
        "--max-number",
        metavar="max_number",
        type=int,
        default=None,
        help="limiting the number of query"

    )

    parser.add_argument(
        "-l",
        "--log-level",
        help="The desired log level (default WARNING).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="WARNING",
    )
    return parser.parse_args(argv)


def main(argv=None):
    """
    Main function
    parse argumnet to get_sample()
    """
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")
    query_data = get_sample(args.filter,args.max_number)
    make_samplesheet(query_data)


if __name__ == "__main__":
    sys.exit(main())


