#!/usr/bin/env python3
import os
import sys
from Bio import SeqIO   # to read bio datasets
import sqlalchemy
'''
    as same as the given snippet code to import the correct version
'''
if sqlalchemy.__version__.startswith('1.4'):
    from sqlalchemy import create_engine
    from sqlalchemy.engine import make_url, URL
    from sqlalchemy.sql import text
elif sqlalchemy.__version__.startswith('2'):
    from sqlalchemy import create_engine
    from sqlalchemy.url import make_url, URL
    from sqlalchemy.sql import text

# as for the assignment mentioned, I read it in hardcode way, otherwise I can read it from CLI
genbank_file = '/homes/aakhlaghi/archaea.1.genomic.gbff'

'''(optional:CLI)
import argparse
def parse_args():
    p = argparse.ArgumentParser(description="Parse GenBank files and store basic stats in MariaDB/MySQL.")
    p.add_argument("files", nargs="+", help="Path(s) to .gb or .gbff file(s)")  # nargs to add one or multiple files
    return p.parse_args()
    '''
    
    
