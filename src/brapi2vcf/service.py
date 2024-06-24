import argparse
import logging
from . import brapi2vcf

parser = argparse.ArgumentParser(description="Construct a VCF-file from BrAPI-endpoint data.")
parser.add_argument("url", type=str, help="url from BrAPI server")
args = parser.parse_args()

logging.basicConfig(format="%(asctime)s | %(name)s |  %(levelname)s: %(message)s",
                    datefmt="%m-%d-%y %H:%M:%S")
logging.getLogger("brapi2vcf").setLevel(logging.DEBUG)

def service():
    entry = brapi2vcf.BrAPI2vcf(args.url)
    return entry.vcf()