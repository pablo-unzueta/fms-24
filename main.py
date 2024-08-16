# this should be a runnable script

import argparse

parser = argparse.ArgumentParser(description='FMS24')
parser.add_argument('--config', type=str, default='config.yaml', help='config file')
args = parser.parse_args()