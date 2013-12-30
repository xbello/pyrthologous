from ConfigParser import ConfigParser
import os

from . import prepare

def config_parse(config_file):
    """Return the module with the values in config."""
   
    print config_file
    if not os.path.isfile(config_file):
        return False

    cfg = ConfigParser()
    cfg.read(config_file)

    d = {}
    for section in cfg.sections():
        d[section] = dict(cfg.items(section))

    return d

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("config", help="Name of the config file.")
    parser.add_argument("pairs", help="Name of the pairs file.")

    args = parser.parse_args()

    # Step 1: get the source of the genomes
    if args.config:
        cfg = config_parse(args.config)
