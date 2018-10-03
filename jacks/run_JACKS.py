import logging
from jacks.jacks_io import runJACKSFromArgs 
from jacks.infer import LOG

if __name__ == '__main__':
    LOG.setLevel(logging.INFO)
    runJACKSFromArgs()
