import os
import sys
import logging

def error(msg):
    logging.error(msg)
    sys.exit(1)

def log(msg):
    logging.log(msg)
