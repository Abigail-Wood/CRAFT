import os
import logging

def error(msg):
    logging.error(msg)
    os.exit(1)

def log(msg):
    logging.log(msg)
