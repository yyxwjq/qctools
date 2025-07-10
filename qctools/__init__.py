# -*- coding: utf-8 -*-
import logging
import os

def qctools_logging(level=logging.INFO, filename='qctools.log'):
    if os.path.exists('qctools.log'):
        os.remove('qctools.log')
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        filename=filename,
        filemode='a'
    )
    
_default_logger = logging.getLogger(__name__)
_default_logger.addHandler(logging.NullHandler())  # Forbid "No handlers"
