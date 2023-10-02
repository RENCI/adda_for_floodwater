#!/usr/bin/env python

# SPDX-FileCopyrightText: 2022 Renaissance Computing Institute. All rights reserved.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-License-Identifier: LicenseRef-RENCI
# SPDX-License-Identifier: MIT

#############################################################
#
# Retain from original utilities, only file IO, URL processing, and logging methods
# Add timing calls
# That can be used by any of the ADCIRC support tools
#
# Rebuilt logging method baseed in part on the work of P. Owen's Supervisor code
#############################################################

import sys,os
import yaml
import logging
import json

LOGGER = None
LOGFILE = None

class utilities:
    """
    Class to manage the logging setup.

    """
    @staticmethod
    def init_logging(subdir=None, config_file=None, log_file_metadata=None):
        """
        Initialize the Utilities class, set up logging
        """
        global LOGGER
        global LOGFILE
        config_data = utilities.load_config(yaml_file=config_file)
        if LOGGER is None:
            log,logfile = utilities.initialize_logging(subdir=subdir, config=config_data, log_file_metadata = log_file_metadata)
            LOGGER = log
            LOGFILE = logfile
        utilities.log = LOGGER
        utilities.LogFile = LOGFILE
        return config_data

    def initialize_logging(subdir=None, config=None, log_file_metadata = None):
        """
        Log file get saved to $LOG_PATH/subdir. LOG_PATH defaults to '.'

        Parameters
            config: dictionary containing logging settings (usually this is main.yml)
            subdir: A subdirectory constructed beneath the value in LOG_PATH 
        Returns
            logger handle
        """
        # logger = logging.getLogger(__name__)
        logger = logging.getLogger("ast_services") # We could simply add the instanceid here as well

        log_level = config["DEFAULT"].get('LOGLEVEL', 'DEBUG')
        logger.setLevel(log_level)

        if subdir is not None:
            Logdir = '/'.join([os.getenv('LOG_PATH','.'),subdir])
        else:
            Logdir = os.getenv('LOG_PATH','.')

        #LogName =os.getenv('LOG_NAME','logs')
        if log_file_metadata is None:
            LogName='AdcircSupportTools.log'
        else:    
            LogName=f'AdcircSupportTools_{log_file_metadata}.log'
        LogFile='/'.join([Logdir,LogName])

        formatter = logging.Formatter('%(asctime)s : %(levelname)s : %(funcName)s : %(module)s : %(name)s : %(message)s ')
        dirname = os.path.dirname(LogFile)
        if dirname and not os.path.exists(dirname):
            os.makedirs(dirname)
        file_handler = logging.FileHandler(LogFile, mode='w')
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
        return logger, LogFile
# YAML
    def load_config(yaml_file=None):
        if yaml_file is None:
            utilities.log.error('Called load_config but didnt specify yaml name: ABORT')
            sys.exit(1)
        if not os.path.exists(yaml_file):
            raise IOError("Failed to load yaml config file {}".format(yaml_file))
        with open(yaml_file, 'r') as stream:
            config = yaml.safe_load(stream)
            print('Opened yaml file {}'.format(yaml_file,))
        return config


