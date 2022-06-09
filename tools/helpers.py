import os
import logging
from configparser import ConfigParser
import subprocess


def setup_logger(name, log_path=None):
    """Enables the logging to be setup correctly"""
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)

    stream_handle = logging.StreamHandler()
    stream_handle.setLevel(logging.DEBUG)
    stream_handle.setFormatter(
        logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    )
    logger.addHandler(stream_handle)

    if log_path:
        file_handle = logging.FileHandler(log_path, "a")
        file_handle.setLevel(logging.DEBUG)
        file_handle.setFormatter(
            logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
        )
        logger.addHandler(file_handle)

    return logger


def get_config():
    converters = {"list": lambda x: [i.strip() for i in x.split(",")]}
    config = ConfigParser(converters=converters)

    root_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    config_path = os.path.join(root_path, "config.ini")

    if not os.path.exists(config_path):
        raise FileNotFoundError(f"No config found at {config_path}")

    config.read(config_path)
    return config


def dir_to_samplesheet(scriptpath):
    subprocess.run(scriptpath)
