# -*- coding: utf-8 -*-
# Copyright (c) 2019 Center for Neuroengineering and Therapeutics.
# Copyright (c) 2019 University of Pennsylvania.
# License: MIT, see https://mit-license.org/

import os
from dotenv import load_dotenv

env_path = os.path.join(os.getcwd(), '.env')
load_dotenv(dotenv_path=env_path, verbose=True)

COMP_DIR = os.getenv("COMP_DIR")
DATA_DIR = os.getenv("DATA_DIR")