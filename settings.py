# settings.py
import os
from dotenv import load_dotenv

env_path = Path('.') / '.env'
load_dotenv(dotenv_path=env_path, verbose=True)

DATA_DIR = os.getenv("DATA_DIR")
COMP_DIR = os.getenv("COMP_DIR")
