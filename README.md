# CODE REFACTOR IN PROGRESS
Create virtualenv virtual environment and install the necessary requisites using the following:
```virtualenv .env && source .env/bin/activate && pip install -r requirements.txt```

Test using `pytest` and `coverage.py` using:
```pytest app/ --cov=app/ --cov-report=html```
