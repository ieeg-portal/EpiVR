## Instructions as following steps
```python
# Check which branch you are on
git branch
# Create a new branch to start work or switch out of TRUNK
git checkout -b NAME-OF-YOUR-BRANCH
git pull origin NAME-OF-YOUR-BRANCH

# Do your work
...

# Ready to push your changes to the remote BRANCH
git pull origin NAME-OF-YOUR-BRANCH
git add -A
#OR
git add path/to/file1 path/to/file2
git commit -m "meaningful commit message"
git push origin NAME-OF-YOUR-BRANCH
```