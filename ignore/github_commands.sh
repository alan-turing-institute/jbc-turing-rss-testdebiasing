# Can configure username
git config --global user.name "georgenicholson"
git config --global user.password "My password here"

# # Clone files from origin (not sure how often need to do this)
# git clone https://github.com/alan-turing-institute/jbc-turing-rss-data.git
# git clone https://github.com/alan-turing-institute/jbc-turing-rss-testdebiasing.git

# cd to Git controlled directory
# cd /mnt/c/Users/nicho/Documents/bauer_sync/projects/covid/covid_github/jbc-turing-rss-data
# mkdir /mnt/c/Users/nicho/Documents/bauer_sync/projects/covid/covid_github/jbc-turing-rss-testdebiasing-1
# cd /mnt/c/Users/nicho/Documents/bauer_sync/projects/covid/covid_github/jbc-turing-rss-testdebiasing-1
cd /mnt/c/Users/nicho/Documents/bauer_sync/projects/covid/covid_github/jbc-turing-rss-testdebiasing

# cd /mnt/c/Users/nicho/Documents/bauer_sync/projects/covid/covid_github


# Look at local branches
git branch
# switch to local main branch
git checkout develop
# check everything up to date
git remote update && git status 
# Make the local master sync with remote master
git pull origin develop

git merge develop

# Checkout a local branch to work on
git checkout -b george
# Switch to local george branch
git checkout george

# git checkout main

#################################
# Make changes here
#################################
git pull origin george
# Switch to local george branch
git checkout george
git add -A # Add all files to be committed
git commit -m 'Updating data sources, now with vaccination updated to 13 May' # Commit files with message
git push --set-upstream origin george

# switch to local main branch
git checkout develop
# check everything up to date
git status
# delete local george branch
git branch -d george  
# Look at local branches
git branch

git fetch origin
git checkout -b george origin/george
git merge main



# Clone REACT repo
# # Clone files from origin (not sure how often need to do this)
# cd to Git controlled directory
cd /mnt/c/Users/nicho/Documents/bauer_sync/projects/covid/covid_github
git clone https://github.com/georgenicholson/reactidd.git


