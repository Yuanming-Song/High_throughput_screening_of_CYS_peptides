#!/bin/bash
# Check for any unstaged or uncommitted changes
if [[ -n $(git status --porcelain) ]]; then
    echo "Unstaged or uncommitted changes detected."
    read -p "Do you want to commit your changes? (y/n): " answer
    if [[ "$answer" =~ ^[yY]$ ]]; then
        git add .
        read -p "Enter commit message: " msg
        git commit -m "$msg"
    else
        echo "Stashing changes..."
        git stash
    fi
fi

# Pull changes with fast-forward only
git pull --ff-only

# Push changes to remote
git push

# If changes were stashed (i.e., you didn't commit), pop the stash to restore them locally
if [[ "$answer" =~ ^[nN]$ ]]; then
    git stash pop
fi
