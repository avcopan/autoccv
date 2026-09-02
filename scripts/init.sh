#!/usr/bin/env bash

# Set up direnv, if present
if command -v direnv >/dev/null 2>&1; then
  direnv allow
else
  echo "Warning: direnv is not installed."
  echo "You can install direnv and then run: direnv allow"
fi

# Initialize git repository
git init
git add --all
if git config --get user.name >/dev/null 2>&1 && git config --get user.email >/dev/null 2>&1; then
  git commit -m "Initial commit"
else
  echo "Warning: Git author identity is not configured."
  echo "You can configure it as follows:"
  echo '  git config --global user.name "Your Name"'
  echo '  git config --global user.email "you@example.com"'
  echo "Afterwards, you can make the initial commit with:"
  echo '  git commit -m "Initial commit"'
fi

# Install git hooks with lefthook
lefthook install

# Install all pixi environments
pixi install --all

USER=$(git config --global user.name)
REPO=$(basename "$PWD")
echo "Setup complete."
echo
echo "Now, run the following commands to set up your remote repository:"
echo "  git remote add origin git@github.com:$USER/$REPO.git"
echo "  git push -u origin main"
echo "(Double-check the URL to make sure the user/org and repo names are what you want.)"
