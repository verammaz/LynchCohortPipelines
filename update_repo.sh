#!/bin/bash

# File to be protected
PROTECTED_FILE="config.sh"

# Function to reset changes to all files except the protected one
reset_except_protected() {
    # Stash untracked files as well
    git stash -u

    # Pull the latest changes from the remote repository
    git pull

    # Apply the stashed changes
    git stash pop

    # Revert any changes to the protected file to the latest commit state
    git checkout HEAD -- "$PROTECTED_FILE"
}

# Check if the protected file has local changes
if git diff --name-only | grep -q "^$PROTECTED_FILE$"; then
    echo "Local changes to $PROTECTED_FILE found. Stashing..."
    git stash push -m "Protected changes" -- "$PROTECTED_FILE"
fi

# Reset changes to all files except the protected one
reset_except_protected

# Make all scripts in the repository executable
find . -type f -name "*.sh" -exec chmod +x {} \;

echo "Update complete."
