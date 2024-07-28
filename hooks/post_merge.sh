mkdir -p .git/hooks
cat << 'EOF' > .git/hooks/post-merge
#!/bin/bash
if [ ! -f config.sh ]; then
    echo "It looks like config.sh is missing. Please copy config.template.sh to config.sh and customize it."
fi
EOF
chmod +x .git/hooks/post-merge
