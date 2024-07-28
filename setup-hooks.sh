cat << 'EOF' > setup-hooks.sh
#!/bin/bash
cp hooks/post-merge .git/hooks/post-merge
chmod +x .git/hooks/post-merge
EOF
chmod +x setup-hooks.sh
