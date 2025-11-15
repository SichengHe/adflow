# Moving ADflow with Resolvent Analysis to Private Repository

## Overview

This guide explains how to move your current ADflow repository (with the new resolvent analysis implementation) to a private repository while preserving all history and new features.

## Quick Steps

### Option 1: Create New Private Repo and Push

```bash
# Navigate to your ADflow directory
cd /home/sicheng/repo/adflow

# Check current remote
git remote -v

# Create a new private repository on GitHub/GitLab (via web interface)
# Then add it as a new remote
git remote add private git@github.com:YOUR_USERNAME/adflow-private.git

# Push everything to private repo
git push private --all
git push private --tags

# Optional: Set private as default remote
git remote set-url origin git@github.com:YOUR_USERNAME/adflow-private.git
```

### Option 2: Fork and Make Private (GitHub)

```bash
# If the original is on GitHub:
# 1. Go to original repo on GitHub
# 2. Click "Fork"
# 3. Go to Settings > Change visibility > Make private

# Then update your local remote
cd /home/sicheng/repo/adflow
git remote set-url origin git@github.com:YOUR_USERNAME/adflow-private.git
```

### Option 3: Mirror to New Repository

```bash
cd /home/sicheng/repo/adflow

# Create bare clone
cd ..
git clone --bare adflow adflow-bare

# Create new private repo on GitHub/GitLab
# Then mirror push
cd adflow-bare
git push --mirror git@github.com:YOUR_USERNAME/adflow-private.git

# Clone the new private repo
cd ..
git clone git@github.com:YOUR_USERNAME/adflow-private.git adflow-private
cd adflow-private

# Verify everything is there
git log --oneline | head -20
ls -la src/modalAnalysis/
```

## Detailed Step-by-Step Guide

### Step 1: Save Current State

```bash
cd /home/sicheng/repo/adflow

# Make sure everything is committed
git status

# Commit any uncommitted changes
git add -A
git commit -m "Add resolvent analysis implementation

- New modalAnalysis module with resolventAPI.F90
- Python integration in pyADflow.py and pyResolventAnalysis.py
- Complete documentation and examples
- Based on resolvent analysis paper by He et al."

# Create a tag for this version
git tag -a v2.13.0-resolvent -m "ADflow with resolvent analysis"
```

### Step 2: Create Private Repository

**On GitHub:**
1. Go to https://github.com/new
2. Repository name: `adflow-private` (or your choice)
3. Description: "ADflow CFD solver with resolvent analysis (private)"
4. Select **Private**
5. Don't initialize with README/license (we're pushing existing repo)
6. Click "Create repository"

**On GitLab:**
1. Go to https://gitlab.com/projects/new
2. Project name: `adflow-private`
3. Visibility Level: **Private**
4. Don't initialize with README
5. Click "Create project"

### Step 3: Push to Private Repository

```bash
cd /home/sicheng/repo/adflow

# Add the new private remote
git remote add private git@github.com:YOUR_USERNAME/adflow-private.git
# Or for GitLab:
# git remote add private git@gitlab.com:YOUR_USERNAME/adflow-private.git

# Push all branches
git push private --all

# Push all tags
git push private --tags

# Verify
git remote -v
```

### Step 4: Update Local Configuration

```bash
# Option A: Keep both remotes (recommended)
git remote -v
# origin  git@github.com:mdolab/adflow.git (fetch)
# origin  git@github.com:mdolab/adflow.git (push)
# private git@github.com:YOUR_USERNAME/adflow-private.git (fetch)
# private git@github.com:YOUR_USERNAME/adflow-private.git (push)

# Push to private by default, pull from origin
git config branch.master.remote private

# Option B: Replace origin with private
git remote remove origin
git remote rename private origin
```

### Step 5: Verify Everything Transferred

```bash
# Clone the private repo in a new location to verify
cd /tmp
git clone git@github.com:YOUR_USERNAME/adflow-private.git test-clone
cd test-clone

# Check recent commits
git log --oneline -10

# Check new files are there
ls -la src/modalAnalysis/
ls -la adflow/pyResolventAnalysis.py
cat RESOLVENT_IMPLEMENTATION_COMPLETE.md

# Check tags
git tag

# Clean up test
cd ..
rm -rf test-clone
```

## What Gets Transferred

### All New Resolvent Files ✓
```
src/modalAnalysis/
├── README.md
└── resolventAPI.F90

adflow/
├── __init__.py (modified)
├── pyADflow.py (modified)
└── pyResolventAnalysis.py (new)

doc/
└── resolvent_analysis.md

examples/
└── resolvent_analysis_example.py

Project root:
├── RESOLVENT_IMPLEMENTATION.md
├── RESOLVENT_INTEGRATION_PLAN.md
├── RESOLVENT_JACOBIAN_NOTES.md
├── RESOLVENT_IMPLEMENTATION_COMPLETE.md
└── RESOLVENT_QUICKSTART.md
```

### All Git History ✓
- Complete commit history
- All branches
- All tags
- Author information
- Timestamps

## Maintaining Relationship with Public ADflow

### Strategy 1: Dual Remotes (Recommended)

Keep connection to both public and private repos:

```bash
# Your setup
git remote -v
# origin  git@github.com:mdolab/adflow.git (public - for updates)
# private git@github.com:YOUR_USERNAME/adflow-private.git (your work)

# Pull updates from public ADflow
git fetch origin
git merge origin/master

# Push your work to private
git push private master

# When ready to contribute back
git push origin feature-branch
# Then create PR on GitHub
```

### Strategy 2: Fork Workflow

```bash
# Public repo is "upstream"
git remote add upstream git@github.com:mdolab/adflow.git
git remote set-url origin git@github.com:YOUR_USERNAME/adflow-private.git

# Get updates
git fetch upstream
git merge upstream/master

# Your work stays private
git push origin master
```

## Handling Sensitive Information

### Before Making Public (If Ever)

If you plan to eventually contribute resolvent analysis back:

```bash
# Check for any sensitive data
git log --all --full-history --source -- config/*.secret

# Remove sensitive files from history (if any)
git filter-branch --tree-filter 'rm -f path/to/sensitive/file' HEAD

# Or use git-filter-repo (recommended)
pip install git-filter-repo
git filter-repo --path path/to/sensitive/file --invert-paths
```

### Files to Consider

Check these for sensitive info before going public:
- Configuration files with credentials
- API keys or tokens
- Internal server addresses
- Proprietary data or results

## Backup Strategy

### Create Local Backup

```bash
# Full repository backup
cd /home/sicheng/repo
tar -czf adflow-resolvent-backup-$(date +%Y%m%d).tar.gz adflow/

# Move to safe location
mv adflow-resolvent-backup-*.tar.gz ~/backups/
```

### Create Archive on GitHub

```bash
# Create a release with source code archive
git tag -a v2.13.0-resolvent -m "Resolvent analysis implementation"
git push private --tags

# On GitHub: Go to Releases > Create new release
# This creates downloadable archives
```

## Collaboration Setup

### Add Collaborators to Private Repo

**On GitHub:**
1. Go to repo Settings > Collaborators
2. Add team members by username
3. Choose permission level (Write, Admin, etc.)

### Clone for Team Members

```bash
# Team members clone private repo
git clone git@github.com:YOUR_USERNAME/adflow-private.git

# Set up environment
cd adflow-private
# ... follow ADflow setup instructions
```

## CI/CD Considerations

### Update CI Configuration

If you have CI/CD (Travis, GitHub Actions, etc.):

```yaml
# .github/workflows/tests.yml
name: Tests

on:
  push:
    branches: [ master, develop ]
  pull_request:
    branches: [ master ]

jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v2
      - name: Run tests
        run: |
          # Your test commands
          python -m pytest tests/
```

## Documentation Updates

### Update README

Add a section to README.md:

```markdown
## Resolvent Analysis Feature

This private fork includes resolvent analysis capabilities for
frequency-domain stability analysis. See:
- `doc/resolvent_analysis.md` for theory and usage
- `RESOLVENT_QUICKSTART.md` for quick start guide
- `examples/resolvent_analysis_example.py` for example

Based on: "Large-Scale Flow Control Performance Optimization via
Differentiable Resolvent Analysis" by He et al.
```

## Troubleshooting

### Large File Issues

```bash
# If repo is too large
git count-objects -vH

# Find large files
git rev-list --objects --all | \
  git cat-file --batch-check='%(objecttype) %(objectname) %(objectsize) %(rest)' | \
  awk '/^blob/ {print substr($0,6)}' | \
  sort --numeric-sort --key=2 | \
  tail -20
```

### Permission Errors

```bash
# If git push fails with permission error
ssh -T git@github.com  # Test SSH connection

# Generate new SSH key if needed
ssh-keygen -t ed25519 -C "your_email@example.com"

# Add to GitHub: Settings > SSH and GPG keys
```

## Summary Checklist

- [ ] Commit all resolvent analysis changes
- [ ] Create annotated tag for this version
- [ ] Create private repository on GitHub/GitLab
- [ ] Add private remote to local repo
- [ ] Push all branches and tags
- [ ] Verify files transferred correctly
- [ ] Update local git configuration
- [ ] Add collaborators (if any)
- [ ] Update README/documentation
- [ ] Create backup archive
- [ ] Test building and running

## Commands Summary

```bash
# Quick migration
cd /home/sicheng/repo/adflow
git add -A
git commit -m "Add resolvent analysis implementation"
git tag -a v2.13.0-resolvent -m "Resolvent analysis"
git remote add private git@github.com:YOUR_USERNAME/adflow-private.git
git push private --all
git push private --tags

# Verify
git clone git@github.com:YOUR_USERNAME/adflow-private.git /tmp/test
ls -la /tmp/test/src/modalAnalysis/
```

## Next Steps After Migration

1. **Build System Integration**: Add `resolventAPI.F90` to Makefile
2. **Testing**: Run tests on private repo
3. **Team Setup**: Add collaborators and set permissions
4. **CI/CD**: Configure automated testing
5. **Development**: Continue development in private repo

---

**Note**: This guide assumes you have write access to create private repositories
on your chosen platform (GitHub, GitLab, etc.).
