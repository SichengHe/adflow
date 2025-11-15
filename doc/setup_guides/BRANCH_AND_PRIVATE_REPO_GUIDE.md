# Complete Guide: Branch + Private Repository

## Quick Start (Recommended Approach)

```bash
cd /home/sicheng/repo/adflow

# Option 1: Automated setup
./setup_resolvent_branch.sh

# Option 2: Manual setup
git checkout -b feature/resolvent-analysis
git add src/modalAnalysis/ adflow/ doc/ examples/ RESOLVENT_*.md
git commit -m "Add resolvent analysis implementation"
git tag -a v2.13.0-resolvent -m "Resolvent analysis"

# Then push to private repo
git remote add private git@github.com:YOUR_USERNAME/adflow-private.git
git push -u private feature/resolvent-analysis
git push private master
git push private --tags
```

## The Complete Workflow

### Step 1: Create Feature Branch ✓

```bash
cd /home/sicheng/repo/adflow

# Create branch for resolvent work
git checkout -b feature/resolvent-analysis

# Add all resolvent files
git add src/modalAnalysis/
git add adflow/pyResolventAnalysis.py
git add adflow/__init__.py adflow/pyADflow.py
git add doc/resolvent_analysis.md
git add examples/resolvent_analysis_example.py
git add RESOLVENT_*.md

# Commit
git commit -m "Add resolvent analysis implementation

[Full message with details...]"

# Tag
git tag -a v2.13.0-resolvent -m "Resolvent analysis implementation"
```

### Step 2: Create Private Repository ✓

On GitHub/GitLab:
1. Create new repository: `adflow-private`
2. Set visibility: **Private**
3. Don't initialize with README

### Step 3: Push Branch to Private Repo ✓

```bash
# Add private remote
git remote add private git@github.com:YOUR_USERNAME/adflow-private.git

# Push feature branch
git push -u private feature/resolvent-analysis

# Also push master (to keep in sync with public ADflow)
git push private master

# Push tags
git push private --tags
```

### Step 4: Set Up Remotes ✓

```bash
# Verify remotes
git remote -v
# origin  git@github.com:mdolab/adflow.git (public ADflow)
# private git@github.com:YOUR_USERNAME/adflow-private.git (your work)

# Configure branch tracking
git branch -vv
# * feature/resolvent-analysis  [private/feature/resolvent-analysis]
#   master                       [origin/master]
```

## Repository Structure

### Your Private Repo

```
private: git@github.com:YOUR_USERNAME/adflow-private.git
├── master                      # Synced with public ADflow
└── feature/resolvent-analysis  # Your resolvent work
    └── All resolvent files ✓
```

### Public ADflow Repo

```
origin: git@github.com:mdolab/adflow.git
└── master                      # Official ADflow (no resolvent)
```

## Day-to-Day Workflow

### Working on Resolvent

```bash
# Switch to feature branch
git checkout feature/resolvent-analysis

# Make changes
vim src/modalAnalysis/resolventAPI.F90

# Commit
git add src/modalAnalysis/resolventAPI.F90
git commit -m "Improve Jacobian assembly performance"

# Push to private repo
git push private feature/resolvent-analysis
```

### Syncing with Public ADflow

```bash
# Update master from public repo
git checkout master
git fetch origin
git merge origin/master

# Push updated master to private
git push private master

# Update feature branch with master changes
git checkout feature/resolvent-analysis
git merge master

# Push updated feature branch
git push private feature/resolvent-analysis
```

### Switching Between Branches

```bash
# Work on resolvent
git checkout feature/resolvent-analysis
ls src/modalAnalysis/  # ✓ Resolvent files present
make  # Build with resolvent

# Switch to clean ADflow
git checkout master
ls src/modalAnalysis/  # Directory doesn't exist
make  # Build standard ADflow

# Back to resolvent
git checkout feature/resolvent-analysis
ls src/modalAnalysis/  # ✓ Files back
```

## Automated Setup Script

Use the provided script:

```bash
./setup_resolvent_branch.sh
```

This script:
- ✓ Creates `feature/resolvent-analysis` branch
- ✓ Stages all resolvent files
- ✓ Creates comprehensive commit
- ✓ Creates version tag
- ✓ Provides next-step instructions

## Benefits of This Approach

### 1. Clean Separation ✓

```
master branch:
- Standard ADflow (synced with public)
- No resolvent files
- Can always build standard version

feature/resolvent-analysis:
- All resolvent implementation
- Your development work
- Isolated from master
```

### 2. Easy Updates ✓

```bash
# Get public ADflow updates
git checkout master
git pull origin master
git push private master

# Incorporate updates into your work
git checkout feature/resolvent-analysis
git merge master
```

### 3. Flexible Merging ✓

```bash
# Merge to private master when ready
git checkout master
git merge feature/resolvent-analysis

# Or keep separate indefinitely
# (develop on branch, master stays clean)
```

### 4. Future Contributions ✓

```bash
# Easy to create PR to public ADflow
git push origin feature/resolvent-analysis
# Then create PR on GitHub: mdolab/adflow ← feature/resolvent-analysis
```

## Advanced Workflows

### Multiple Sub-Branches

```bash
feature/resolvent-analysis      # Main branch
├── feature/resolvent-adjoint   # Adjoint derivatives (future)
├── feature/resolvent-matfree   # Matrix-free (future)
└── feature/resolvent-slepc     # SLEPc integration (future)

# Create sub-branch
git checkout feature/resolvent-analysis
git checkout -b feature/resolvent-adjoint

# Merge back when done
git checkout feature/resolvent-analysis
git merge feature/resolvent-adjoint
```

### Release Branches

```bash
# Create release branch
git checkout feature/resolvent-analysis
git checkout -b release/2.13.0-resolvent

# Final testing, documentation updates
# ...

# Tag release
git tag -a v2.13.0-resolvent-release -m "Production release"

# Merge to master
git checkout master
git merge release/2.13.0-resolvent
```

## Collaboration

### Add Team Members

On GitHub private repo:
1. Settings > Collaborators
2. Add by username
3. Choose permission level

### Clone for Team

```bash
# Team member clones
git clone git@github.com:YOUR_USERNAME/adflow-private.git
cd adflow-private

# Work on feature branch
git checkout feature/resolvent-analysis

# Make changes and push
git add <files>
git commit -m "Description"
git push origin feature/resolvent-analysis
```

## Comparison Table

| Aspect | Master Branch | Feature Branch | Private Repo |
|--------|--------------|----------------|--------------|
| Content | Standard ADflow | +Resolvent | All branches |
| Synced with | Public ADflow | Master + resolvent | Your work |
| Development | No changes | Active work | Container |
| Build | Standard ADflow | With resolvent | Both versions |
| Purpose | Stable baseline | Development | Storage/sharing |

## Commands Cheat Sheet

```bash
# Setup
git checkout -b feature/resolvent-analysis
git add <resolvent-files>
git commit -m "Add resolvent"
git remote add private git@github.com:USER/adflow-private.git
git push -u private feature/resolvent-analysis

# Daily work
git checkout feature/resolvent-analysis
# ... edit files ...
git add <files>
git commit -m "Description"
git push private feature/resolvent-analysis

# Sync with public
git checkout master
git pull origin master
git push private master
git checkout feature/resolvent-analysis
git merge master

# View differences
git diff master..feature/resolvent-analysis --stat

# Switch branches
git checkout master                    # Clean ADflow
git checkout feature/resolvent-analysis # With resolvent
```

## Troubleshooting

### Branch Already Exists

```bash
# If branch exists locally
git checkout feature/resolvent-analysis

# If exists remotely but not locally
git fetch private
git checkout -b feature/resolvent-analysis private/feature/resolvent-analysis
```

### Merge Conflicts

```bash
# When merging master into feature branch
git checkout feature/resolvent-analysis
git merge master
# CONFLICT in <file>

# Resolve conflicts
vim <conflicted-file>
git add <resolved-file>
git commit
```

### Wrong Branch

```bash
# If you committed to master by mistake
git checkout master
git log --oneline -5  # Find the commit

# Move commit to feature branch
git checkout feature/resolvent-analysis
git cherry-pick <commit-hash>

# Remove from master
git checkout master
git reset --hard HEAD~1
```

## Summary

**You now have:**
- ✅ Feature branch: `feature/resolvent-analysis`
- ✅ Private repository: All your work
- ✅ Clean master: Synced with public ADflow
- ✅ Tagged version: `v2.13.0-resolvent`

**Benefits:**
- Clean separation of features
- Easy to sync with public ADflow
- Flexible development workflow
- Simple collaboration
- Ready for future contributions

**Next Steps:**
1. Push to private repo (if not done)
2. Continue development on feature branch
3. Test and refine implementation
4. Eventually merge to master (your private master)
5. Optionally contribute back to public ADflow

See individual guides for details:
- `RESOLVENT_BRANCH_SETUP.md` - Branch details
- `MIGRATION_TO_PRIVATE_REPO.md` - Private repo setup
- `RESOLVENT_QUICKSTART.md` - Usage guide
