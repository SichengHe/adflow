# Resolvent Analysis Branch Setup

## Overview

This guide shows how to organize the resolvent analysis implementation on a dedicated branch. This is recommended for:
- Clean separation of features
- Easier code review
- Organized development workflow
- Simple merging later

## Quick Setup

```bash
cd /home/sicheng/repo/adflow

# Create and switch to resolvent branch
git checkout -b feature/resolvent-analysis

# Check what files have changed
git status

# Add all resolvent-related files
git add src/modalAnalysis/
git add adflow/pyResolventAnalysis.py
git add adflow/__init__.py
git add adflow/pyADflow.py
git add doc/resolvent_analysis.md
git add examples/resolvent_analysis_example.py
git add RESOLVENT_*.md

# Commit on the branch
git commit -m "Add resolvent analysis implementation

Complete implementation of resolvent analysis for ADflow including:

- Fortran API (src/modalAnalysis/resolventAPI.F90)
  * setupResolventMatrix() - Assemble Jacobian
  * getResolventMatrixDense() - Export J = ∂R/∂w
  * getResolventMatrixInfo() - Matrix information
  * exportResolventMatrixToFile() - PETSc binary export

- Python integration (adflow/)
  * ResolventAnalysis class with solveExplicit()
  * ResolventAnalysisMatrixFree for large problems
  * Complex and real doubled formulations
  * setupResolventJacobian() in ADFLOW class
  * getJacobianMatrix() method

- Documentation
  * doc/resolvent_analysis.md - Theory and usage
  * RESOLVENT_*.md files - Implementation guides
  * examples/resolvent_analysis_example.py

Based on: 'Large-Scale Flow Control Performance Optimization via
Differentiable Resolvent Analysis' by He et al.

Addresses:
- Proper Jacobian transpose handling
- Modal analysis organization (not adjoint)
- Both complex and real formulations"

# Tag this version
git tag -a v2.13.0-resolvent -m "Resolvent analysis implementation"

# Push to remote (after setting up private repo)
git push -u private feature/resolvent-analysis
git push private --tags
```

## Recommended Branch Structure

```
master (or main)
├── develop                          # Main development branch
└── feature/resolvent-analysis       # Resolvent implementation ← You are here
    ├── docs/resolvent-theory       # (Optional) Documentation updates
    ├── feature/resolvent-adjoint   # (Future) Adjoint derivatives
    └── feature/resolvent-matfree   # (Future) Matrix-free methods
```

## Detailed Workflow

### Step 1: Create Feature Branch

```bash
cd /home/sicheng/repo/adflow

# Make sure master is clean
git checkout master
git status

# Create feature branch from master
git checkout -b feature/resolvent-analysis

# Verify you're on the new branch
git branch
# * feature/resolvent-analysis
#   master
```

### Step 2: Organize Your Changes

Check what files you've modified:

```bash
# See modified files
git status

# See detailed changes
git diff master --stat

# Review specific file changes
git diff master src/modalAnalysis/resolventAPI.F90
```

### Step 3: Stage Files Strategically

```bash
# Stage by category for clear commits

# 1. Fortran implementation
git add src/modalAnalysis/

# 2. Python implementation
git add adflow/pyResolventAnalysis.py
git add adflow/__init__.py

# 3. ADflow integration
git add adflow/pyADflow.py

# 4. Documentation
git add doc/resolvent_analysis.md
git add RESOLVENT_IMPLEMENTATION.md
git add RESOLVENT_INTEGRATION_PLAN.md
git add RESOLVENT_JACOBIAN_NOTES.md
git add RESOLVENT_IMPLEMENTATION_COMPLETE.md
git add RESOLVENT_QUICKSTART.md

# 5. Examples
git add examples/resolvent_analysis_example.py

# 6. Migration docs (if needed)
git add MIGRATION_TO_PRIVATE_REPO.md
git add RESOLVENT_BRANCH_SETUP.md

# Check staged files
git status
```

### Step 4: Make Structured Commits

Option A: Single comprehensive commit (simpler)
```bash
git commit -m "Add resolvent analysis implementation

[Detailed message as shown in Quick Setup]"
```

Option B: Multiple logical commits (cleaner history)
```bash
# Commit 1: Fortran API
git add src/modalAnalysis/
git commit -m "Add Fortran API for resolvent analysis

- New modalAnalysis module (not adjoint)
- setupResolventMatrix() to assemble Jacobian
- getResolventMatrixDense() with correct transpose
- Matrix export utilities"

# Commit 2: Python implementation
git add adflow/pyResolventAnalysis.py
git commit -m "Add Python ResolventAnalysis class

- ResolventAnalysis base class
- solveExplicit() with SVD computation
- Complex and real doubled formulations
- Frequency sweep capability"

# Commit 3: ADflow integration
git add adflow/__init__.py adflow/pyADflow.py
git commit -m "Integrate resolvent analysis with ADflow

- Add setupResolventJacobian() method
- Add getJacobianMatrix() method
- Export resolvent classes in __init__"

# Commit 4: Documentation
git add doc/ RESOLVENT_*.md examples/
git commit -m "Add resolvent analysis documentation and examples

- Complete theory and usage guide
- Implementation notes
- Quick start guide
- Working example script"
```

### Step 5: Push to Remote

```bash
# Push branch to private repo
git push -u private feature/resolvent-analysis

# Push tags
git push private --tags
```

## Working on the Branch

### Make Additional Changes

```bash
# Make sure you're on the feature branch
git checkout feature/resolvent-analysis

# Make changes
vim src/modalAnalysis/resolventAPI.F90

# Commit
git add src/modalAnalysis/resolventAPI.F90
git commit -m "Fix: Correct matrix transpose in dense export"

# Push updates
git push private feature/resolvent-analysis
```

### Keep Branch Updated with Master

```bash
# Fetch latest from master
git fetch origin master

# Option 1: Rebase (cleaner history)
git checkout feature/resolvent-analysis
git rebase origin/master

# Option 2: Merge (preserves branch history)
git checkout feature/resolvent-analysis
git merge origin/master

# Push updated branch
git push private feature/resolvent-analysis --force-with-lease  # if rebased
# or
git push private feature/resolvent-analysis  # if merged
```

## Branch Management Strategies

### Strategy 1: Long-lived Feature Branch

Keep the branch separate for ongoing development:

```bash
feature/resolvent-analysis
├── master (stable ADflow)
└── feature branch (resolvent work)

# Periodically sync with master
git checkout feature/resolvent-analysis
git merge master
```

**Pros:**
- Clear separation
- Easy to maintain
- Can develop without affecting master

**Cons:**
- May diverge from master over time
- Need to periodically merge

### Strategy 2: Sub-branches for Components

Create sub-branches for different aspects:

```bash
feature/resolvent-analysis      # Main feature branch
├── feature/resolvent-core      # Core implementation (current work)
├── feature/resolvent-adjoint   # Adjoint derivatives (future)
└── feature/resolvent-matfree   # Matrix-free methods (future)

# Create sub-branch
git checkout feature/resolvent-analysis
git checkout -b feature/resolvent-adjoint

# Merge sub-branch back when done
git checkout feature/resolvent-analysis
git merge feature/resolvent-adjoint
```

### Strategy 3: Release Branch

For production-ready versions:

```bash
# Create release branch
git checkout feature/resolvent-analysis
git checkout -b release/2.13.0-resolvent

# Final testing and documentation
# ... make final tweaks ...

# Tag release
git tag -a v2.13.0-resolvent -m "Resolvent analysis production release"

# Merge to master (when ready)
git checkout master
git merge release/2.13.0-resolvent
```

## Comparison: Branch vs Master

### On Feature Branch ✓

```bash
git checkout feature/resolvent-analysis

# Your work (new files)
src/modalAnalysis/resolventAPI.F90
adflow/pyResolventAnalysis.py
doc/resolvent_analysis.md
...

# Standard ADflow files (unchanged)
src/adjoint/
adflow/pyADflow.py (with your additions)
...
```

### On Master

```bash
git checkout master

# Only standard ADflow
# No resolvent files
ls src/modalAnalysis/
# Directory doesn't exist

ls adflow/pyResolventAnalysis.py
# File not found
```

## Pull Request Workflow

When ready to merge (to your private master or contribute back):

### For Private Repo

```bash
# On GitHub/GitLab web interface:
# 1. Go to your private repo
# 2. Click "Pull Requests" or "Merge Requests"
# 3. Click "New Pull Request"
# 4. Base: master, Compare: feature/resolvent-analysis
# 5. Review changes
# 6. Merge when ready

# Or via command line:
git checkout master
git merge feature/resolvent-analysis
git push private master
```

### For Contributing Back to Public ADflow

```bash
# 1. Fork mdolab/adflow on GitHub
# 2. Push your branch
git remote add upstream git@github.com:mdolab/adflow.git
git push origin feature/resolvent-analysis

# 3. Create Pull Request on GitHub
# Base repository: mdolab/adflow, base: master
# Head repository: YOUR_USERNAME/adflow, compare: feature/resolvent-analysis
```

## Viewing Differences

### See What's Changed

```bash
# Compare branch to master
git diff master..feature/resolvent-analysis --stat

# List new files
git diff master..feature/resolvent-analysis --name-status | grep "^A"

# See specific file changes
git diff master..feature/resolvent-analysis src/modalAnalysis/resolventAPI.F90
```

### Generate Patch

```bash
# Create patch file for sharing
git diff master..feature/resolvent-analysis > resolvent-analysis.patch

# Apply patch elsewhere
git apply resolvent-analysis.patch
```

## Testing on Branch

### Build and Test

```bash
# Make sure you're on the feature branch
git checkout feature/resolvent-analysis

# Build
make clean
make

# Run tests
python -m pytest tests/

# Test resolvent specifically
cd article_resolvent_opt/code
python example.py  # Should still work
```

### Continuous Integration

If you have CI/CD, it will test the branch:

```yaml
# .github/workflows/tests.yml
on:
  push:
    branches: [ master, feature/** ]  # Test feature branches too
  pull_request:
    branches: [ master ]
```

## Switching Between Branches

```bash
# Work on resolvent
git checkout feature/resolvent-analysis
# ... make changes ...
git commit -m "Update resolvent implementation"

# Switch back to master for other work
git checkout master
# ... resolvent files disappear ...

# Switch back to resolvent
git checkout feature/resolvent-analysis
# ... resolvent files reappear ...
```

## Archive/Cleanup

### Archive a Branch

```bash
# Create archive tag before deleting
git tag archive/resolvent-analysis-v1 feature/resolvent-analysis

# Delete local branch (if merged)
git branch -d feature/resolvent-analysis

# Delete remote branch
git push private --delete feature/resolvent-analysis

# Restore from archive if needed
git checkout -b feature/resolvent-analysis archive/resolvent-analysis-v1
```

## Quick Reference

```bash
# Create branch
git checkout -b feature/resolvent-analysis

# Commit work
git add <files>
git commit -m "Description"

# Push branch
git push -u private feature/resolvent-analysis

# Update from master
git merge master

# See changes
git diff master --stat

# Switch branches
git checkout master                    # Go to master
git checkout feature/resolvent-analysis # Back to feature

# Delete branch (if done)
git branch -d feature/resolvent-analysis
```

## Recommended Approach

For your situation, I recommend:

1. **Create feature branch** ✓
   ```bash
   git checkout -b feature/resolvent-analysis
   ```

2. **Commit all resolvent work** ✓
   ```bash
   git add src/modalAnalysis/ adflow/ doc/ examples/ RESOLVENT_*.md
   git commit -m "Add resolvent analysis implementation"
   ```

3. **Push to private repo** ✓
   ```bash
   git remote add private git@github.com:YOUR_USERNAME/adflow-private.git
   git push -u private feature/resolvent-analysis
   git push private master  # Also push master
   ```

4. **Develop on branch** ✓
   - All resolvent work on `feature/resolvent-analysis`
   - Keep `master` clean for syncing with public ADflow
   - Merge to master when ready for production

This gives you maximum flexibility while keeping things organized!
