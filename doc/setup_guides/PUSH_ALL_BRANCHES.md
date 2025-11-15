# Push All Branches to Private Repository

## Quick Commands

### Push All Branches at Once

```bash
cd /home/sicheng/repo/adflow

# Push all local branches
git push private --all

# Push all tags
git push private --tags
```

### Or Push Specific Remote Branches

```bash
# First, create local tracking branches for the ones you want
git checkout -b ANK_modifications origin/ANK_modifications
git checkout -b Dev origin/Dev
git checkout -b Matrix_Free_ANK origin/Matrix_Free_ANK
# ... etc for any branches you want

# Then push all local branches
git push private --all
```

## Mirror Everything (Recommended)

To create an exact mirror with ALL branches from origin:

```bash
# Method 1: Using git push --mirror (copies everything)
cd /home/sicheng/repo/adflow

# Push mirror (includes all branches and tags)
git push private --mirror
```

**⚠️ Warning**: `--mirror` will make private an exact copy, including all remote-tracking branches.

## Selective Branch Migration

If you only want specific branches:

### Option 1: Mirror Then Add Resolvent

```bash
# 1. Mirror all public branches
git push private --mirror

# 2. Add your resolvent branch
git checkout feature/resolvent-analysis
git push -u private feature/resolvent-analysis
```

### Option 2: Push Only What You Need

```bash
# Push just master and your feature branch
git push private master
git push private feature/resolvent-analysis
git push private --tags

# Later, if you need other branches, push them individually
git push private ANK_modifications
git push private Dev
# etc.
```

## Recommended Approach

For your situation, I recommend:

```bash
cd /home/sicheng/repo/adflow

# 1. First, fix the authentication (HTTPS or SSH as discussed)

# 2. Mirror everything to preserve all ADflow branches
git push private --mirror

# 3. Your resolvent branch is already pushed, verify it
git ls-remote private
```

## Verify All Branches Were Pushed

```bash
# List remote branches on private
git ls-remote --heads private

# Or clone in a new location to check
cd /tmp
git clone https://github.com/YOUR_USERNAME/adflow-private.git test-mirror
cd test-mirror
git branch -a
# Should see all branches
```

## Updating Branches Later

```bash
# Pull updates from public origin
git fetch origin

# Push updates to private
git push private --all
git push private --tags

# Or for specific branch
git fetch origin master:master
git push private master
```

## Complete Setup Script

Here's a complete script to set up your private repo with all branches:

```bash
#!/bin/bash

cd /home/sicheng/repo/adflow

# 1. Create private repo on GitHub first (via web interface)
#    https://github.com/new
#    Name: adflow-private
#    Private: Yes

# 2. Add remote (choose HTTPS or SSH)
# For HTTPS:
git remote add private https://github.com/YOUR_USERNAME/adflow-private.git
# For SSH:
# git remote add private git@github.com:YOUR_USERNAME/adflow-private.git

# 3. Push everything as a mirror
git push private --mirror

# 4. Verify
echo "Verifying remote branches..."
git ls-remote --heads private | head -20

echo "Done! All branches have been pushed to private repository."
echo ""
echo "Branches on private repo:"
git ls-remote --heads private | awk '{print $2}' | sed 's/refs\/heads\///'
```

## What Gets Pushed

With `git push private --mirror`, you get:

```
private repository will contain:
├── master
├── feature/resolvent-analysis (your work)
├── ANK_modifications
├── Dev
├── Matrix_Free_ANK
├── buffet
├── force_coupling
├── jac
├── low_speed
├── matrixfree-adjoint
├── reverse_adjoint
├── timespectral
├── timespectral_adjoint
└── ... (all other branches)
└── All tags
```

## Maintaining Sync

### Regular Updates

```bash
# Weekly or as needed:

# 1. Fetch from public ADflow
git fetch origin

# 2. Update all branches
git push private --all
git push private --tags

# Or use mirror again
git push private --mirror
```

### Two-Way Sync Strategy

```bash
# From public to private (updates)
git fetch origin
git push private --all

# Your work stays on feature branch
git checkout feature/resolvent-analysis
# ... your changes ...
git push private feature/resolvent-analysis
```

## Branch Organization

After pushing all branches, your private repo structure:

```
Private Repository (adflow-private)
├── Public ADflow Branches
│   ├── master                  ← Synced with mdolab/adflow
│   ├── Dev
│   ├── Matrix_Free_ANK
│   └── ... (all others)
│
└── Your Development
    └── feature/resolvent-analysis  ← Your resolvent work
```

## FAQ

### Q: Will this copy all the history?
**A**: Yes! `--mirror` copies complete history for all branches.

### Q: Can I still update from public ADflow?
**A**: Yes! Use `git fetch origin` then `git push private --all`

### Q: What if I don't want all branches?
**A**: Just push selectively:
```bash
git push private master
git push private feature/resolvent-analysis
```

### Q: How much space will this use?
**A**: Same as your local repo. All branches share the same object database, so it's efficient.

## Recommended Commands for Your Situation

```bash
# After fixing SSH/HTTPS authentication:

cd /home/sicheng/repo/adflow

# Option 1: Mirror everything (recommended)
git push private --mirror

# Option 2: Just your work + master
git push private master
git push private feature/resolvent-analysis
git push private --tags

# Verify
git ls-remote private
```

I recommend **Option 1** (mirror) so you have a complete backup of the entire ADflow repository plus your resolvent work.
