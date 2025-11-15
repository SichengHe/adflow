# Fix GitHub SSH Permission Issue

## The Error

```
git@github.com: Permission denied (publickey).
fatal: Could not read from remote repository.
```

This means GitHub doesn't recognize your SSH key.

## Quick Fix

### Option 1: Use HTTPS Instead (Easiest)

```bash
# Remove the SSH remote
git remote remove private

# Add HTTPS remote instead
git remote add private https://github.com/YOUR_USERNAME/adflow-private.git

# Push (will ask for GitHub credentials)
git push -u private feature/resolvent-analysis

# Or use GitHub Personal Access Token
# Settings > Developer settings > Personal access tokens > Generate new token
# Use token as password when prompted
```

### Option 2: Set Up SSH Keys (More Secure)

#### Step 1: Check if you have SSH keys

```bash
ls -la ~/.ssh
# Look for id_rsa.pub, id_ed25519.pub, or similar
```

#### Step 2: Generate SSH key (if needed)

```bash
# Generate new SSH key
ssh-keygen -t ed25519 -C "your_email@example.com"

# Press Enter for default location (~/.ssh/id_ed25519)
# Enter passphrase (or leave empty)

# Start SSH agent
eval "$(ssh-agent -s)"

# Add key to agent
ssh-add ~/.ssh/id_ed25519
```

#### Step 3: Copy public key to clipboard

```bash
# Display public key
cat ~/.ssh/id_ed25519.pub

# Or copy to clipboard (Ubuntu/Linux)
cat ~/.ssh/id_ed25519.pub | xclip -selection clipboard
# If xclip not installed: sudo apt install xclip
```

#### Step 4: Add to GitHub

1. Go to GitHub.com
2. Click your profile picture (top right) → **Settings**
3. Click **SSH and GPG keys** (left sidebar)
4. Click **New SSH key**
5. Title: "ADflow Dev Machine" (or whatever you want)
6. Paste your public key
7. Click **Add SSH key**

#### Step 5: Test connection

```bash
ssh -T git@github.com
# Should see: "Hi YOUR_USERNAME! You've successfully authenticated..."
```

#### Step 6: Push to private repo

```bash
git push -u private feature/resolvent-analysis
```

## Detailed Walkthrough

### Check Current Setup

```bash
# Check what remotes you have
git remote -v

# Check if private repo exists on GitHub
# Visit: https://github.com/YOUR_USERNAME/adflow-private
# If not, create it first!
```

### Create Private Repo First (If Not Done)

1. Go to https://github.com/new
2. Repository name: `adflow-private`
3. Description: "ADflow with resolvent analysis (private)"
4. Visibility: **Private**
5. Don't initialize with README
6. Click **Create repository**

### Choose Authentication Method

#### Method A: HTTPS with Personal Access Token (Recommended for Beginners)

```bash
# 1. Create Personal Access Token on GitHub
#    Settings > Developer settings > Personal access tokens > Tokens (classic)
#    Click "Generate new token (classic)"
#    Select scopes: repo (all)
#    Click Generate token
#    COPY THE TOKEN (you won't see it again!)

# 2. Update remote to use HTTPS
git remote set-url private https://github.com/YOUR_USERNAME/adflow-private.git

# 3. Push (use token as password)
git push -u private feature/resolvent-analysis
# Username: YOUR_USERNAME
# Password: [paste your token]

# 4. Store credentials (optional, so you don't have to enter token every time)
git config --global credential.helper store
# Next time you push, it will save the credentials
```

#### Method B: SSH (More Secure, No Passwords)

```bash
# Follow Steps 1-6 above to set up SSH keys

# Make sure remote uses SSH format
git remote set-url private git@github.com:YOUR_USERNAME/adflow-private.git

# Test
ssh -T git@github.com

# Push
git push -u private feature/resolvent-analysis
```

## Common Issues

### Issue 1: Repository Doesn't Exist

```bash
# Error: Repository not found
# Solution: Create the repo on GitHub first (see above)
```

### Issue 2: Wrong Username

```bash
# Make sure YOUR_USERNAME is correct
git remote -v
# Check the URL is correct

# Fix if wrong
git remote set-url private https://github.com/CORRECT_USERNAME/adflow-private.git
```

### Issue 3: SSH Key Not Recognized

```bash
# Make sure key is added to ssh-agent
ssh-add -l  # List keys

# If empty, add it
ssh-add ~/.ssh/id_ed25519

# Test GitHub connection
ssh -T git@github.com
```

### Issue 4: Multiple SSH Keys

```bash
# If you have multiple keys, create SSH config
vim ~/.ssh/config

# Add this:
Host github.com
    HostName github.com
    User git
    IdentityFile ~/.ssh/id_ed25519
    IdentitiesOnly yes
```

## Step-by-Step Solution

### Complete Setup (HTTPS Method)

```bash
cd /home/sicheng/repo/adflow

# 1. Create Personal Access Token on GitHub
#    Go to: https://github.com/settings/tokens
#    Generate new token (classic)
#    Select: repo (full control)
#    Click: Generate token
#    COPY IT!

# 2. Create private repository on GitHub
#    Go to: https://github.com/new
#    Name: adflow-private
#    Private: Yes
#    Create repository

# 3. Update remote to HTTPS
git remote remove private  # If exists
git remote add private https://github.com/YOUR_USERNAME/adflow-private.git

# 4. Push with token
git push -u private feature/resolvent-analysis
# Username: YOUR_GITHUB_USERNAME
# Password: [paste your token]

# 5. Push master too
git push private master

# 6. Push tags
git push private --tags

# Done! ✓
```

### Complete Setup (SSH Method)

```bash
cd /home/sicheng/repo/adflow

# 1. Generate SSH key
ssh-keygen -t ed25519 -C "your_email@example.com"
# Press Enter for defaults
# Enter passphrase (optional)

# 2. Start SSH agent and add key
eval "$(ssh-agent -s)"
ssh-add ~/.ssh/id_ed25519

# 3. Copy public key
cat ~/.ssh/id_ed25519.pub
# Copy the output

# 4. Add to GitHub
#    Go to: https://github.com/settings/keys
#    Click: New SSH key
#    Paste key and save

# 5. Test connection
ssh -T git@github.com
# Should see: "Hi YOUR_USERNAME! You've successfully authenticated"

# 6. Create repository on GitHub
#    Go to: https://github.com/new
#    Name: adflow-private
#    Private: Yes

# 7. Add remote and push
git remote add private git@github.com:YOUR_USERNAME/adflow-private.git
git push -u private feature/resolvent-analysis
git push private master --tags

# Done! ✓
```

## My Recommendation

**For you, I recommend HTTPS with Personal Access Token:**

1. It's simpler to set up
2. No SSH configuration needed
3. Works immediately
4. Just as secure with a strong token

Here's the exact commands:

```bash
# Create token at: https://github.com/settings/tokens
# Then:

cd /home/sicheng/repo/adflow
git remote remove private
git remote add private https://github.com/YOUR_USERNAME/adflow-private.git
git push -u private feature/resolvent-analysis
# Enter username and token when prompted
```

## Verify It Worked

```bash
# Check remote
git remote -v

# Check branch is pushed
git ls-remote private

# Visit on GitHub
# https://github.com/YOUR_USERNAME/adflow-private
```

## Next Steps After Success

```bash
# Push master
git push private master

# Push tags
git push private --tags

# Set up credential storage (optional)
git config --global credential.helper store
```

Let me know which method you prefer and I can help with any issues!
