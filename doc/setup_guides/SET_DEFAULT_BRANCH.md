# Set Default Branch to Master on GitHub

## Quick Steps

1. Go to your private repository: https://github.com/SichengHe/adflow-private

2. Click **Settings** (top right, near the tabs)

3. In the left sidebar, click **Branches** (under "Code and automation")

4. Under "Default branch", you'll see the current default branch

5. Click the **switch icon** (⇄) or pencil icon next to the branch name

6. Select **master** from the dropdown

7. Click **Update**

8. Confirm by clicking **I understand, update the default branch**

Done! ✅

## Why This Matters

The default branch is what people see when they visit your repo and what gets checked out when someone clones.

**Current default**: Probably `feature/resolvent-analysis` (most recently pushed)
**Want**: `master` (standard main branch)

## Verify It Worked

After changing:

1. Visit: https://github.com/SichengHe/adflow-private
2. You should see the master branch files (standard ADflow)
3. To see resolvent work, switch to: feature/resolvent-analysis

## Alternative: Use GitHub CLI (if installed)

```bash
# If you have gh CLI installed
gh repo edit SichengHe/adflow-private --default-branch master
```

## Or Use GitHub API

```bash
# Using curl (replace YOUR_TOKEN with GitHub Personal Access Token)
curl -X PATCH \
  -H "Accept: application/vnd.github.v3+json" \
  -H "Authorization: token YOUR_TOKEN" \
  https://api.github.com/repos/SichengHe/adflow-private \
  -d '{"default_branch":"master"}'
```

## Recommended

Just use the web interface - it's the easiest and safest way!

**Direct link**: https://github.com/SichengHe/adflow-private/settings/branches
