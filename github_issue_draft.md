Title: Action required: re-clone repo after Git LFS migration

---

## What happened

The repository history was rewritten to migrate all `.rds` files to [Git Large File Storage (LFS)](https://git-lfs.github.com/). This was necessary because several precomputed result files exceeded GitHub's 100 MB file size limit, preventing pushes.

## What you need to do

**You must re-clone the repository.** Because the commit history was rewritten, your existing local clone will have conflicts that cannot be cleanly resolved with `git pull`.

```bash
# 1. Delete your old clone (back up any local uncommitted work first!)
rm -rf DeSurv-paper

# 2. Clone fresh
git clone git@github.com:rashidlab/DeSurv-paper.git

# 3. Ensure Git LFS is installed (one-time setup)
git lfs install
```

### Prerequisites

You need Git LFS installed on your machine. If you don't have it:
- **macOS**: `brew install git-lfs`
- **Ubuntu/Debian**: `sudo apt install git-lfs`
- **Conda**: `conda install -c conda-forge git-lfs`

## Why this was needed

73 `.rds` files (data and precomputed results) were being stored as regular Git objects. Two files exceeded GitHub's hard 100 MB limit, blocking all pushes. The LFS migration converts these to lightweight pointers in the repo, with the actual file contents stored in GitHub's LFS backend.

This is a one-time disruption. Going forward, all `.rds` files will be automatically handled by LFS.
