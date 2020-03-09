# morex_reference

Preparation and processing steps related to all versions of the morex reference genome.

All detailed methods/documentation are within each of the subdirectories `morex_v1`, `morex_v2`, etc.

## Important instructions to get all submodules within this repository

To clone repository, run both of the following commands:

```bash
# Initial clone
git clone https://github.com/MorrellLAB/morex_reference.git

# Go into Github repo
cd morex_reference

# After initial clone, pull submodules
# You only need to run this once
git submodule update --init --recursive
```

After these steps, you can commit and push changes as you would normally.

If you would like to update each submodule so they are at the latest version of the original repository, run the following:

```bash
# Update each submodule
git submodule foreach git pull origin master

# Then commit and push these changes (to the parent repo) so your
# repo containing submodules stays up to date
git commit -am "update submodule directories"
git push
```
