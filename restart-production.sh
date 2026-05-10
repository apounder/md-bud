#!/bin/bash
# Summarize total frames across trajectory files, then move all .nc files into a folder named frames-<TOTAL>

# Fail fast on pipeline errors, treat unmatched globs as empty
shopt -s nullglob
set -o pipefail

echo "Selecting most recent PRMTOP..."
PT=$(ls -t -- *prmtop 2>/dev/null | head -n 1)
if [[ -z "$PT" ]]; then
  echo "Error: No *.prmtop file found in the current directory." >&2
  exit 1
fi

echo "Using $PT as the prmtop"
echo "This may take a moment... please be patient."

# Collect trajectory files (nc and mdcrd)
traj_files=( *.nc *.mdcrd )
if (( ${#traj_files[@]} == 0 )); then
  echo "Error: No trajectory files (*.nc or *.mdcrd) found in the current directory." >&2
  exit 1
fi

total_frames=0
file_count=0

for file in "${traj_files[@]}"; do
  # Extract the last integer reported by cpptraj -tl as the frame count for this file
  frames=$(cpptraj -p "$PT" -y "$file" -tl 2>/dev/null | grep -oE '[0-9]+' | tail -n 1)
  if [[ -n "$frames" && "$frames" =~ ^[0-9]+$ ]]; then
    total_frames=$(( total_frames + frames ))
    file_count=$(( file_count + 1 ))
  else
    echo "Warning: Could not determine frame count for '$file' — skipping its contribution." >&2
  fi
done

echo "${total_frames} <- Total frames over all ${file_count} simulations"

# Create destination directory named frames-<TOTAL> and move all .nc files there
dest_dir="frames-${total_frames}"
if [[ -d "$dest_dir" ]]; then
  echo "Note: Destination directory '$dest_dir' already exists."
else
  mkdir -p "$dest_dir"
fi

nc_files=( *.nc )
if (( ${#nc_files[@]} > 0 )); then
  echo "Moving ${#nc_files[@]} .nc file(s) to '$dest_dir/'..."
  mv -- *.nc "$dest_dir"/
else
  echo "No .nc files found to move."
fi
