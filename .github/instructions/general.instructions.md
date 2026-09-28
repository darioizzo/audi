---
description: "Generic instructions for the Copilot agent working in the audi repository."
applyTo: "**"
---

# General agent instructions

- Always use the conda environment `audi` when running, building, testing, or performing any other code operation in this repository.
- Activate it before running commands in a terminal: `conda activate audi`.
- Tools such as `cmake`, `ninja`, compilers, and Python packages are expected to come from this environment (`$CONDA_PREFIX`).

## Files and git

- Never delete files without first asking the user for permission.
- Never commit or push to git without first asking the user for permission.
- Assume and prefer that the user performs commits, pushes, and file deletions themselves; suggest the commands instead of running them.
