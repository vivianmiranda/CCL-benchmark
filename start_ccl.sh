#!/bin/bash
# Activate the private PR environment after its Conda base.
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
  echo "Use: source start_ccl.sh" >&2
  exit 1
fi
if [[ -z "${CONDA_PREFIX:-}" || -n "${VIRTUAL_ENV:-}" ||
      -n "${ROOTDIR:-}" ]]; then
  echo "Activate only the ccldev Conda base in a fresh shell." >&2
  return 1
fi
source "$(dirname "${BASH_SOURCE[0]}")/.local/bin/activate" || return 1
return 0
