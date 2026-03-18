# Release process

This page documents the current PySME release flow (CI + PyPI).

## Branch and CI model

- Day-to-day development happens on `develop` (and feature branches).
- Pull requests into `develop` or `master` run the main workflow in `.github/workflows/python-app.yml`.
- PR checks are change-aware:
  - Code changes run pytest matrices on Linux and macOS.
  - Docs-only changes run the docs build job.

## What triggers a real release

Only a `v*` tag push triggers package build/publish jobs, and release tags are gated to commits on `origin/master`.

This means:
- Tag on `master` commit: release pipeline runs.
- Tag on non-`master` commit: workflow is skipped by gate.

## Standard release steps

1. Merge validated `develop` changes into `master`.
2. Confirm `master` CI is green.
3. Create and push a SemVer tag on `master` (for example `v1.0.0`).
4. Wait for GitHub Actions release jobs:
   - Build source distribution
   - Build wheels
   - Publish to PyPI
   - Create GitHub Release
5. Verify published package:
   - `pip install pysme-astro==<version>`
   - `python -c "import pysme; print(pysme.__version__)"`

## Repository settings to keep aligned

- In branch protection/rulesets, keep required checks aligned with current job names.
- For PR merging, require test jobs and docs build as appropriate.
- Do not require release-only jobs (`Build source distribution`, wheel build, `Publish to PyPI`, `Create GitHub Release`) for PRs.

