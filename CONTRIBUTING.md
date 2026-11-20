# Contributing to PySME

Thank you for your interest in improving PySME. Contributions of all sizes are welcome: bug reports, fixes, documentation, tests and new features.

## Where to open pull requests

**Please open pull requests against the `develop` branch, not `master`.**
`master` only receives release merges from `develop`. A pull request opened against `master` is retargeted to `develop` automatically.

Urgent fixes to a released version are the only exception: maintainers open them from a `hotfix/...` branch directly against `master`.

## Before larger changes

For a new feature or a larger change, it is recommended to open an issue first. 

## Suggested workflow

1. Fork the repository, or create a branch if you have write access.
2. Use a descriptive branch name, for example `fix/short-description`, `feature/short-description` or `docs/short-description`.
3. Keep each pull request focused on one change.
4. Tests for new features or bug fixes are welcome but not required. Please do not change existing tests or their reference values to make them pass.
5. Add a line to `changelog.md` under the unreleased section.
6. Open the pull request against `develop` and describe what changed and why.

## Pull request titles

A short prefix in the title helps us sort changes and write release notes, for example:

- `fix: ...` for a bug fix
- `feat: ...` for a new feature
- `docs: ...` for documentation
- `ci: ...` for build, workflow and other maintenance changes

This is a recommendation, not a requirement.

## Behaviour changes

- New options should be **off by default**, and with the option off, results should stay identical to the previous version.
- If a change alters results on the default path (for example a numerical fix), please say so explicitly in the pull request and show the size of the change.

## Questions

If you are unsure about anything, open an issue or a draft pull request and ask.

