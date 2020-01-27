---
name: Bug report
about: Create a report to help us improve
title: ''
labels: ''
assignees: ''

---

**Describe the bug**
A clear and concise description of what the problem is.

**OS/Platform**
- OS/Platform: [e.g. Ubuntu 16.04, Google Cloud, Stanford Sherlock/SCG cluster, ...]
- Conda version: If you used Conda (`$ conda --version`).
- Pipeline version: [e.g. v1.5.3]
- Caper version: [e.g. v0.6.0]

**Caper configuration file**
Paste contents of `~/.caper/default.conf`.

**Input JSON file**
Paste contents of your input JSON file.

**Error log**
Caper automatically runs a troubleshooter for failed workflows. If it doesn't then get a `WORKFLOW_ID` of your failed workflow with `caper list` or directly use a `metadata.json` file on Caper's output directory.
```
$ caper debug [WORKFLOW_ID_OR_METADATA_JSON_FILE]
```
