---
name: Bug report
about: Create a report to help us improve
title: ''
labels: ''
assignees: ''

---

## **Describe the bug**
A clear and concise description of what the problem is.

## **OS/Platform**
- OS/Platform: [e.g. Ubuntu 18.04, Google Cloud, Stanford Sherlock/SCG cluster, ...]
- Conda version: If you used Conda (`$ conda --version`).
- Pipeline version: [e.g. v1.8.0]
- Caper version: [e.g. v1.2.0]

## **Caper configuration file**
Paste contents of `~/.caper/default.conf`.
```ini
PASTE CAPER CONF CONTENTS HERE
```

## **Input JSON file**
Paste contents of your input JSON file.
```json
PASTE INPUT JSON CONTENTS HERE
```

## **Troubleshooting result**

If you ran `caper run` without Caper server then Caper automatically runs a troubleshooter for failed workflows. Find troubleshooting result in the bottom of Caper's screen log.

If you ran `caper submit` with a running Caper server then first find your workflow ID (1st column) with `caper list` and run `caper debug [WORKFLOW_ID]`.

Paste troubleshooting result.
```
PASTE TROUBLESHOOTING RESULT HERE
```
