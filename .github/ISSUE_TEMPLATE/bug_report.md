---
name: Bug report
about: Create a report to help us improve
title: ''
labels: ''
assignees: ''

---

**Describe the bug**
A clear and concise description of what the problem is.

**OS/Platform and dependencies**
- OS or Platform: [e.g. Ubuntu 16.04, Google Cloud, Stanford Sherlock/SCG cluster, ...]
- Cromwell/dxWDL version: [e.g. `cromwell-34.jar`, `dxWDL-78.jar`]
- Conda version: If you have used Conda (`$ conda --version`).
- singularity version: If you have used singularity (`$ singularity --version`).

**Attach error logs**
For Cromwell users only.
1) Move to your working directory where you ran a pipeline. You should be able to find a directory named `cromwell-executions/` which includes all outputs and logs for debugging.

2) Run the following command line to print all non-empty STDERR outputs. This will be greatly helpful for developers to figure out the problem. Copy-paste its output to the issue page.
```
$ find -name stderr -not -empty | xargs tail -n +1
```

3) (OPTIONAL) Run the following command to collect all logs. For developer's convenience, please add `[ISSUE_ID]` to the name of the tar ball file. This command will generate a tar ball including all debugging information. Post an issue with the tar ball (`.tar.gz`) attached.
```
$ find . -type f -name 'stdout' -or -name 'stderr' -or -name 'script' -or \
-name '*.qc' -or -name '*.txt' -or -name '*.log' -or -name '*.png' -or -name '*.pdf' \
| xargs tar -zcvf debug_[ISSUE_ID].tar.gz
```
