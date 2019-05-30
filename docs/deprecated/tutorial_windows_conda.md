Tutorial for Windows 10 Pro/Home computers with Conda
===============================================

1. [Enable Developer mode](https://www.howtogeek.com/292914/what-is-developer-mode-in-windows-10/).

2. [Turn on `Windows Subsystem for Linux`](https://www.pcworld.com/article/3106463/windows/how-to-get-bash-on-windows-10-with-the-anniversary-update.html)

3. Sign up for a Microsoft account.

4. [Install Ubuntu from Microsoft Store](https://www.microsoft.com/en-us/p/ubuntu/9nblggh4msv6?activetab=pivot:overviewtab).

5. Reboot your computer.

6. Disable/turn off any VPN clients (Cisco VPN, ...) and antivirus softwares (AVG, AVAST, Kaspersky, ...) on your system. It blocks any IPv6 network traffic from `apt` in Ubuntu. You can re-enable antivirus softwares later after installation.

7. Install basic packages for developers.
```
$ sudo apt-get update
$ sudo apt-get install build-essential python python3 git wget default-jre
```

8. Proceed to [Local system with `Conda`](tutorial_local_conda.md).