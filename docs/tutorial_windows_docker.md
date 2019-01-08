Tutorial for Windows 10 Pro computers with docker
===============================================

## Windows 10 Home

This instruction is for Windows 10 Pro only. It is a bit tricky to install `docker` for Windows 10 Home, which does not support virtualization (Hyper-V). For users with Windows 10 Home, read through [this](https://gist.github.com/strarsis/44ded0d254066d9cb125ebbb04650d6c) to install `docker toolbox` instead. You may need to configure environment variables (`DOCKER_TLS_VERIFY` and `DOCKER_HOST`) and docker communication port (2376) correctly with TSL and copy encryption keys into Ubuntu. Here is an example of `~/.bashrc` and encryption keys inside Ubuntu.
```
export DOCKER_HOST="tcp://192.168.99.100:2376"
export DOCKER_TLS_VERIFY=1
```
```
$ ls ~/.docker/
ca-key.pem  ca.pem  cert.pem  config.json  key.pem
```

## Windows 10 Pro

1. [Enable Developer mode](https://www.howtogeek.com/292914/what-is-developer-mode-in-windows-10/).

2. [Turn on `Windows Subsystem for Linux`](https://www.pcworld.com/article/3106463/windows/how-to-get-bash-on-windows-10-with-the-anniversary-update.html)

3. Sign up for a Microsoft account.

4. [Install Ubuntu from Microsoft Store](https://www.microsoft.com/en-us/p/ubuntu/9nblggh4msv6?activetab=pivot:overviewtab).

5. Reboot your computer.

6. Disable/turn off any VPN clients (Cisco VPN, ...) and antivirus softwares (AVG, AVAST, Kaspersky, ...) on your system. It blocks any IPv6 network traffic from `apt` in Ubuntu. You can re-enable antivirus softwares later after installation.

7. Sign up for [docker hub](https://hub.docker.com/).

8. [Install `docker` on your Windows 10 Pro](https://hub.docker.com/editions/community/docker-ce-desktop-windows).

9. Press a Window button and type `ubuntu`. A BASH terminal will pop up. Press Enter to install Ubuntu on your Windows system.

10. Install basic packages for developers.
```
$ sudo apt-get update
$ sudo apt-get install build-essential python python3 git wget default-jre
```

11. [Follow this instruction to install docker inside Ubuntu](https://medium.com/@sebagomez/installing-the-docker-client-on-ubuntus-windows-subsystem-for-linux-612b392a44c4). If you see any `Connection refused` errors during `apt` installation then disable VPN and antivirus softwares.

12. Remove `~/.bash_profile` which was created during previous step. Run the following command line.
```bash
$ echo "export DOCKER_HOST=localhost:2375" >> ~/.bash_rc
```

13. Close BASH terminal and re-open it.

14. Check if docker works by running `hello-world`.
```bash
$ docker run hello-world
```

15. Proceed to [Local system with `docker`](tutorial_local_docker.md).