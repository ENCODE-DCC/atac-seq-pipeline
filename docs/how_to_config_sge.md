How to configure SGE for pipeline
=================================

1. List all parallel environments (PE) on your SGE.
    ```bash
      $ qconf -spl
    ```

2. If you don't have one then ask your system admin to add a new one with name `shm`.
    ```bash
      $ sudo qconf -ap shm
    ```

3. Give a large number to `slots` for your PE.
    ```bash
      $ sudo qconf -mp shm
      pe_name            shm
      slots              999
      ...
    ```

4. List all queues on your SGE.
    ```bash
      $ qconf -sql
    ```

5. Ask your system admin to connect PE to your queue.
    ```bash
      $ sudo qconf -mq [QUEUE_NAME]
      ...
      pe_list               make shm
      ...
    ```
