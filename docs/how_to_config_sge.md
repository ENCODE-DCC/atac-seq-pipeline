How to configure SGE for pipeline
=================================

1. List all parallel environments (PE) on your SGE.
    ```
      $ qconf -spl
    ```

2. If you don't have one then ask your system admin to add a new one with name `shm`.
    ```
      $ sudo qconf -ap shm
    ```

3. Give a large number to `slots` for your PE.
    ```
      $ sudo qconf -mp shm
      pe_name            shm
      slots              999
      ...
    ```

4. List all queues on your SGE.
    ```
      $ qconf -sql
    ```

5. Ask your system admin to connect PE to your queue.
    ```
      $ sudo qconf -mq [QUEUE_NAME]
      ...
      pe_list               make shm
      ...
    ```
