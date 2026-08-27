# Samples

Samples demonstrate supported public APIs with complete error handling. They must not access private structures or use benchmark-only
shortcuts.

`make sample` builds and runs `basic_usage.c` as `bin/basic_usage`.

`make client-sample` builds `client_usage.c` as `bin/client_usage`. Start `geoboltd`, then run it with the token outside the process
argument list:

```bash
GEOBOLT_TOKEN='secret' bin/client_usage 127.0.0.1 7447
```
