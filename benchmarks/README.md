# Benchmarks

`benchmark.py` writes BED lines at compression level 6, then reads them back, and compares pybgzf with the standard library's `gzip` and, if installed, htslib's `bgzip` and `tabix`.

```console
uv run python benchmarks/benchmark.py 200
```

The argument is the size of the data in megabytes, 200 by default.

Results for 200 MB on an Apple M3 Max:

| Command                             | Seconds |  MB/s |
|-------------------------------------|--------:|------:|
| `gzip` (stdlib)                     |    4.80 |    42 |
| pybgzf, 1 thread, `.tbi`            |    1.84 |   109 |
| pybgzf, 4 threads, `.tbi`           |    0.43 |   466 |
| pybgzf, 16 threads, `.tbi`          |    0.31 |   641 |
| `bgzip -@1` then `tabix`            |    3.01 |    66 |
| `bgzip -@4` then `tabix`            |    1.31 |   153 |
| `bgzip -@16` then `tabix`           |    0.90 |   222 |
| lines from `gzip.open` (stdlib)     |    0.66 |   303 |
| lines from `reader`, 1 thread       |    0.44 |   456 |
| lines from `reader`, 4 threads      |    0.30 |   671 |
| bytes from `BgzfReader`, 16 threads |    0.03 |  7608 |
| `bgzip -d -@16`                     |    0.03 |  7652 |
