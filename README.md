# CS124Fun

This is a project from winter quarter 2020. I adapted the EM algorithm I'd
written when I took the class in 2019. Out of curiosity, I ported the Python
code to C++ (and added slightly better command line parsing). To no one's
surprise, the C++ code ran dramatically faster (~17x faster single-threaded
and far faster multi-threaded).

The Python version (ip3.py) can be invoked as follows:
```
python3 ip3.py assignment/example_data_1
```

The C++ versions can be run as follows:
```
./ip6 assignment/example_data_1 8 10 8
```

The Makefile currently compiles ip6.cpp; ip5.cpp is included here primarily as
a point of comparison. The newer version utilizes shared memory rather than
writing partial work to files and reading them back in needlessly. Most of the
transpose() calls were also removed.

Brief, non-definitive performance comparisons:
```
  (missing/heterozygous: 14, overlap: 8, file: assignment/h1000):
    ip3.py:
      real  0m51.663s
      user  0m51.350s
      sys   0m0.624s
    ip5.cpp (single-threaded):
      real  0m3.313s
      user  0m3.304s
      sys   0m0.009s
    ip6.cpp (single-threaded):
      real  0m2.966s
      user  0m2.953s
      sys   0m0.012s

  (missing/heterozygous: 14, overlap: 8, file: assignment/example_data_1):
    ip3.py:
      real  37m22.822s
      user  37m7.940s
      sys   0m19.197s
    ip5.cpp (8 threads):
      real  1m23.682s
      user  10m38.838s
      sys   0m0.480s
    ip6.cpp (8 threads):
      real  1m19.145s
      user  10m10.054s
      sys   0m0.212s
```

Accuracy comparisons for various hyperparameters:
```
  (missing/heterozygous: 14, overlap: 8, file: assignment/example_data_1):
    ip6.cpp (8 threads): 0.9117316
  (missing/heterozygous: 10, overlap: 8, file: assignment/example_data_1):
    ip6.cpp (8 threads): 0.9048697
```

Unless you really need the slight accuracy bump, setting missing/heterozygous
to 10 is a much better use of time (runs in <5 sec on my machine).
