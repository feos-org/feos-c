# Compilation

Do this once to generate dynamic link to library.

## LINUX

```
cargo build
cd src/c
ln -s ../../target/debug/libfeos_c.so .
```

To compile the C code from within `src/c` directory:

```
gcc test.c -L. -l:libfeos_c.so -o test
``` 

Run with:

```
LD_LIBRARY_PATH=.:$LD_LIBRARY_PATH  ./test
```

## MAC

```
cargo build
cd src/c
ln -s ../../target/debug/libfeos_c.dylib .
```

To compile the C code from within `src/c` directory:

```
clang test.c -o test libfeos_c.dylib -Wl,-rpath,/global-path-to/src/c
``` 

(adjust global path)

Run with:

```
./test
```
