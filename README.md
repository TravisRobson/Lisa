
## Best Practices

- All header files must be documented fully so that a user of the functions and data types can understand what is being done, what inputs and outputs are, and what assumptions are made.
- Snake case is preferred.
- Forward declare inside of header files whenever possible.
- In *.c/cpp files be sure to include the same named header file first.
- Files shall not include header files they don't need. However, header files themselves must include all other header files necessary to fully parse itself.


## Commands

To run the tests and view HTML output.

```bash
$ ./bin/tests_lisa --gtest_output=xml:test_detail.xml
$ xsltproc src/tests/gtest_to_html.xslt test_detail.xml > test_detail.html
```

## Todo

- Doxygen documentation (i.e. doc Makefile target)
- Replicate functionality of old Galacy code. By that I mean the unresolved galaxy + instrument noise curve.
- Profile the code.
- Assess whether I have the bad galaxy file (i.e. the file that Tyson figured out had many duplicates of the same GB).
- Galaxy_LEGACY.c is the first exectuable I need to reproduce the functionality for.
- I have created an output file by running the LEGACY code. I need to be able to produce it in a binary equivalent manner. This will serve as a regression test that I keep around forever basically. I might end up changing it as I modify things (e.g. using the more modern noise curve approximation).
- Create a version number that can be obtained from the code.
- Try to improve the compile time. Currently, it's about 4 seconds to make and install.
