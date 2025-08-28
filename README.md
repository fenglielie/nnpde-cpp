# NNPDE-CPP

For the 1D Burgers equation problem, examples based on the following high-precision numerical method are provided:

- RKDG
- FV-WENO
- FD-WENO


C++23 is required, with the following feature support:

- `std::format` (optional)
- `deducing this` (required for some examples, skip if unsupported)


Test environment

|   Platform    | Compiler |  Version   | Required flags |
| :-----------: | :------: | :--------: | :------------: |
| Ubuntu-latest |   GCC    |     13     |  `-std=c++23`  |
| Ubuntu-latest |  Clang   |     18     |  `-std=c++23`  |
| Windows-2022  |   MSVC   | 14(VS2022) |  `/std:c++23`  |
