# Linear-Programming-Problem

## Overview
Linear-Programming-Problem is a C++ console application that demonstrates solving linear optimization tasks using simplex-based methods.

The program includes two predefined scenarios:
- Simplex method flow
- Gomory (integer-focused) flow using additional dual-simplex style iterations

## Why This Project Matters
- Demonstrates core optimization workflow step by step in a terminal application
- Prints simplex tables at each iteration for learning and debugging
- Provides a compact educational implementation in one source file

## Quick Start
1. Open `Linear Programming Problem/Linear Programming Problem.sln` in Visual Studio 2022.
2. Build the solution.
3. Run the program in a terminal.
4. Enter:
	- `0` for simplex scenario
	- `1` for Gomory scenario

## Features
- Iterative simplex table printing
- Pivot selection and row operations
- Optional Gomory-style integer refinement branch
- Final variable and objective output (`F` value)

## Technology Stack
- Language: C++
- Project type: Visual Studio C++ Console Application (`.sln` / `.vcxproj`)
- Toolset: v143 (Visual Studio 2022)
- Target platforms: Win32 and x64

## Prerequisites
### System Requirements
- Windows
- Visual Studio 2022 with Desktop Development for C++

### Dependencies
- Standard C++ library only
- No external solver library

## Build Instructions
### Visual Studio (recommended)
1. Open `Linear Programming Problem/Linear Programming Problem.sln`.
2. Select configuration (`Debug` or `Release`) and platform (`x64` recommended).
3. Build and run.

### Command line
Use `msbuild` from a Visual Studio Developer Command Prompt.

## Running the Application
When started, the program asks:

```text
You need to use Gomori method? 1 \ 0
```

Input options:
- `0`: run simplex demonstration on predefined matrix
- `1`: run Gomory-related flow on predefined matrix

The app prints intermediate tables, pivot values, iteration markers, and final results.

## Project Structure
- `Linear Programming Problem/Linear Programming Problem.sln`: solution file
- `Linear Programming Problem/Linear Programming Problem/Linear Programming Problem.cpp`: algorithm implementation
- `Linear Programming Problem/Linear Programming Problem/Linear Programming Problem.vcxproj`: project configuration

## Algorithm Notes
- Input data is currently hardcoded in `main()`.
- The program is designed for demonstration and learning, not as a general-purpose LP parser.

## Known Limitations
- No file or interactive coefficient input beyond method selection
- No validation for arbitrary LP model formats
- Rounding/printing behavior is fixed in source
- Build automation is not included in this repository

## Release Package
This repository can publish a Windows executable package as a zip asset.

Current environment note:
- `msbuild` is unavailable in this environment, so Release rebuild was not possible here.
- Existing local artifact path detected: `Linear Programming Problem/x64/Debug/Linear Programming Problem.exe`.
- Prepared uploadable package: `ReleaseAssets/Linear-Programming-Problem-v1.0.0-win-x64-debug.zip`.

## License
MIT License. See `LICENSE`.