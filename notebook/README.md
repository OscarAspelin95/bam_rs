## Purpose
This interactive [marimo](https://marimo.io/) notebook is mean to be an exploratory complement to the output file from `bam_rs`.

## Features
- File browser for navigating to the `bam_rs` output file.
- Dataframe explorer.
- Various interactive plots.


## Dependencies
[uv](https://github.com/astral-sh/uv) package manager


## Usage
The easiest way to get started is by running
```bash
uv venv
uv run marimo edit
```
and running the notebook interactively.

NOTE - using `uv run marimo run notebook.py` will most likely generate an error until the output tsv file has been chosen.
