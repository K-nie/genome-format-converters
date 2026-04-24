# Contributing

Thanks for considering a contribution. This is a small, single-author
package, so the workflow is lightweight.

## Development setup

```bash
git clone https://github.com/K-nie/genome-format-converters.git
cd genome-format-converters
python -m venv .venv && source .venv/bin/activate
pip install -e ".[test]"
pytest -v
```

## Adding a new converter

1. Drop a `src/genome_format_converters/converters/<name>.py` module that
   exposes a `batch_convert(input_dir: str, output_dir: str, **kwargs)`
   entry point. Use helpers from `converters/_common.py`
   (`log_info`, `log_warn`, `prepare_output_dir`, `iter_input_files`,
   `strip_compound_suffix`, `load_two_col_map`) rather than rolling your
   own.
2. Register the module in `converters/__init__.py` and add a subparser +
   dispatch branch in `cli.py`. Use the `_add_io_args` / `_add_batch_flags`
   helpers so the new subcommand inherits the standard `--force`,
   `--pattern`, `--dry-run`, `--input-dir`/`--input` flag surface.
3. Add a test in `tests/test_converters.py` that exercises the happy path
   with a fixture from `tests/test_data/`. Where behaviour is subtle (file
   formats with deterministic output, EIGENSTRAT encodings, etc.), pin a
   golden file.
4. Update the `Command Reference` and `Scripts Overview` tables in
   `README.md` and add a line to `CHANGELOG.md` under `## [Unreleased]`.

## Style

- Format: keep modules short and single-purpose; helpers go in
  `converters/_common.py` or a purpose-specific `_*_common.py`.
- Logging: progress lines via `log_info`, warnings via `log_warn`,
  fatal errors via `die`. **Never** write progress to stdout — stdout is
  reserved for converter output where it makes sense.
- Filenames: use `strip_compound_suffix` for any format with a
  `.something.gz` compound extension.
- Errors: raise loudly on malformed input or exit with a non-zero code
  via `die`. Silent failures are harder to notice in a pipeline.
- Python: target ≥3.9; prefer `pathlib.Path` over `os.path`, prefer
  f-strings, use type hints where they help.

## Versioning

Semantic versioning. Breaking CLI changes require a minor bump; bug
fixes and non-breaking additions go out as patch releases. Update
`pyproject.toml` and `CHANGELOG.md` together.

## License

By contributing you agree to license your contribution under the MIT
License — the same terms the rest of the project ships under.
