"""The CLI must not import converter dependencies until a subcommand runs."""
import subprocess
import sys

HEAVY = ("Bio", "BCBio", "pandas", "numpy", "cyvcf2", "pysam")


def _loaded(code: str):
    out = subprocess.run(
        [sys.executable, "-c",
         code + f"; import sys; print(','.join(m for m in {HEAVY!r} if m in sys.modules))"],
        capture_output=True, text=True, check=True).stdout.strip()
    return set(filter(None, out.split(",")))


def test_cli_import_loads_no_heavy_libraries():
    assert _loaded("import genome_format_converters.cli") == set()


def test_annotation_converter_stays_light():
    assert _loaded("import genome_format_converters.cli; "
                   "from genome_format_converters.cli import hmmer_tblout_to_tsv as m; "
                   "m.batch_convert") == set()


def test_lazy_module_resolves_attributes():
    from genome_format_converters import cli
    assert callable(cli.gff3_to_gtf.batch_convert)
