import pytest
from click.testing import CliRunner

from skadi.cli import _validate_path, cli


class TestValidatePath:
    def test_safe_path(self):
        assert _validate_path("/home/user/data", "input") == "/home/user/data"

    def test_safe_relative_path(self):
        assert _validate_path("data/samples", "input") == "data/samples"

    def test_rejects_semicolon(self):
        with pytest.raises(Exception):
            _validate_path("/tmp; rm -rf /", "input")

    def test_rejects_pipe(self):
        with pytest.raises(Exception):
            _validate_path("/tmp | cat /etc/passwd", "input")

    def test_rejects_backtick(self):
        with pytest.raises(Exception):
            _validate_path("/tmp/`whoami`", "input")

    def test_rejects_dollar(self):
        with pytest.raises(Exception):
            _validate_path("/tmp/$HOME", "input")


class TestCLIValidation:
    def test_contigs_rejects_malicious_input(self):
        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "contigs",
                "-i", "/tmp; rm -rf /",
                "-o", "/tmp/out",
            ],
        )
        assert result.exit_code != 0
        assert "unsafe" in result.output.lower() or result.exit_code == 2
