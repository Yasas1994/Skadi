import pytest
from click.testing import CliRunner

from skadi.cli import cli, _run_command


class TestRunCommand:
    def test_run_command_success(self, monkeypatch):
        calls = []

        def fake_run(cmd, **kwargs):
            calls.append(cmd)
            return 0

        monkeypatch.setattr("skadi.cli.subprocess.run", fake_run)
        _run_command(["echo", "hello"])
        assert calls == [["echo", "hello"]]

    def test_run_command_failure(self, monkeypatch):
        import subprocess

        def fake_run(cmd, **kwargs):
            raise subprocess.CalledProcessError(1, cmd)

        monkeypatch.setattr("skadi.cli.subprocess.run", fake_run)
        with pytest.raises(SystemExit) as exc:
            _run_command(["false"])
        assert exc.value.code == 1


class TestCLI:
    def test_cli_help(self):
        runner = CliRunner()
        result = runner.invoke(cli, ["--help"])
        assert result.exit_code == 0
        assert "SKADI" in result.output

    def test_utils_help(self):
        runner = CliRunner()
        result = runner.invoke(cli, ["utils", "--help"])
        assert result.exit_code == 0
        assert "ani" in result.output
        assert "aai" in result.output
