"""
Unit and CLI tests for ini_template_generator.py
"""

import sys
import tempfile
import shutil
from pathlib import Path
import subprocess
import pytest
from unittest.mock import patch, MagicMock
import scripts.ini_template_generator as itg


SCRIPT = Path(itg.__file__).resolve()

# ---------------------------------------------------------------------
# Fixtures and helpers
# ---------------------------------------------------------------------

@pytest.fixture
def temp_dir():
    """Temporary working directory for test isolation."""
    d = Path(tempfile.mkdtemp())
    yield d
    shutil.rmtree(d)


@pytest.fixture
def mock_logger():
    """Provide a dummy logger that won't write to disk."""
    mock_log = MagicMock()
    return mock_log


def run_script(*args, cwd=None, input_data=None, timeout=10):
    """Run the ini_template_generator as a subprocess."""
    cmd = [sys.executable, str(SCRIPT), *args]
    result = subprocess.run(
        cmd,
        cwd=cwd,
        input=input_data,
        capture_output=True,
        text=True,
        timeout=timeout,
    )
    return result

# ---------------------------------------------------------------------
# Main CLI Tests
# ---------------------------------------------------------------------

class TestIniTemplateGenerator:
    """Tests for ini_template_generator.py main CLI behavior."""

    def test_creates_file_and_log(self, temp_dir):
        """Normal run: should create .ini and .log files."""
        ini_path = temp_dir / "run01.ini"
        result = run_script(str(ini_path), cwd=temp_dir)

        assert result.returncode == 0, result.stderr
        assert ini_path.exists(), "INI file not created"

        log_dir = temp_dir / "logs"
        log_files = list(log_dir.glob("*.log"))
        assert len(log_files) > 0, "No log files created"

        content = ini_path.read_text(encoding="utf-8")
        assert "[general]" in content, "Missing INI content"
        assert "input_dir" in content, "Template not written properly"

    def test_force_overwrite_skips_prompt(self, temp_dir):
        """--force should overwrite without prompting."""
        ini_path = temp_dir / "overwrite.ini"
        run_script(str(ini_path), cwd=temp_dir)
        assert ini_path.exists()

        result = run_script(str(ini_path), "--force", cwd=temp_dir)
        assert result.returncode == 0
        assert ini_path.exists()

    def test_prompt_abort(self, temp_dir):
        """If file exists and user types 'n', should abort without overwrite."""
        ini_path = temp_dir / "abort.ini"
        run_script(str(ini_path), cwd=temp_dir)

        result = run_script(str(ini_path), input_data="n\n", cwd=temp_dir)

        assert "Aborted" in result.stdout or "Aborted" in result.stderr
        assert result.returncode == 0
        assert ini_path.exists()

    def test_no_log_file_option(self, temp_dir):
        """--no-log-file should skip writing to logs/"""
        ini_path = temp_dir / "nolog.ini"
        result = run_script(str(ini_path), "--no-log-file", cwd=temp_dir)
        assert result.returncode == 0
        assert ini_path.exists()
        assert not (temp_dir / "logs").exists(), "logs/ should not exist"

    def test_help_flag(self, temp_dir):
        """-h or --help should display usage and exit 0."""
        for flag in ("-h", "--help"):
            result = run_script(flag, cwd=temp_dir)
            assert result.returncode == 0
            assert "Usage:" in result.stdout

    def test_missing_filename(self, temp_dir):
        """Missing <output.ini> argument should print usage and return 1."""
        result = run_script(cwd=temp_dir)
        assert result.returncode == 1
        assert "Usage:" in result.stdout or "Usage:" in result.stderr


# ---------------------------------------------------------------------
# Edge / Unit-level Tests
# ---------------------------------------------------------------------

class TestEdgeCases:
    """Unit-level tests for edge cases and argument parsing."""

    def test_empty_argv(self):
        """Test handling of completely empty argv list."""
        result = itg.main([])
        assert result == 1

    def test_only_flags_no_filename(self):
        """Test that flags alone without filename fail gracefully."""
        with patch.object(itg, "setup_logging") as mock_setup:
            mock_setup.return_value = MagicMock()
            result = itg.main(["script.py", "--force"])
            assert result == 1

    def test_mixed_flags_and_multiple_files(self):
        """Test handling of multiple files with flags."""
        with patch.object(itg, "setup_logging") as mock_setup:
            mock_setup.return_value = MagicMock()
            result = itg.main(["script.py", "file1.ini", "--force", "file2.ini"])
            assert result == 1

    def test_file_in_nonexistent_deep_directory(self, tmp_path, mock_logger):
        """Test creating file in deeply nested non-existent directory."""
        output_file = tmp_path / "a" / "b" / "c" / "d" / "e" / "config.ini"
        itg.write_template_config(output_file, log=mock_logger, force=False)
        assert output_file.exists()
        assert output_file.parent.exists()