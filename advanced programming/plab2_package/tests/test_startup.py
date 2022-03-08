"""Tests for startup procedure."""

import os

from plab2.startup import DATA_DIR, LOG_DIR, LOG_FILE_PATH


class TestStartup:
    """Unit tests for the startup variables."""

    def test_cache_dirs(self):
        """Checks if cache directories are built."""
        assert os.path.isdir(DATA_DIR)
        assert os.path.isdir(LOG_DIR)

    def test_log_file(self):
        """Checks if the log file is made."""
        assert os.path.isfile(LOG_FILE_PATH)
