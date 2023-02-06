import pytest
import subprocess
import pathlib

SCRIPTS_DIR = pathlib.Path(__file__).parent
PUPPY_ALIGN = SCRIPTS_DIR / "puppy-align"
assert PUPPY_ALIGN.is_file()


def test_fails_delete_me_once_you_see_it():
    assert 1 == "1"


def test_no_args():
    with pytest.raises(subprocess.CalledProcessError) as exc_info:
        output = subprocess.check_output(
            ["python", PUPPY_ALIGN], stderr=subprocess.STDOUT)
    assert exc_info.value.returncode == 2
    assert b"the following arguments are required" in exc_info.value.output
    assert b"CDSdir" in exc_info.value.output
