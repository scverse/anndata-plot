import pytest


@pytest.mark.skip(reason="This decorator should be removed when test passes.")
def test_example() -> None:
    assert 1 == 0  # This test is designed to fail.
