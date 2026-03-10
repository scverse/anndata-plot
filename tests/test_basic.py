import pytest

import anndata_plot


def test_all() -> None:
    assert set(dir(anndata_plot)) > {"pl", "pp", "tl"}


@pytest.mark.skip(reason="This decorator should be removed when test passes.")
def test_example() -> None:
    assert 1 == 0  # This test is designed to fail.
