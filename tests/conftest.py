import pytest


@pytest.fixture
def sample_header():
    return [
        "query", "target", "theader", "fident", "qlen", "tlen",
        "alnlen", "mismatch", "gapopen", "qstart", "qend",
        "tstart", "tend", "evalue", "bits", "taxid", "taxname", "taxlineage",
    ]
