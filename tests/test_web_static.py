"""Tests for static SPA serving."""
from __future__ import annotations
from pathlib import Path

import pytest

pytest.importorskip("fastapi")
from fastapi.testclient import TestClient

from familias.web.app import create_app, _dist_dir


@pytest.fixture(scope="module")
def client() -> TestClient:
    return TestClient(create_app())


def test_root_returns_html(client: TestClient) -> None:
    r = client.get("/")
    assert r.status_code == 200
    assert "text/html" in r.headers.get("content-type", "")


def test_spa_fallback_for_client_route(client: TestClient) -> None:
    r = client.get("/one-parent")
    assert r.status_code == 200
    assert "text/html" in r.headers.get("content-type", "")


def test_unknown_api_path_404(client: TestClient) -> None:
    r = client.get("/api/this-does-not-exist")
    assert r.status_code == 404


def test_dist_dir_resolves() -> None:
    p = _dist_dir()
    assert isinstance(p, Path)
