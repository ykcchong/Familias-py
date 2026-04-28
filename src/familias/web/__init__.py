"""Web interface for Familias-py.

Provides a FastAPI REST API and a bundled React + AG Grid frontend
served by the same uvicorn process. Mirrors the three TUI modes
(one-parent, two-parent, arbitrary).
"""
from .app import create_app

__all__ = ["create_app"]
