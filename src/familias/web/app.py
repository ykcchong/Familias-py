"""FastAPI app: REST API + static SPA serving."""
from __future__ import annotations
from importlib import resources
from pathlib import Path

from fastapi import FastAPI, HTTPException, Request
from fastapi.responses import FileResponse, HTMLResponse, JSONResponse
from fastapi.staticfiles import StaticFiles

from . import compute_api as ca
from .models import (
    ArbitraryCase,
    ComputeResponse,
    DatabaseFullResponse,
    DatabaseListResponse,
    DatabaseLociResponse,
    MendelianRequest,
    MendelianResponse,
    OneParentCase,
    SingleLocusRequest,
    SingleLocusResponse,
    TwoParentCase,
)
from ..tui import compute as tc
from ..tui import freq_loader


# ---------------------------------------------------------------------------
# App factory
# ---------------------------------------------------------------------------
def create_app() -> FastAPI:
    app = FastAPI(
        title="Familias-py Web",
        description="REST API + bundled UI for kinship/paternity LRs.",
        version="1.0.0",
    )

    # -------- API: meta --------
    @app.get("/api/health")
    def health() -> dict:
        return {"ok": True}

    # -------- API: frequencies --------
    @app.get("/api/frequencies/databases", response_model=DatabaseListResponse)
    def list_databases() -> DatabaseListResponse:
        dbs = list(freq_loader.builtin_databases().keys())
        # Prefer the bundled Hong Kong Chinese DB; fall back to the first
        # available entry.
        default = next(
            (d for d in dbs if d.lower().startswith("hong kong")),
            dbs[0] if dbs else "",
        )
        return DatabaseListResponse(databases=dbs, default=default)

    def _get_db(name: str) -> dict:
        dbs = freq_loader.builtin_databases()
        if name not in dbs:
            raise HTTPException(404, f"Unknown database {name!r}.")
        return dbs[name]

    @app.get("/api/frequencies/{db_name}", response_model=DatabaseFullResponse)
    def get_database(db_name: str) -> DatabaseFullResponse:
        db = _get_db(db_name)
        return DatabaseFullResponse(database=db_name, loci=db)

    @app.get("/api/frequencies/{db_name}/loci",
             response_model=DatabaseLociResponse)
    def get_database_loci(db_name: str) -> DatabaseLociResponse:
        db = _get_db(db_name)
        return DatabaseLociResponse(database=db_name, loci=sorted(db.keys()))

    # -------- API: defaults --------
    @app.get("/api/defaults/loci")
    def default_loci() -> dict:
        from ..tui.defaults import DEFAULT_LOCI
        return {"loci": list(DEFAULT_LOCI)}

    # -------- API: compute --------
    @app.post("/api/compute/single-locus", response_model=SingleLocusResponse)
    def single_locus(req: SingleLocusRequest) -> SingleLocusResponse:
        try:
            return ca.compute_single_locus(req)
        except Exception as e:
            return SingleLocusResponse(LR=None, reason=str(e))

    @app.post("/api/compute/one-parent", response_model=ComputeResponse)
    def one_parent(case: OneParentCase) -> ComputeResponse:
        try:
            return ca.compute_one_parent(case)
        except ValueError as e:
            raise HTTPException(400, str(e))

    @app.post("/api/compute/two-parent", response_model=ComputeResponse)
    def two_parent(case: TwoParentCase) -> ComputeResponse:
        try:
            return ca.compute_two_parent(case)
        except ValueError as e:
            raise HTTPException(400, str(e))

    @app.post("/api/compute/arbitrary", response_model=ComputeResponse)
    def arbitrary(case: ArbitraryCase) -> ComputeResponse:
        try:
            return ca.compute_arbitrary(case)
        except ValueError as e:
            raise HTTPException(400, str(e))

    # -------- API: mendelian check --------
    @app.post("/api/check/mendelian", response_model=MendelianResponse)
    def mendelian(req: MendelianRequest) -> MendelianResponse:
        if req.mode == "one-parent":
            mismatch = tc.is_mismatch_one_parent(req.parent, req.child)
            pe = None
            if req.child and req.frequencies:
                pe = tc.power_of_exclusion_one_parent(
                    req.child, req.frequencies)
            return MendelianResponse(mismatch=mismatch, power_of_exclusion=pe)
        else:
            mismatch = tc.is_mismatch_two_parents(
                req.known, req.alleged, req.child)
            return MendelianResponse(mismatch=mismatch)

    # -------- Static SPA --------
    _mount_spa(app)

    return app


# ---------------------------------------------------------------------------
# SPA serving
# ---------------------------------------------------------------------------
def _dist_dir() -> Path:
    """Locate the bundled frontend dist directory."""
    try:
        with resources.as_file(resources.files("familias.web") / "dist") as p:
            return Path(p)
    except (ModuleNotFoundError, FileNotFoundError):
        return Path(__file__).parent / "dist"


_PLACEHOLDER_HTML = """<!doctype html>
<html><head><meta charset="utf-8"><title>Familias web — frontend not built</title>
<style>body{font-family:system-ui,sans-serif;max-width:48rem;margin:3rem auto;
padding:0 1rem;line-height:1.5}code{background:#f3f3f3;padding:.1em .3em;
border-radius:3px}</style></head>
<body><h1>Familias web</h1>
<p>The REST API is running, but the React frontend bundle was not found at
<code>src/familias/web/dist/</code>.</p>
<p>To build it:</p>
<pre><code>cd frontend
npm install
npm run build</code></pre>
<p>Or run the dev server with Vite proxy:</p>
<pre><code>uvicorn familias.web.app:app --reload &amp;
cd frontend &amp;&amp; npm run dev</code></pre>
<p>API endpoints are listed at <a href="/docs">/docs</a>.</p>
</body></html>"""


def _mount_spa(app: FastAPI) -> None:
    dist = _dist_dir()
    index_html = dist / "index.html"
    assets_dir = dist / "assets"

    if assets_dir.is_dir():
        app.mount("/assets", StaticFiles(directory=str(assets_dir)),
                  name="assets")

    # Common single-file static assets at root (favicon, etc.)
    @app.get("/favicon.ico", include_in_schema=False, response_model=None)
    def favicon() -> FileResponse | JSONResponse:
        f = dist / "favicon.ico"
        if f.is_file():
            return FileResponse(str(f))
        return JSONResponse(status_code=404, content={"detail": "no favicon"})

    @app.get("/", include_in_schema=False, response_model=None)
    def root() -> HTMLResponse | FileResponse:
        if index_html.is_file():
            return FileResponse(str(index_html))
        return HTMLResponse(_PLACEHOLDER_HTML)

    # SPA fallback for client-side routes (anything not /api, /assets, /docs).
    @app.get("/{full_path:path}", include_in_schema=False, response_model=None)
    def spa_fallback(full_path: str, request: Request
                     ) -> HTMLResponse | FileResponse | JSONResponse:
        if (full_path.startswith("api/") or full_path.startswith("assets/")
                or full_path in ("openapi.json", "docs", "redoc")):
            return JSONResponse(status_code=404,
                                content={"detail": "Not Found"})
        if index_html.is_file():
            return FileResponse(str(index_html))
        return HTMLResponse(_PLACEHOLDER_HTML)


# Module-level app instance (used by ``uvicorn familias.web.app:app``).
app = create_app()
