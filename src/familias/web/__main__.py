"""``familias-web`` console entry point.

Launches uvicorn serving the FastAPI app and (by default) opens the
browser to the home page.
"""
from __future__ import annotations
import argparse
import threading
import time
import webbrowser
from pathlib import Path


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(prog="familias-web",
                                     description="Run the Familias web UI.")
    parser.add_argument("--host", default="127.0.0.1",
                        help="Bind host (default 127.0.0.1; use 0.0.0.0 to "
                             "expose on the LAN — note: no auth).")
    parser.add_argument("--port", type=int, default=8765,
                        help="Bind port (default 8765).")
    parser.add_argument("--no-browser", action="store_true",
                        help="Don't open the browser automatically.")
    parser.add_argument("--reload", action="store_true",
                        help="Enable uvicorn auto-reload (development).")
    parser.add_argument("--ssl-certfile", default=None,
                        help="Path to a PEM-encoded TLS certificate. "
                             "Pass with --ssl-keyfile to enable HTTPS.")
    parser.add_argument("--ssl-keyfile", default=None,
                        help="Path to a PEM-encoded TLS private key. "
                             "Required together with --ssl-certfile.")
    parser.add_argument("--ssl-keyfile-password", default=None,
                        help="Optional passphrase for the TLS private key.")
    args = parser.parse_args(argv)

    try:
        import uvicorn
    except ImportError as e:                                    # pragma: no cover
        raise SystemExit(
            "The 'web' extra is required: pip install 'familias[web]'."
        ) from e

    # ---- TLS option validation ----
    ssl_kwargs: dict = {}
    if bool(args.ssl_certfile) ^ bool(args.ssl_keyfile):
        raise SystemExit(
            "--ssl-certfile and --ssl-keyfile must be provided together.")
    if args.ssl_certfile and args.ssl_keyfile:
        for label, p in (("--ssl-certfile", args.ssl_certfile),
                         ("--ssl-keyfile", args.ssl_keyfile)):
            if not Path(p).is_file():
                raise SystemExit(f"{label}: file not found: {p}")
        ssl_kwargs = {
            "ssl_certfile": args.ssl_certfile,
            "ssl_keyfile": args.ssl_keyfile,
        }
        if args.ssl_keyfile_password:
            ssl_kwargs["ssl_keyfile_password"] = args.ssl_keyfile_password

    scheme = "https" if ssl_kwargs else "http"
    display_host = args.host if args.host != "0.0.0.0" else "127.0.0.1"
    url = f"{scheme}://{display_host}:{args.port}/"

    if not args.no_browser:
        def _open_later() -> None:
            time.sleep(1.0)
            try:
                webbrowser.open(url)
            except Exception:
                pass
        threading.Thread(target=_open_later, daemon=True).start()

    print(f"Familias web UI on {url}")
    if args.reload:
        uvicorn.run(
            "familias.web.app:app",
            host=args.host,
            port=args.port,
            reload=True,
            **ssl_kwargs,
        )
    else:
        from .app import app
        uvicorn.run(app, host=args.host, port=args.port, **ssl_kwargs)


if __name__ == "__main__":                                      # pragma: no cover
    main()
