# Familias Web Frontend

React 18 + TypeScript + Vite + AG Grid Community.

## Dev workflow

In one terminal, start the FastAPI backend:

```bash
pip install -e ".[web,tui]"
familias-web --no-browser --port 8765
# or with auto-reload:
uvicorn familias.web.app:app --reload --port 8765
```

In another, start the Vite dev server (proxies `/api` to 8765):

```bash
cd frontend
npm install
npm run dev   # http://localhost:5173
```

## Production build

```bash
cd frontend
npm install
npm run build
```

This writes the bundle into [`../src/familias/web/dist/`](../src/familias/web/dist/),
which is shipped inside the Python wheel. After building, commit the
`dist/` directory so end users do not need Node installed.
