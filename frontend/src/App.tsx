import { Link, Outlet } from "react-router-dom";

export default function App() {
  return (
    <div className="app">
      <header className="topbar">
        <Link to="/" className="brand">Familias</Link>
        <nav>
          <Link to="/one-parent">One-parent</Link>
          <Link to="/two-parent">Two-parent</Link>
          <Link to="/arbitrary">Arbitrary</Link>
        </nav>
      </header>
      <main className="content">
        <Outlet />
      </main>
    </div>
  );
}
