import { Link } from "react-router-dom";

export default function HomePage() {
  return (
    <div>
      <h1>Familias</h1>
      <p style={{ color: "var(--muted)", maxWidth: "60ch" }}>
        Compute likelihood ratios and posterior probabilities for kinship
        and parentage cases. Choose a workflow:
      </p>
      <div className="cards">
        <Link to="/one-parent" className="card">
          <h3>One-parent</h3>
          <p>
            Tested Parent vs. unrelated. Two persons (Tested Parent, Child).
            Per-locus and combined LR with live updates.
          </p>
        </Link>
        <Link to="/two-parent" className="card">
          <h3>Two-parent</h3>
          <p>
            Known Parent + Tested Parent vs. Known Parent only.
            Three persons (Known Parent, Tested Parent, Child).
          </p>
        </Link>
        <Link to="/arbitrary" className="card">
          <h3>Arbitrary relationships</h3>
          <p>
            Up to 8 persons with custom parent→child relations.
            Mark relations as &quot;test&quot; to define H<sub>1</sub> vs. H<sub>2</sub>.
          </p>
        </Link>
      </div>
    </div>
  );
}
