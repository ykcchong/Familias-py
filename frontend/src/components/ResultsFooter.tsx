import type { ComputeResponse } from "../apiClient";

interface Props {
  result?: ComputeResponse;
  warnings: string[];
}

function fmt(n: number, digits = 4): string {
  if (!Number.isFinite(n)) return String(n);
  if (n === 0) return "0";
  if (Math.abs(n) >= 1e4 || Math.abs(n) < 1e-3) return n.toExponential(3);
  return n.toFixed(digits);
}

// LR is shown as a plain decimal (no scientific notation), with up to
// 4 fractional digits.
function fmtLR(n: number): string {
  if (!Number.isFinite(n)) return String(n);
  return n.toFixed(4);
}

// Posterior probability as a percentage, e.g. "99.9876%". Up to 4
// fractional digits.
function fmtPercent(p: number): string {
  if (!Number.isFinite(p)) return String(p);
  return `${(p * 100).toFixed(4)}%`;
}

export default function ResultsFooter({ result, warnings }: Props) {
  if (!result && warnings.length === 0) return null;
  return (
    <div>
      {result && (
        <div className="results">
          <div className="metric">
            <span className="lbl">LR (H₁ / H₂)</span>
            <span className="val">{fmtLR(result.LR)}</span>
          </div>
          <div className="metric">
            <span className="lbl">log₁₀ LR</span>
            <span className="val">{fmt(result.log10_LR, 4)}</span>
          </div>
          <div className="metric">
            <span className="lbl">Posterior P(H₁)</span>
            <span className="val">{fmtPercent(result.posterior_h1)}</span>
          </div>
          <div className="metric">
            <span className="lbl">L(H₁)</span>
            <span className="val">{fmt(result.likelihoods[0])}</span>
          </div>
          <div className="metric">
            <span className="lbl">L(H₂)</span>
            <span className="val">{fmt(result.likelihoods[1])}</span>
          </div>
        </div>
      )}
      {warnings.length > 0 && (
        <div className="warnings">
          <strong>Warnings:</strong>
          <ul>{warnings.map((w, i) => <li key={i}>{w}</li>)}</ul>
        </div>
      )}
    </div>
  );
}
