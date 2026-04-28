import type { ModeState } from "../state/caseState";

interface Props {
  state: ModeState;
  databases: string[];
  onChange: (patch: Partial<ModeState["params"]>) => void;
  onApplyDatabase: () => void;
  onCompute: () => void;
  computing?: boolean;
}

export default function ParamsBar({
  state, databases, onChange, onApplyDatabase, onCompute, computing,
}: Props) {
  const p = state.params;
  return (
    <div className="params">
      <label>
        Database
        <select
          value={p.database}
          onChange={e => onChange({ database: e.target.value })}
        >
          {databases.map(d => <option key={d} value={d}>{d}</option>)}
        </select>
      </label>
      <button className="secondary" onClick={onApplyDatabase}>Apply DB</button>

      <label>
        Mutation model
        <select
          value={p.mutation_model}
          onChange={e => onChange({ mutation_model: e.target.value })}
        >
          <option value="Equal">Equal</option>
          <option value="Proportional">Proportional</option>
          <option value="Stepwise">Stepwise</option>
        </select>
      </label>
      <label>
        Mutation rate
        <input
          type="number" step="0.0001" min={0} max={1}
          value={p.mutation_rate}
          onChange={e => onChange({ mutation_rate: parseFloat(e.target.value) || 0 })}
        />
      </label>
      <label>
        Kinship θ
        <input
          type="number" step="0.001" min={0} max={1}
          value={p.kinship}
          onChange={e => onChange({ kinship: parseFloat(e.target.value) || 0 })}
        />
      </label>
      <label>
        Prior P(H₁)
        <input
          type="number" step="0.01" min={0} max={1}
          value={p.prior}
          onChange={e => onChange({ prior: parseFloat(e.target.value) || 0 })}
        />
      </label>

      <button onClick={onCompute} disabled={computing}>
        {computing ? "Computing…" : "Recompute"}
      </button>
    </div>
  );
}
