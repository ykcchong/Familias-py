import type { ModeState } from "../state/caseState";
import { exportCase, importCase, type TuiCaseFile } from "../state/caseState";

interface Props {
  state: ModeState;
  onLoad: (next: ModeState) => void;
}

export default function CaseIO({ state, onLoad }: Props) {
  const onSave = () => {
    const data = exportCase(state);
    const blob = new Blob([JSON.stringify(data, null, 2)],
      { type: "application/json" });
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = `${state.mode}-case.json`;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    URL.revokeObjectURL(url);
  };

  const onLoadFile = (e: React.ChangeEvent<HTMLInputElement>) => {
    const f = e.target.files?.[0];
    e.target.value = "";
    if (!f) return;
    const reader = new FileReader();
    reader.onload = () => {
      try {
        const obj = JSON.parse(String(reader.result)) as TuiCaseFile;
        if (obj.schema !== "familias-py-tui-case/1") {
          alert(`Unrecognised schema: ${(obj as any).schema ?? "(none)"}`);
          return;
        }
        onLoad(importCase(obj, state));
      } catch (err: any) {
        alert("Failed to load case: " + (err?.message ?? err));
      }
    };
    reader.readAsText(f);
  };

  return (
    <div className="actions">
      <button className="secondary" onClick={onSave}>Save case</button>
      <label className="secondary" style={{
        display: "inline-block", padding: "0.4rem 0.85rem",
        border: "1px solid var(--border)", borderRadius: 4, cursor: "pointer",
      }}>
        Load case
        <input type="file" accept="application/json,.json"
          style={{ display: "none" }} onChange={onLoadFile} />
      </label>
    </div>
  );
}
