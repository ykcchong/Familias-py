import { useMemo } from "react";
import { AgGridReact } from "ag-grid-react";
import type { ColDef } from "ag-grid-community";
import type { Person, Relation } from "../apiClient";

interface Props {
  persons: Person[];
  relations: Relation[];
  onPersonsChange: (next: Person[]) => void;
  onRelationsChange: (next: Relation[]) => void;
}

export default function RelationsGrid({
  persons, relations, onPersonsChange, onRelationsChange,
}: Props) {
  const ids = useMemo(() => persons.map(p => p.id), [persons]);

  const personCols = useMemo<ColDef<Person>[]>(() => [
    { headerName: "ID", field: "id", editable: true, width: 120 },
    {
      headerName: "Sex", field: "sex", editable: true, width: 110,
      cellEditor: "agSelectCellEditor",
      cellEditorParams: { values: ["male", "female"] },
    },
  ], []);

  const relationCols = useMemo<ColDef<Relation>[]>(() => [
    {
      headerName: "Parent", field: "parent", editable: true, width: 130,
      cellEditor: "agSelectCellEditor",
      cellEditorParams: { values: ids },
    },
    {
      headerName: "Child", field: "child", editable: true, width: 130,
      cellEditor: "agSelectCellEditor",
      cellEditorParams: { values: ids },
    },
    {
      headerName: "Test (defines H₁)", field: "flag", editable: true, width: 160,
      cellEditor: "agSelectCellEditor",
      cellEditorParams: { values: ["fixed", "test"] },
    },
  ], [ids]);

  const numPersons = persons.length;
  const setNumPersons = (n: number) => {
    n = Math.max(2, Math.min(8, n | 0));
    if (n === numPersons) return;
    if (n > numPersons) {
      const extra: Person[] = [];
      for (let i = numPersons; i < n; i++) {
        extra.push({ id: `P${i + 1}`, sex: "male" });
      }
      onPersonsChange([...persons, ...extra]);
    } else {
      const keep = persons.slice(0, n);
      const keepIds = new Set(keep.map(p => p.id));
      onPersonsChange(keep);
      onRelationsChange(relations.filter(
        r => keepIds.has(r.parent) && keepIds.has(r.child)));
    }
  };

  const addRelation = () => {
    onRelationsChange([
      ...relations,
      { parent: ids[0] ?? "", child: ids[1] ?? "", flag: "test" },
    ]);
  };
  const removeLast = () => {
    onRelationsChange(relations.slice(0, -1));
  };

  return (
    <div>
      <h2>Persons</h2>
      <label className="params" style={{ marginBottom: "0.5rem" }}>
        <span style={{ fontSize: "0.8rem", color: "var(--muted)" }}>
          Number of individuals (2–8)
        </span>
        <input
          type="number" min={2} max={8} value={numPersons}
          onChange={e => setNumPersons(parseInt(e.target.value, 10) || 2)}
          style={{ width: "5rem" }}
        />
      </label>
      <div className="ag-theme-quartz grid-wrap" style={{ height: 200 }}>
        <AgGridReact<Person>
          rowData={persons}
          columnDefs={personCols}
          getRowId={p => p.data.id}
          stopEditingWhenCellsLoseFocus
          onCellValueChanged={() => onPersonsChange([...persons])}
          defaultColDef={{ resizable: true }}
        />
      </div>

      <h2>Relations</h2>
      <div className="actions" style={{ marginBottom: "0.5rem" }}>
        <button onClick={addRelation}>+ Add relation</button>
        <button className="secondary" onClick={removeLast}
          disabled={relations.length === 0}>− Remove last</button>
      </div>
      <div className="ag-theme-quartz grid-wrap" style={{ height: 220 }}>
        <AgGridReact<Relation>
          rowData={relations}
          columnDefs={relationCols}
          stopEditingWhenCellsLoseFocus
          onCellValueChanged={() => onRelationsChange([...relations])}
          defaultColDef={{ resizable: true }}
        />
      </div>
    </div>
  );
}
