import { useMemo, useRef, useEffect } from "react";
import { AgGridReact } from "ag-grid-react";
import type { ColDef, ICellRendererParams, ValueGetterParams, ValueSetterParams } from "ag-grid-community";
import type { LocusRow } from "../state/caseState";

interface Props {
  rows: LocusRow[];
  personIds: string[];
  /** Optional display labels keyed by person id. Falls back to id. */
  personLabels?: Record<string, string>;
  onChange: (next: LocusRow[]) => void;
  height?: number;
}

function genoStr(ab: [string, string]): string {
  if (!ab[0] && !ab[1]) return "";
  if (ab[0] === ab[1]) return ab[0];
  return `${ab[0]}/${ab[1]}`;
}

function parseGeno(s: string): [string, string] {
  const parts = s.split(/[\/,\s]+/).map(x => x.trim()).filter(Boolean);
  if (parts.length === 0) return ["", ""];
  if (parts.length === 1) return [parts[0], parts[0]];
  return [parts[0], parts[1]];
}

function freqsStr(f: Record<string, number>): string {
  return Object.entries(f)
    .map(([a, v]) => `${a}:${v}`)
    .join(", ");
}

function parseFreqs(s: string): Record<string, number> {
  const out: Record<string, number> = {};
  for (const tok of s.split(/[,;]+/)) {
    const t = tok.trim();
    if (!t) continue;
    const m = t.match(/^([^:=\s]+)\s*[:=]\s*([\d.eE+-]+)$/);
    if (!m) continue;
    const v = parseFloat(m[2]);
    if (Number.isFinite(v) && v > 0) out[m[1]] = v;
  }
  return out;
}

function fmtLR(v: number | null | undefined, reason?: string | null): string {
  if (v === null || v === undefined) return reason ? "—" : "";
  if (!Number.isFinite(v)) return String(v);
  // Plain decimal display (no scientific notation): up to 2 fractional digits.
  return v.toFixed(2);
}

export default function LocusGrid({ rows, personIds, personLabels, onChange, height = 480 }: Props) {
  const gridRef = useRef<AgGridReact<LocusRow>>(null);

  const colDefs = useMemo<ColDef<LocusRow>[]>(() => {
    const cols: ColDef<LocusRow>[] = [
      {
        headerName: "Locus",
        field: "name",
        width: 130,
        editable: true,
        pinned: "left",
      },
    ];
    for (const pid of personIds) {
      cols.push({
        headerName: personLabels?.[pid] ?? pid,
        colId: `p_${pid}`,
        editable: true,
        width: 155,
        valueGetter: (p: ValueGetterParams<LocusRow>) =>
          p.data ? genoStr(p.data.genos[pid] ?? ["", ""]) : "",
        valueSetter: (p: ValueSetterParams<LocusRow>) => {
          if (!p.data) return false;
          p.data.genos = { ...p.data.genos, [pid]: parseGeno(String(p.newValue ?? "")) };
          return true;
        },
      });
    }
    cols.push({
      headerName: "Frequencies (allele:freq, …)",
      colId: "freqs",
      editable: true,
      flex: 1,
      minWidth: 220,
      valueGetter: (p: ValueGetterParams<LocusRow>) =>
        p.data ? freqsStr(p.data.freqs) : "",
      valueSetter: (p: ValueSetterParams<LocusRow>) => {
        if (!p.data) return false;
        p.data.freqs = parseFreqs(String(p.newValue ?? ""));
        return true;
      },
    });
    cols.push({
      headerName: "LR",
      colId: "lr",
      width: 110,
      editable: false,
      valueGetter: (p: ValueGetterParams<LocusRow>) =>
        p.data ? fmtLR(p.data.liveLR, p.data.liveReason) : "",
      cellRenderer: (p: ICellRendererParams<LocusRow>) => {
        const r = p.data;
        if (!r) return "";
        const v = r.liveLR;
        const txt = fmtLR(v, r.liveReason);
        const cls = v !== null && v !== undefined && v < 1 ? "mismatch" : "";
        return <span className={cls} title={r.liveReason ?? ""}>{txt || "—"}</span>;
      },
    });
    return cols;
  }, [personIds, personLabels]);

  // When external rows change, refresh grid display.
  useEffect(() => {
    gridRef.current?.api?.refreshCells({ force: true });
  }, [rows]);

  return (
    <div className="ag-theme-quartz grid-wrap" style={{ height, width: "100%" }}>
      <AgGridReact<LocusRow>
        ref={gridRef}
        rowData={rows}
        columnDefs={colDefs}
        defaultColDef={{ resizable: true, sortable: false }}
        getRowId={p => p.data.name}
        singleClickEdit={false}
        stopEditingWhenCellsLoseFocus
        onCellValueChanged={() => {
          // The grid mutated row data in place; surface a fresh array.
          onChange([...rows]);
        }}
        animateRows={false}
      />
    </div>
  );
}
