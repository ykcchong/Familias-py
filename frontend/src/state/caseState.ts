// Shared mode types and case-state helpers.

import type { CaseBase, ComputeResponse, Genotype, LocusInput, Person, Relation } from "../apiClient";

export type Mode = "one-parent" | "two-parent" | "arbitrary";

export interface LocusRow {
  name: string;
  // For each personId -> ["a","b"] or ["",""] when blank.
  genos: Record<string, [string, string]>;
  freqs: Record<string, number>;
  // Live per-locus LR, set asynchronously by debounced single-locus calls.
  liveLR?: number | null;
  liveReason?: string | null;
}

export interface ModeState {
  mode: Mode;
  persons: Person[];                // for arbitrary; otherwise canonical
  relations: Relation[];            // arbitrary only
  loci: LocusRow[];
  params: {
    mutation_model: string;
    mutation_rate: number;
    kinship: number;
    prior: number;
    database: string;
  };
  combined?: ComputeResponse;
}

/** Sentinel value for the database picker meaning "don't auto-apply
 * a database; preserve whatever frequencies are currently in each row."
 * Selecting this also makes the Apply DB button a no-op. */
export const MANUAL_DB = "__manual__";

export const DEFAULT_PARAMS = {
  mutation_model: "Equal",
  mutation_rate: 0,
  kinship: 0.0,
  prior: 0.5,
  database: "",
};

export const PERSONS_ONE_PARENT: Person[] = [
  { id: "AF", sex: "male" },
  { id: "CH", sex: "male" },
];

export const PERSONS_TWO_PARENT: Person[] = [
  { id: "AF", sex: "male" },
  { id: "MO", sex: "female" },
  { id: "CH", sex: "male" },
];

/**
 * Display labels for the fixed person ids in the one- and two-parent
 * workflows. The wire-protocol ids remain "AF"/"MO"/"CH" because the
 * backend pedigrees are hard-coded to those identifiers; only the UI
 * presentation differs (this app is geared toward generic parent/child
 * testing rather than strictly paternity).
 */
export const PERSON_LABELS_ONE_PARENT: Record<string, string> = {
  AF: "Tested Parent",
  CH: "Child",
};

export const PERSON_LABELS_TWO_PARENT: Record<string, string> = {
  AF: "Tested Parent",
  MO: "Known Parent",
  CH: "Child",
};

export function rowToLocusInput(row: LocusRow): LocusInput | null {
  const genos: Genotype[] = [];
  for (const [pid, ab] of Object.entries(row.genos)) {
    if (ab[0] && ab[1]) genos.push({ person: pid, a1: ab[0], a2: ab[1] });
  }
  if (genos.length === 0) return null;
  return { name: row.name, frequencies: row.freqs, genotypes: genos };
}

export function buildCase(state: ModeState): CaseBase {
  const loci = state.loci
    .map(rowToLocusInput)
    .filter((x): x is LocusInput => x !== null);
  return {
    mutation_model: state.params.mutation_model,
    mutation_rate: state.params.mutation_rate,
    kinship: state.params.kinship,
    prior: state.params.prior,
    loci,
  };
}

// --- TUI case JSON schema interop ---------------------------------------
// {
//   schema: "familias-py-tui-case/1",
//   mode: "One-parent mode" | "Two-parent mode" | "Arbitrary mode",
//   persons: [...],
//   mut_rate, kinship, database,
//   loci: { name: { genotypes: {pid: "a" | "a,b"}, freqs: {a: "0.12"} } },
//   amelo: { pid: "XX"|"XY" },
//   relations?: [{parent, child, flag}]   (arbitrary only)
// }

export interface TuiCaseFile {
  schema: "familias-py-tui-case/1";
  mode: string;
  persons: string[];
  mut_rate: number;
  kinship: number;
  database: string;
  loci: Record<string, {
    genotypes: Record<string, string>;
    freqs: Record<string, string | number>;
  }>;
  amelo?: Record<string, string>;
  relations?: { parent: string; child: string; flag: string }[];
}

export function exportCase(state: ModeState): TuiCaseFile {
  const modeLabel = state.mode === "one-parent" ? "One-parent mode"
    : state.mode === "two-parent" ? "Two-parent mode"
    : "Arbitrary mode";
  const loci: TuiCaseFile["loci"] = {};
  for (const r of state.loci) {
    const g: Record<string, string> = {};
    for (const [pid, ab] of Object.entries(r.genos)) {
      if (!ab[0] && !ab[1]) continue;
      g[pid] = ab[0] === ab[1] ? ab[0] : `${ab[0]},${ab[1]}`;
    }
    const f: Record<string, string> = {};
    for (const [a, v] of Object.entries(r.freqs)) f[a] = String(v);
    loci[r.name] = { genotypes: g, freqs: f };
  }
  const out: TuiCaseFile = {
    schema: "familias-py-tui-case/1",
    mode: modeLabel,
    persons: state.persons.map(p => p.id),
    mut_rate: state.params.mutation_rate,
    kinship: state.params.kinship,
    database: state.params.database,
    loci,
  };
  if (state.mode === "arbitrary") {
    out.relations = state.relations.map(r => ({ ...r }));
  }
  return out;
}

export function importCase(file: TuiCaseFile, current: ModeState): ModeState {
  const persons = current.persons;  // mode-fixed for one/two; arbitrary handled below
  const loci: LocusRow[] = Object.entries(file.loci).map(([name, e]) => {
    const genos: Record<string, [string, string]> = {};
    for (const pid of (file.persons ?? persons.map(p => p.id))) {
      const v = e.genotypes[pid] ?? "";
      const parts = v.split(",").map(s => s.trim()).filter(Boolean);
      if (parts.length === 0) genos[pid] = ["", ""];
      else if (parts.length === 1) genos[pid] = [parts[0], parts[0]];
      else genos[pid] = [parts[0], parts[1]];
    }
    const freqs: Record<string, number> = {};
    for (const [a, v] of Object.entries(e.freqs)) {
      const n = typeof v === "number" ? v : parseFloat(v);
      if (Number.isFinite(n) && n > 0) freqs[a] = n;
    }
    return { name, genos, freqs };
  });
  let nextPersons = persons;
  let nextRelations = current.relations;
  if (current.mode === "arbitrary" && Array.isArray(file.persons)) {
    // Best-effort: assume all persons male unless amelo says otherwise.
    nextPersons = file.persons.map(id => {
      const a = file.amelo?.[id];
      const sex: "male" | "female" = a === "XX" ? "female" : "male";
      return { id, sex };
    });
    if (file.relations) {
      nextRelations = file.relations.map(r => ({
        parent: r.parent,
        child: r.child,
        flag: r.flag === "test" ? "test" : "fixed",
      }));
    }
  }
  return {
    ...current,
    persons: nextPersons,
    relations: nextRelations,
    loci,
    params: {
      ...current.params,
      mutation_rate: file.mut_rate ?? current.params.mutation_rate,
      kinship: file.kinship ?? current.params.kinship,
      // After loading a case file, force Manual mode so the auto-sync
      // does not overwrite the file's verbatim frequencies. The user can
      // still pick a DB and click "Apply DB" to override.
      database: MANUAL_DB,
    },
  };
}
