import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import { api, type ComputeResponse, type Person, type Relation } from "../apiClient";
import {
  buildCase,
  DEFAULT_PARAMS,
  type LocusRow,
  type Mode,
  type ModeState,
  rowToLocusInput,
} from "./caseState";

interface InitArgs {
  mode: Mode;
  initialPersons: Person[];
  initialRelations?: Relation[];
}

export function useModeCase({ mode, initialPersons, initialRelations }: InitArgs) {
  const [state, setState] = useState<ModeState>({
    mode,
    persons: initialPersons,
    relations: initialRelations ?? [],
    loci: [],
    params: { ...DEFAULT_PARAMS },
  });
  const [databases, setDatabases] = useState<string[]>([]);
  const [databaseCache, setDatabaseCache] = useState<
    Record<string, Record<string, Record<string, number>>>
  >({});
  const [computing, setComputing] = useState(false);
  const [combined, setCombined] = useState<ComputeResponse | undefined>();
  const [topWarnings, setTopWarnings] = useState<string[]>([]);

  // ---- bootstrap: fetch DB list + default loci (rows start with empty freqs)
  useEffect(() => {
    let cancelled = false;
    (async () => {
      try {
        const [dbs, defLoci] = await Promise.all([
          api.databases(),
          api.defaultLoci(),
        ]);
        if (cancelled) return;
        setDatabases(dbs.databases);
        const personIds = initialPersons.map(p => p.id);
        const empty: LocusRow[] = defLoci.loci.map(name => ({
          name,
          genos: Object.fromEntries(personIds.map(id => [id, ["", ""]])),
          freqs: {},
        }));
        const dbName = dbs.default;
        if (dbName) {
          const got = await api.database(dbName);
          if (cancelled) return;
          setDatabaseCache(c => ({ ...c, [dbName]: got.loci }));
        }
        setState(s => ({
          ...s,
          loci: empty,
          params: { ...s.params, database: dbName },
        }));
      } catch (e) {
        // eslint-disable-next-line no-console
        console.warn("Bootstrap failed:", e);
      }
    })();
    return () => { cancelled = true; };
  // eslint-disable-next-line react-hooks/exhaustive-deps
  }, []);

  // ---- helpers ----
  const personIds = useMemo(() => state.persons.map(p => p.id), [state.persons]);

  const setLoci = useCallback((next: LocusRow[]) => {
    setState(s => ({ ...s, loci: next }));
  }, []);

  const setParams = useCallback((patch: Partial<ModeState["params"]>) => {
    setState(s => ({ ...s, params: { ...s.params, ...patch } }));
  }, []);

  const applyDatabase = useCallback(async () => {
    const name = state.params.database;
    if (!name) return;
    let body = databaseCache[name];
    if (!body) {
      const got = await api.database(name);
      body = got.loci;
      setDatabaseCache(c => ({ ...c, [name]: body }));
    }
    setState(s => {
      const next = s.loci.map(r => {
        const lookup = body[r.name] ?? findCaseInsensitive(body, r.name);
        if (!lookup) return r;
        // Reset freqs to DB values for currently-observed alleles only.
        const observed = observedAlleles(r);
        const freqs: Record<string, number> = {};
        for (const a of observed) {
          if (lookup[a] !== undefined) freqs[a] = lookup[a];
        }
        return { ...r, freqs };
      });
      return { ...s, loci: next };
    });
  }, [state.params.database, databaseCache]);

  // Auto-sync per-row freqs to currently-observed alleles:
  //   - drop frequencies for alleles no longer observed
  //   - fill missing observed alleles from active DB cache (when present)
  // Manual user entries for observed alleles are preserved.
  useEffect(() => {
    const dbName = state.params.database;
    const body = dbName ? databaseCache[dbName] : undefined;
    let changed = false;
    const next = state.loci.map(r => {
      const observed = observedAlleles(r);
      const observedSet = new Set(observed);
      const lookup = body
        ? (body[r.name] ?? findCaseInsensitive(body, r.name))
        : undefined;
      const newFreqs: Record<string, number> = {};
      // Keep existing entries that are still observed.
      for (const a of observed) {
        if (r.freqs[a] !== undefined) newFreqs[a] = r.freqs[a];
        else if (lookup && lookup[a] !== undefined) newFreqs[a] = lookup[a];
      }
      // Detect change vs current.
      const oldKeys = Object.keys(r.freqs);
      const same = oldKeys.length === Object.keys(newFreqs).length
        && oldKeys.every(k => observedSet.has(k) && r.freqs[k] === newFreqs[k]);
      if (same) return r;
      changed = true;
      return { ...r, freqs: newFreqs };
    });
    if (changed) setState(s => ({ ...s, loci: next }));
  // Re-run whenever genotypes change or the active DB changes.
  // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [JSON.stringify(state.loci.map(r => [r.name, r.genos])),
      state.params.database, databaseCache]);

  // ---- combined (full-case) compute ----
  const compute = useCallback(async () => {
    setComputing(true);
    setTopWarnings([]);
    try {
      const c = buildCase(state);
      let res: ComputeResponse;
      if (mode === "one-parent") res = await api.computeOneParent(c);
      else if (mode === "two-parent") res = await api.computeTwoParent(c);
      else res = await api.computeArbitrary({
        ...c,
        persons: state.persons,
        relations: state.relations,
      });
      setCombined(res);
      setTopWarnings(res.warnings ?? []);
    } catch (e: any) {
      setCombined(undefined);
      setTopWarnings([e?.message ?? String(e)]);
    } finally {
      setComputing(false);
    }
  }, [mode, state]);

  // ---- live per-locus LR (debounced per row) ----
  const debouncers = useRef<Record<string, number>>({});
  const aborters = useRef<Record<string, AbortController>>({});

  const lociRef = useRef(state.loci);
  useEffect(() => { lociRef.current = state.loci; }, [state.loci]);

  useEffect(() => {
    // For each row, schedule a debounced single-locus call.
    state.loci.forEach((row, idx) => {
      const li = rowToLocusInput(row);
      if (!li) return;
      const key = `${row.name}#${idx}`;
      if (debouncers.current[key]) window.clearTimeout(debouncers.current[key]);
      debouncers.current[key] = window.setTimeout(async () => {
        aborters.current[key]?.abort();
        const ctrl = new AbortController();
        aborters.current[key] = ctrl;
        try {
          const r = await api.singleLocus({
            mode,
            locus: row.name,
            frequencies: row.freqs,
            genotypes: li.genotypes,
            mutation_model: state.params.mutation_model,
            mutation_rate: state.params.mutation_rate,
            kinship: state.params.kinship,
            persons: mode === "arbitrary" ? state.persons : undefined,
            relations: mode === "arbitrary" ? state.relations : undefined,
          }, ctrl.signal);
          // Update only the row at the same name (avoid stomping on edits).
          setState(s => {
            const next = s.loci.map(rr =>
              rr.name === row.name
                ? { ...rr, liveLR: r.LR ?? null, liveReason: r.reason ?? null }
                : rr,
            );
            return { ...s, loci: next };
          });
        } catch (e: any) {
          if (e?.name === "AbortError") return;
          setState(s => {
            const next = s.loci.map(rr =>
              rr.name === row.name
                ? { ...rr, liveLR: null, liveReason: e?.message ?? "error" }
                : rr,
            );
            return { ...s, loci: next };
          });
        }
      }, 200);
    });
    return () => {
      // No-op: stale entries get aborted on the next scheduling pass.
    };
  // Note: depending on `state.loci` re-runs whenever any row's value changes.
  // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [JSON.stringify(state.loci.map(r =>
    [r.name, r.genos, r.freqs])), state.params.mutation_model,
      state.params.mutation_rate, state.params.kinship,
      JSON.stringify(state.persons), JSON.stringify(state.relations), mode]);

  return {
    state, setState, personIds,
    databases, computing, combined, topWarnings,
    setLoci, setParams, applyDatabase, compute,
  };
}

function findCaseInsensitive<T>(
  obj: Record<string, T>, key: string,
): T | undefined {
  const lower = key.toLowerCase();
  for (const k of Object.keys(obj)) {
    if (k.toLowerCase() === lower) return obj[k];
  }
  return undefined;
}

function observedAlleles(r: LocusRow): string[] {
  const out: string[] = [];
  for (const ab of Object.values(r.genos)) {
    for (const a of ab) {
      if (a && !out.includes(a)) out.push(a);
    }
  }
  return out;
}
