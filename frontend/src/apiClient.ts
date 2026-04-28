// Minimal fetch wrapper for the Familias REST API.

export interface Genotype { person: string; a1: string; a2: string }
export interface LocusInput {
  name: string;
  frequencies: Record<string, number>;
  genotypes: Genotype[];
}
export interface Person { id: string; sex: "male" | "female" }
export interface Relation { parent: string; child: string; flag: "fixed" | "test" }

export interface CaseBase {
  mutation_model: string;
  mutation_rate: number;
  kinship: number;
  prior: number;
  loci: LocusInput[];
}

export interface ArbitraryCase extends CaseBase {
  persons: Person[];
  relations: Relation[];
}

export interface PerLocusResult { locus: string; L1: number; L2: number; LR: number }
export interface ComputeResponse {
  LR: number;
  log10_LR: number;
  likelihoods: [number, number];
  posterior: [number, number];
  posterior_h1: number;
  per_locus: PerLocusResult[];
  warnings: string[];
}

export interface SingleLocusResponse {
  LR: number | null;
  L1?: number | null;
  L2?: number | null;
  reason?: string | null;
}

async function jsonFetch<T>(url: string, init?: RequestInit): Promise<T> {
  const r = await fetch(url, {
    ...init,
    headers: { "Content-Type": "application/json", ...(init?.headers ?? {}) },
  });
  if (!r.ok) {
    let msg = `${r.status} ${r.statusText}`;
    try { const j = await r.json(); msg = j.detail ?? msg; } catch { /* ignore */ }
    throw new Error(msg);
  }
  return r.json() as Promise<T>;
}

export const api = {
  health: () => jsonFetch<{ ok: boolean }>("/api/health"),
  databases: () => jsonFetch<{ databases: string[]; default: string }>(
    "/api/frequencies/databases"),
  database: (name: string) => jsonFetch<{
    database: string; loci: Record<string, Record<string, number>>;
  }>(`/api/frequencies/${encodeURIComponent(name)}`),
  defaultLoci: () => jsonFetch<{ loci: string[] }>("/api/defaults/loci"),
  computeOneParent: (c: CaseBase) => jsonFetch<ComputeResponse>(
    "/api/compute/one-parent", { method: "POST", body: JSON.stringify(c) }),
  computeTwoParent: (c: CaseBase) => jsonFetch<ComputeResponse>(
    "/api/compute/two-parent", { method: "POST", body: JSON.stringify(c) }),
  computeArbitrary: (c: ArbitraryCase) => jsonFetch<ComputeResponse>(
    "/api/compute/arbitrary", { method: "POST", body: JSON.stringify(c) }),
  singleLocus: (
    body: {
      mode: "one-parent" | "two-parent" | "arbitrary";
      locus: string;
      frequencies: Record<string, number>;
      genotypes: Genotype[];
      mutation_model: string;
      mutation_rate: number;
      kinship: number;
      persons?: Person[];
      relations?: Relation[];
    },
    signal?: AbortSignal,
  ) => jsonFetch<SingleLocusResponse>("/api/compute/single-locus",
    { method: "POST", body: JSON.stringify(body), signal }),
};

// Promise-aware debouncer; cancels stale in-flight requests via AbortController.
export function makeDebouncer<T>(delayMs = 150) {
  let timer: number | null = null;
  let controller: AbortController | null = null;
  return (fn: (signal: AbortSignal) => Promise<T>): Promise<T | null> => {
    if (timer !== null) window.clearTimeout(timer);
    if (controller) controller.abort();
    controller = new AbortController();
    const c = controller;
    return new Promise((resolve) => {
      timer = window.setTimeout(async () => {
        try { resolve(await fn(c.signal)); }
        catch (e: any) { if (e?.name === "AbortError") resolve(null); else resolve(null); }
      }, delayMs);
    });
  };
}
