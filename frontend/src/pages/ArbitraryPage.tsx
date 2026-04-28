import { useModeCase } from "../state/useModeCase";
import LocusGrid from "../components/LocusGrid";
import ParamsBar from "../components/ParamsBar";
import RelationsGrid from "../components/RelationsGrid";
import ResultsFooter from "../components/ResultsFooter";
import CaseIO from "../components/CaseIO";

const INITIAL_PERSONS = [
  { id: "P1", sex: "male" as const },
  { id: "P2", sex: "male" as const },
];

export default function ArbitraryPage() {
  const m = useModeCase({
    mode: "arbitrary",
    initialPersons: INITIAL_PERSONS,
    initialRelations: [{ parent: "P1", child: "P2", flag: "test" }],
  });

  return (
    <div>
      <h1>Arbitrary relationships</h1>
      <p style={{ color: "var(--muted)", marginTop: 0 }}>
        Define 2–8 individuals and parent→child relations. Mark at least one
        as <strong>test</strong>: the LR contrasts H<sub>1</sub> (with the
        tested relations) vs. H<sub>2</sub> (without them).
      </p>
      <RelationsGrid
        persons={m.state.persons}
        relations={m.state.relations}
        onPersonsChange={(persons) => {
          const ids = persons.map(p => p.id);
          m.setState(s => ({
            ...s,
            persons,
            // Keep loci row genotype dicts in sync with current person ids.
            loci: s.loci.map(r => {
              const genos = { ...r.genos };
              for (const id of ids) if (!genos[id]) genos[id] = ["", ""];
              for (const k of Object.keys(genos)) {
                if (!ids.includes(k)) delete genos[k];
              }
              return { ...r, genos };
            }),
          }));
        }}
        onRelationsChange={(relations) =>
          m.setState(s => ({ ...s, relations }))}
      />

      <h2>DNA</h2>
      <ParamsBar
        state={m.state}
        databases={m.databases}
        onChange={m.setParams}
        onApplyDatabase={m.applyDatabase}
        onCompute={m.compute}
        computing={m.computing}
      />
      <LocusGrid
        rows={m.state.loci}
        personIds={m.personIds}
        onChange={m.setLoci}
      />
      <CaseIO state={m.state} onLoad={m.setState} />
      <ResultsFooter result={m.combined} warnings={m.topWarnings} />
    </div>
  );
}
