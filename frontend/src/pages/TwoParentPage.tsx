import { useModeCase } from "../state/useModeCase";
import { PERSONS_TWO_PARENT, PERSON_LABELS_TWO_PARENT } from "../state/caseState";
import LocusGrid from "../components/LocusGrid";
import ParamsBar from "../components/ParamsBar";
import ResultsFooter from "../components/ResultsFooter";
import CaseIO from "../components/CaseIO";

export default function TwoParentPage() {
  const m = useModeCase({ mode: "two-parent", initialPersons: PERSONS_TWO_PARENT });
  return (
    <div>
      <h1>Two-parent kinship</h1>
      <p style={{ color: "var(--muted)", marginTop: 0 }}>
        H<sub>1</sub>: <code>Tested Parent</code> + <code>Known Parent</code> are the parents of <code>Child</code>.
        &nbsp;H<sub>2</sub>: only <code>Known Parent</code> is a parent of <code>Child</code>.
      </p>
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
        personLabels={PERSON_LABELS_TWO_PARENT}
        onChange={m.setLoci}
      />
      <CaseIO state={m.state} onLoad={m.setState} />
      <ResultsFooter result={m.combined} warnings={m.topWarnings} />
    </div>
  );
}
