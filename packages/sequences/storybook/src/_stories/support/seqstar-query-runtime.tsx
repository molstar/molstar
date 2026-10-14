/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The runtime gallery story's application: a gallery of example queries beside a
// mounted viewer, and a query box. A query runs on the viewer's current set through `@molstar/seq-core/model/query`; its result
// becomes the viewer's selection (or a highlight channel, the toggle), and "Show me" focuses on its extent, switching to
// the view that shows it when the current one shows none (`showExtents`). The viewer has no commands for queries yet: here the result goes through the existing channels, `selection`, `highlight` and focus. The result's groups
// can be stepped through, one at a time.

import type { SeqViewSpecDoc } from '@molstar/seq-core/format';
import { Extents, type SeqSet } from '@molstar/seq-core/model';
import { compileQuery, createQueryEnv, SeqSelection } from '@molstar/seq-core/model/query';
import {
  bindQuery,
  checkText,
  type Expression,
  printQuery,
  type QueryIssue,
  readQuery,
  renderIssue,
} from '@molstar/seq-core/query';
import { type SeqViewerContext, showExtents } from '@molstar/seq-viewer';
import { SeqViewerLayout, SeqViewerProvider, useSeqState } from '@molstar/seq-viewer/react';
import { type KeyboardEvent, type ReactNode, useEffect, useMemo, useRef, useState } from 'react';
import { readPromise } from './seqstar-read.js';
import { PrismPage } from '../support/seqstar-prism.js';
import { ScrollArea, Switch, TextAction, TextTabs, Tooltip } from '../ui/index.js';
import { type CompiledCase, compiledCase } from './seqstar-cases.js';
import { familyDoc } from './seqstar-docs-views.js';
import f from './seqstar-format.module.css';
import { PrismCanvas } from './seqstar-prism.js';
import q from './seqstar-query-runtime.module.css';
import { caseParts, useCaseContext, useStoryUi } from './seqstar-viewer.js';

/** The story's app registers its letter-class tables as data: `search :classes`. */
const ENV = createQueryEnv([
  {
    namespace: 'story',
    symbols: [],
    classes: { protein: { X: 'ACDEFGHIKLMNPQRSTVWY' }, iupac: { N: 'ACGT', R: 'AG', Y: 'CT', S: 'CG', W: 'AT' } },
  },
]);
const ORIGIN = 'host:query';
/** Highlight channels, one per loci the result gives (a unit's and an alignment's are two). */
const CHANNEL = 'query';

interface Entry {
  readonly text: string;
  readonly says: string;
}
interface CaseConfig {
  readonly name: string;
  readonly tab: string;
  readonly entries: readonly Entry[];
  /** A document of the story's own, in place of the case's. */
  readonly doc?: () => SeqViewSpecDoc;
}

let family: SeqViewSpecDoc | null = null;

// Cases 7, 1, 3 and 14, and the family: a chain range, an index pointer, an annotation's items, a union, values above
// a threshold; a pair's partners, positions onto the alignment and back, across codons, picked items, widened
// runs; letters, motifs, position and system tests, group tests and an empty result that says why.
const CASES: readonly CaseConfig[] = [
  {
    name: 'case7',
    tab: 'Barnase–barstar',
    entries: [
      { text: '(sel "1brs/A[10-20]")', says: 'residues 10 to 20 of barnase, chain A: a chain range in its labels' },
      { text: '(sel "1brs/A[10-20], 1brs/D[30-40]")', says: 'a site on each chain: a union, one group each' },
      { text: '(items :of contacts)', says: "each contact with both its ends: an annotation's items" },
      { text: '(and 1brs/A (not iface))', says: 'barnase chain A outside the interface' },
      {
        text: '(and 1brs/A (partners 1brs/D :via contacts))',
        says: "interface residues with their partner: barnase's residues in contact with barstar (chain D)",
      },
      { text: '(partners 1brs/A :via contacts)', says: "barstar's side of the same contacts, one group per run" },
      {
        text: '(systems :of 1brs :where (> (length) 100))',
        says: "the complex's chains longer than 100: barnase's three",
      },
      {
        text: '(pick (items :of contacts) :where (within (sel "1brs/A[80-90]") 2))',
        says: 'contacts within 2 residues of barnase 80-90: a group test, its selection read once',
      },
    ],
  },
  {
    name: 'family',
    tab: 'Barnase family',
    doc: () => {
      family ??= familyDoc();
      return family;
    },
    entries: [
      {
        text: '(onto (and 1brs/A (partners 1brs/D :via contacts)) msa)',
        says: "barnase's interface on the family's columns: Show me switches to the Family view",
      },
      {
        text: '(onto (values :of cons :above 0.8) 1brs/A)',
        says: "the family's conserved columns on chain A, through row 0's binding: back to the Complex view",
      },
      { text: '(onto (sel "msa[40-60]") msa :kind row)', says: "columns 40 to 60 on every row's residues" },
      { text: '(onto msa/3 1brs/A)', says: 'a row bound to no sequence: nothing, and a warning that says why' },
      { text: '(rows :of msa :where (matches (label) "HOM0[1-3]"))', says: 'rows by their labels: a system test' },
    ],
  },
  {
    name: 'case1',
    tab: 'p53',
    entries: [
      { text: '(sel "p53@0:50")', says: 'the first 50 positions, as indices a tool returned: an index pointer' },
      { text: '(items :of regions)', says: "p53's regions, one group each" },
      { text: '(and p53 (not regions))', says: 'what lies between the regions' },
      { text: '(range p53 175 245-249 :numbering "uniprot:P04637")', says: 'hotspots by their UniProt numbers' },
      {
        text: '(pick (items :of regions) :where (matches (label) binding))',
        says: "an item by its label, as kinase domains are picked: p53's DNA-binding region",
      },
      { text: '(first (items :of regions) 2)', says: 'the first two regions' },
      { text: '(expand (sel "p53[175]") 20)', says: '20 residues either side of residue 175' },
      { text: '(each (items :of regions) (expand (current) 5))', says: 'every region widened by 5, one group each' },
      { text: '(letters C :in p53)', says: "p53's cysteines, one group per run: step through them" },
      { text: '(positions :in p53 :where (in (letter) C H))', says: 'cysteines and histidines: a test per position' },
      {
        text: '(each (items :of regions) (and (current) (letters C)))',
        says: 'the cysteines of each region, one group each',
      },
    ],
  },
  {
    name: 'case3',
    tab: 'SH2 family',
    entries: [
      { text: '(values :of cons :above 0.8)', says: 'columns with conservation above 0.8, one group per run' },
      { text: '(rows :of sh2)', says: "the family's rows, one group each" },
      { text: '(sel "sh2/1[10-30]")', says: 'residues 10 to 30 of the second member' },
      { text: '(and (rows :of sh2) (letters W))', says: 'tryptophans in every row: a scan of the rows only' },
      { text: '(values :of cons :above 1)', says: 'nothing passes: the warning says what cons holds' },
    ],
  },
  {
    name: 'case13',
    tab: 'p53 features',
    entries: [
      {
        text: '(pick (items :of p53-variants) :where (overlaps (items :of p53-features :where (= (category) Domain))))',
        says: 'variants overlapping a domain: the hotspots in the DNA-binding domain',
      },
      { text: '(search PXXP :classes protein :in p53)', says: "a motif: PXXP, the proline-rich region's SH3 sites" },
      {
        text: '(search "[KR]{2}.{10,12}[KR]{3,}" :in p53 :regex true)',
        says: 'a regular expression (quoted): the bipartite nuclear localization signal, residues 305 to 321',
      },
      { text: '(items :of p53-features :where (= (category) Domain))', says: 'the features that are domains' },
      {
        text: '(search PPXP :classes protein :edits 1 :in p53)',
        says: 'PPXP with one edit allowed: each hit scored by its edits, as you step through them',
      },
    ],
  },
  {
    name: 'case14',
    tab: 'Gene and protein',
    entries: [
      { text: '(or exons domain)', says: "the exons and the protein's domain: a union of two annotations' items" },
      { text: '(and gene (not exons))', says: 'the gene outside its exons: the intron and the untranslated ends' },
      { text: '(sel "gene[1-12]")', says: 'the 5′ untranslated region by its labels' },
      { text: '(onto (sel "protein[5-10]") gene)', says: 'the codons of residues 5 to 10, across the codon mapping' },
      {
        text: '(onto (sel "gene[40-60]") protein)',
        says: 'the residues whose codons touch nucleotides 40 to 60 (by overlap)',
      },
      { text: '(search ATG :in gene)', says: 'the start codon, after the 5′ untranslated region' },
      { text: '(search TRA :classes iupac :in gene)', says: 'stop codons TAA and TGA, by IUPAC classes' },
    ],
  },
];

/** The text with each issue under its line: a caret at its column, then the issue. */
function inline(text: string, issues: readonly QueryIssue[]): string {
  const lines = text.split('\n');
  const out = lines.flatMap((line, k) => [
    line,
    ...issues
      .filter((i) => i.at?.line === k + 1)
      .map((i) => `${' '.repeat((i.at?.column ?? 1) - 1)}^ ${renderIssue(i).replace(/^\S+ /, '').split('\n')[0]}`),
  ]);
  return [...out, ...issues.filter((i) => !i.at).map((i) => renderIssue(i).split('\n')[0])].join('\n');
}

/** The expansion `sel` stands for, printed, when the text uses it. */
function expansion(text: string): string | null {
  const c = checkText(text);
  return text.includes('(sel ') && c.expanded ? printQuery(c.expanded as Expression) : null;
}

/** The person's selection as a query input (`(input selection)`). */
function selectionInput(ctx: SeqViewerContext, set: SeqSet): SeqSelection {
  const sel = ctx.selection.current.value;
  const loci = [sel.seq, sel.rows, sel.columns, sel.aligned, sel.items].flatMap((l) => (l && l.set === set ? [l] : []));
  return SeqSelection.fromLoci(set, loci);
}

interface Ran {
  readonly text: string;
  readonly selection: SeqSelection | null;
  readonly issues: readonly QueryIssue[];
}
/** Reads, binds, checks, compiles and runs a text on the set: its selection, or the issues that stopped it. */
function run(ctx: SeqViewerContext, set: SeqSet, text: string): Ran {
  const read = readQuery(text);
  const bound = read.node ? bindQuery(read.node, ENV.language) : null;
  if (!bound?.expr) return { text, selection: null, issues: [...read.issues, ...(bound?.issues ?? [])] };
  const r = compileQuery(bound.expr, ENV, bound.source).run({ set, inputs: { selection: selectionInput(ctx, set) } });
  const failed = r.issues.some((i) => i.severity === 'error');
  return { text, selection: failed ? null : r.selection, issues: [...bound.issues, ...r.issues] };
}

const count = (n: number, what: string) => `${n.toLocaleString('en-US')} ${what}${n === 1 ? '' : 's'}`;

/** Its systems: distinct (kind, id, row) of its runs, in canonical order. */
function systemCount(sel: SeqSelection): number {
  const u = SeqSelection.union(sel);
  let n = 0;
  for (let k = 0; k < u.runCount; k++) {
    if (k === 0 || u.kind[k] !== u.kind[k - 1] || u.id[k] !== u.id[k - 1] || u.row[k] !== u.row[k - 1]) n++;
  }
  return n;
}

/** Group g of a selection as one of its own, tagged as it was (the step's summary and "Show me"). */
function groupOf(sel: SeqSelection, g: number): SeqSelection {
  const ext = Extents.create(16);
  SeqSelection.toExtents(sel, ext, g);
  return SeqSelection.fromExtents(sel.set, ext);
}

/** The result into the viewer: the selection (the first loci set, the rest added), or one highlight channel each;
 * with `group`, that group only. */
function apply(
  ctx: SeqViewerContext,
  sel: SeqSelection | null,
  highlight: boolean,
  lit: number,
  group?: number,
): number {
  for (let i = 0; i < lit; i++) ctx.highlight.clear(`${CHANNEL}:${i}`);
  const loci = sel ? SeqSelection.toLoci(sel, { group }) : [];
  if (highlight) {
    for (const [i, l] of loci.entries()) ctx.highlight.set(`${CHANNEL}:${i}`, l, { style: 'highlight' });
    return loci.length;
  }
  ctx.selection.set(loci[0] ?? null, ORIGIN);
  for (const l of loci.slice(1)) ctx.selection.add(l, ORIGIN);
  return 0;
}

/** "Show me": focused on the result's extent, in the view that shows it (switched to); false when none does. */
function showMe(ctx: SeqViewerContext, sel: SeqSelection): Promise<boolean> {
  const ext = Extents.create(64);
  SeqSelection.toExtents(sel, ext);
  return showExtents(ctx, ext, { origin: ORIGIN });
}

function QueryPanel({ ctx, config }: { ctx: SeqViewerContext; config: CaseConfig }) {
  const set = useSeqState(ctx.graph.assembled)?.set ?? null;
  const [text, setText] = useState(config.entries[0].text);
  const [ran, setRan] = useState<Ran>({ text: '', selection: null, issues: [] });
  const [highlight, setHighlight] = useState(false);
  const [lit, setLit] = useState(0);
  const [shown, setShown] = useState<boolean | null>(null);
  /** The group stepped to; null: every group. */
  const [step, setStep] = useState<number | null>(null);
  const checked = useMemo(() => checkText(text, ENV.language), [text]);
  const go = (t: string) => {
    setText(t);
    setShown(null);
    setStep(null);
    if (!set) return;
    const r = run(ctx, set, t);
    setRan(r);
    setLit(apply(ctx, r.selection, highlight, lit));
  };
  const stepTo = (g: number | null) => {
    setStep(g);
    setShown(null);
    setLit(apply(ctx, ran.selection, highlight, lit, g ?? undefined));
  };
  // The first entry runs once the viewer's set is there.
  useEffect(() => {
    if (set) go(text);
  }, [set]);
  const toggle = (on: boolean) => {
    setHighlight(on);
    if (on) ctx.selection.clear(ORIGIN);
    setLit(apply(ctx, ran.selection, on, lit, step ?? undefined));
  };
  // The run's issues (the runtime's included) for the text it ran; the checker's while the text is being edited.
  const issues = ran.text === text ? ran.issues : checked.issues;
  const sel = ran.selection;
  const shownSel = sel && step !== null ? groupOf(sel, step) : sel;
  const groups = sel?.groupCount ?? 0;
  const expanded = useMemo(() => expansion(text), [text]);
  return (
    <div className={q.panel}>
      <Gallery entries={config.entries} active={text} onRun={go} />
      <section className={q.work}>
        <textarea
          className={q.box}
          value={text}
          onChange={(e) => setText(e.target.value)}
          onKeyDown={(e) => {
            if (e.key === 'Enter' && (e.metaKey || e.ctrlKey)) go(text);
          }}
          spellCheck={false}
          rows={2}
          aria-label="SeqQL query"
        />
        <div className={f.row}>
          <TextAction onClick={() => go(text)} disabled={!set || checked.expanded === null} title="Run (Ctrl+Enter)">
            Run
          </TextAction>
          <TextAction
            tone="quiet"
            disabled={!sel || sel.groupCount === 0}
            onClick={() => shownSel && void showMe(ctx, shownSel).then(setShown)}
            title="Focus on the result's extent, in the view that shows it"
          >
            Show me
          </TextAction>
          <Switch label="Highlight instead of select" checked={highlight} onCheckedChange={toggle} />
        </div>
        {expanded && <code className={q.expansion}>= {expanded}</code>}
        {issues.length > 0 && <pre className={q.issues}>{inline(text, issues)}</pre>}
        {sel && shownSel && (
          <div className={q.summary} aria-live="polite">
            <span className={q.line}>{SeqSelection.summary(shownSel, { indices: true, maxParts: 6 })}</span>
            {sel.groupCount > 0 && <span className={q.note}>in {count(systemCount(shownSel), 'system')}</span>}
            {step !== null && !Number.isNaN(sel.groupScore[step]) && (
              <span className={q.note}>score {sel.groupScore[step]}</span>
            )}
            {shown === false && <span className={q.note}>no view shows it</span>}
          </div>
        )}
        {groups > 1 && (
          <div className={f.row}>
            <TextAction tone="quiet" onClick={() => stepTo(step === null ? groups - 1 : (step + groups - 1) % groups)}>
              ‹ Previous
            </TextAction>
            <span className={q.note}>{step === null ? `all ${groups} groups` : `group ${step + 1} of ${groups}`}</span>
            <TextAction tone="quiet" onClick={() => stepTo(step === null ? 0 : (step + 1) % groups)}>
              Next ›
            </TextAction>
            {step !== null && (
              <TextAction tone="quiet" onClick={() => stepTo(null)}>
                All
              </TextAction>
            )}
          </div>
        )}
      </section>
    </div>
  );
}

/**
 * The example queries as compact cards in a grid (3 or 4 a row on a desktop, 1 on a phone), capped in height by a
 * scroll area so the box and the viewer stay in view. One roving group: the active card takes Tab, the arrows move
 * between cards, Home and End to the ends; a card's tooltip holds its whole text and what `sel` expands to.
 */
function Gallery(props: { entries: readonly Entry[]; active: string; onRun: (text: string) => void }) {
  const { entries, active, onRun } = props;
  const refs = useRef<(HTMLButtonElement | null)[]>([]);
  const current = Math.max(
    0,
    entries.findIndex((e) => e.text === active),
  );
  const move = (e: KeyboardEvent<HTMLButtonElement>, i: number) => {
    const n = entries.length;
    const to = { ArrowRight: i + 1, ArrowDown: i + 1, ArrowLeft: i - 1, ArrowUp: i - 1, Home: 0, End: n - 1 }[e.key];
    if (to === undefined) return;
    e.preventDefault();
    refs.current[(to + n) % n]?.focus();
  };
  return (
    <ScrollArea className={q.galleryArea} viewportClassName={q.gallery}>
      <div className={q.cards} role="toolbar" aria-label="Example queries">
        {entries.map((e, i) => {
          const sel = expansion(e.text);
          const tip = (
            <span className={q.tip}>
              <code>{e.text}</code>
              {sel && <code>= {sel}</code>}
              <span>{e.says}</span>
            </span>
          );
          return (
            <Tooltip key={e.text} content={tip} side="bottom" align="start">
              <button
                ref={(el) => {
                  refs.current[i] = el;
                }}
                type="button"
                className={q.entry}
                aria-current={e.text === active}
                tabIndex={i === current ? 0 : -1}
                onClick={() => onRun(e.text)}
                onKeyDown={(k) => move(k, i)}
              >
                <code className={q.code}>{e.text}</code>
                <span className={q.says}>{e.says}</span>
              </button>
            </Tooltip>
          );
        })}
      </div>
    </ScrollArea>
  );
}

/** One case: its viewer on the fixtures, the gallery and the box above it. */
function CaseQueries({ config }: { config: CaseConfig }) {
  const compiled = config.doc ? null : readPromise(compiledCase(config.name));
  const doc = config.doc ? config.doc() : (compiled as CompiledCase).doc;
  const m = useCaseContext(doc, { parts: compiled ? caseParts(compiled) : undefined });
  const ui = useStoryUi(m);
  return (
    <section className={f.card}>
      {m ? (
        <SeqViewerProvider context={m.ctx} ui={ui}>
          <QueryPanel ctx={m.ctx} config={config} />
          <PrismCanvas>
            <SeqViewerLayout />
          </PrismCanvas>
        </SeqViewerProvider>
      ) : (
        <div className={f.status}>Loading…</div>
      )}
    </section>
  );
}

/** The story: the gallery and the box over cases 1, 3, 7 and 14, and the barnase family. */
export function QueryStory(): ReactNode {
  const [name, setName] = useState(CASES[0].name);
  const config = CASES.find((c) => c.name === name) ?? CASES[0];
  return (
    <PrismPage
      title="Queries"
      caption={
        <>
          SeqQL queries run on the viewer's set through <code>@molstar/seq-core/model/query</code>: click one, or write
          your own (Ctrl+Enter runs it). The result becomes the selection, or a highlight; "Show me" focuses the view on
          it.
        </>
      }
    >
      <ScrollArea orientation="horizontal">
        <TextTabs tabs={CASES.map((c) => ({ value: c.name, label: c.tab }))} value={name} onValueChange={setName} />
      </ScrollArea>
      <CaseQueries key={name} config={config} />
    </PrismPage>
  );
}
