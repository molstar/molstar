/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The Seq*/Query and search/Search story's application: a search box beside a mounted
// viewer. Everything a search means is decided here, not in the library or the viewer: which sequences are searched, what
// the letters of a degenerate query stand for, what a hit is lit as (its coordinates in the sequence's coordinate system,
// through the viewer's highlight, selection and focus), and how a hit is ranked and listed.

import {
  type ApproximateEnds,
  approximateMatcher,
  classMatcher,
  compileMotif,
  exactMatcher,
  type Matcher,
  regexMatcher,
  toCodes,
} from '@molstar/seq-core/algorithms/search';
import { KeyKind, type SeqSet } from '@molstar/seq-core/model';
import { type Coordinates, SeqCommands, type SeqViewerContext } from '@molstar/seq-viewer';
import { SeqViewerLayout, SeqViewerProvider } from '@molstar/seq-viewer/react';
import { type KeyboardEvent, type ReactNode, use, useEffect, useMemo, useState } from 'react';
import { PrismPage } from '../support/seqstar-prism.js';
import { ScrollArea, TextAction, type TextTab, TextTabs } from '../ui/index.js';
import { compiledCase } from './seqstar-cases.js';
import f from './seqstar-format.module.css';
import { PrismCanvas } from './seqstar-prism.js';
import q from './seqstar-search.module.css';
import { caseParts, useCaseContext, useStoryUi } from './seqstar-viewer.js';

type Mode = 'exact' | 'degenerate' | 'regex' | 'approximate';
type Alphabet = 'dna' | 'protein';

/** IUPAC's nucleotide codes: the application's data, passed to `compileMotif`. */
const IUPAC: Readonly<Record<string, string>> = {
  R: 'AG',
  Y: 'CT',
  S: 'CG',
  W: 'AT',
  K: 'GT',
  M: 'AC',
  B: 'CGT',
  D: 'AGT',
  H: 'ACT',
  V: 'ACG',
  N: 'ACGT',
};
/** A protein grouping of the story's own: B, Z, J and X as in IUPAC, and signs for charge and hydrophobicity. */
const GROUPS: Readonly<Record<string, string>> = {
  B: 'DN',
  Z: 'EQ',
  J: 'IL',
  X: 'ACDEFGHIKLMNPQRSTVWY',
  '+': 'KRH',
  '-': 'DE',
  h: 'AVILMFWC',
};
const CLASSES: Record<Alphabet, Readonly<Record<string, string>>> = { dna: IUPAC, protein: GROUPS };

interface Target {
  readonly key: string;
  readonly label: string;
  readonly alphabet: Alphabet;
  readonly examples: Readonly<Record<Mode, string>>;
}
interface CaseConfig {
  readonly name: string;
  readonly tab: string;
  readonly targets: readonly Target[];
}
const CASES: readonly CaseConfig[] = [
  {
    name: 'case1',
    tab: 'p53 (protein)',
    targets: [
      {
        key: 'p53',
        label: 'p53',
        alphabet: 'protein',
        examples: { exact: 'RHSVV', degenerate: '+HSVV', regex: '[KR]{2}.{10,12}[KR]{3,}', approximate: 'RHSWV' },
      },
    ],
  },
  {
    name: 'case14',
    tab: 'Gene and protein',
    targets: [
      {
        key: 'gene',
        label: 'DNA',
        alphabet: 'dna',
        examples: { exact: 'ATGGCT', degenerate: 'GCNNCC', regex: 'GC..CC', approximate: 'ATGGCC' },
      },
      {
        key: 'protein',
        label: 'Protein',
        alphabet: 'protein',
        examples: { exact: 'ELT', degenerate: 'K-LT', regex: 'K[DE]LT', approximate: 'KDLTP' },
      },
    ],
  },
];
const MODES: readonly TextTab[] = [
  { value: 'exact', label: 'Exact', hint: "A Horspool search for the pattern's letters." },
  {
    value: 'degenerate',
    label: 'Degenerate',
    hint: "A set of letters per position, from the application's table: IUPAC codes for DNA (R is A or G), a grouping for proteins (+ basic, - acidic, h hydrophobic, B, Z, J, X).",
  },
  {
    value: 'regex',
    label: 'Regular expression',
    hint: "A JavaScript regular expression over the letters, as SeqQL's (search P :regex true): [KR] one of two, . any, {2,4} a count, ^ not in a class. Every start is tried, so hits overlap.",
  },
  {
    value: 'approximate',
    label: 'Approximate',
    hint: "Myers' bit-parallel search: the pattern with at most k edits (substitutions, insertions or deletions).",
  },
];

/** The molecule's letters in the set, or null. */
function moleculeCodes(set: SeqSet, key: string): Uint8Array | null {
  const out = { kind: 0, id: 0, row: 0 };
  return set.byName(key, out) && out.kind === KeyKind.molecule ? (set.molecules[out.id]?.code ?? null) : null;
}

type HitPolicy = 'best' | 'all';
type Tie = 'last' | 'first';
const ENDS: readonly TextTab[] = [
  {
    value: 'best',
    label: 'Best per site',
    hint: 'The ends whose edit count is a local minimum along the sequence: one hit per site, every exact match included.',
  },
  {
    value: 'all',
    label: 'Every end',
    hint: 'Every end with at most k edits: an edit at a site makes its neighbors hits too.',
  },
];
const TIES: readonly TextTab[] = [
  {
    value: 'last',
    label: 'Last end',
    hint: "Of a run of ends with the same edit count, the last: a mismatch at the pattern's end reads as a substitution (RHSWV finds RHSVV).",
  },
  {
    value: 'first',
    label: 'First end',
    hint: "Of a run of ends with the same edit count, the first: a mismatch at the pattern's end reads as a deletion (RHSWV finds RHSV).",
  },
];
/** The story's two controls as the library's one option. */
const endsOf = (hits: HitPolicy, tie: Tie): ApproximateEnds =>
  hits === 'all' ? 'all' : tie === 'last' ? 'local_best_last' : 'local_best_first';

/** The matcher for the application's reading of the box, or a word on why there is none. */
function matcherFor(
  mode: Mode,
  text: string,
  alphabet: Alphabet,
  edits: number,
  ends: ApproximateEnds,
): Matcher | string | null {
  if (text.length === 0) return null;
  try {
    if (mode === 'degenerate') {
      // Lower case stands for itself unless the table gives the letter a meaning ("h"); letters are upper case otherwise.
      const motif = Array.from(text, (c) => (Object.hasOwn(CLASSES[alphabet], c) ? c : c.toUpperCase())).join('');
      return classMatcher(compileMotif(motif, CLASSES[alphabet]));
    }
    if (mode === 'regex') return regexMatcher(new RegExp(text, 'i'));
    const pattern = toCodes(text.toUpperCase());
    if (mode === 'exact') return exactMatcher(pattern);
    return pattern.length <= edits
      ? `${edits} edits need a pattern of ${edits + 1} letters or more`
      : approximateMatcher(pattern, edits, { ends });
  } catch (e) {
    return e instanceof Error ? e.message : String(e);
  }
}

/** The hits of one search: triples (start, stop, score) in the sequence's indices. */
interface Hits {
  readonly runs: Int32Array;
  readonly count: number;
  readonly note: string | null;
}
const NO_HITS: Hits = { runs: new Int32Array(0), count: 0, note: null };

function search(codes: Uint8Array | null, matcher: Matcher | string | null): Hits {
  if (matcher === null || !codes) return NO_HITS;
  if (typeof matcher === 'string') return { ...NO_HITS, note: matcher };
  const found: number[] = [];
  matcher.find(codes, 0, codes.length, false, (start, stop, score) => found.push(start, stop, score));
  return { runs: Int32Array.from(found), count: found.length / 3, note: null };
}

/** The hits as `Coordinates` in the molecule's coordinate system the application's reading of a hit. */
function coordinatesOf(key: string, hits: Hits, only?: number): Coordinates {
  if (only !== undefined)
    return { systems: [`seqstar:${key}`], runs: Int32Array.of(0, hits.runs[only * 3], hits.runs[only * 3 + 1]) };
  const runs = new Int32Array(hits.count * 3);
  for (let i = 0; i < hits.count; i++) runs.set([0, hits.runs[i * 3], hits.runs[i * 3 + 1]], i * 3);
  return { systems: [`seqstar:${key}`], runs };
}

const CHANNEL = 'search';
const CURRENT = 'search:current';

/** Lights the hits (a highlight), the current one (an accent), and scrolls the viewer to it. */
function light(ctx: SeqViewerContext, key: string, hits: Hits, current: number, focus: boolean) {
  if (hits.count === 0) {
    ctx.highlight.clear(CHANNEL);
    ctx.highlight.clear(CURRENT);
    return;
  }
  ctx.highlight.set(CHANNEL, coordinatesOf(key, hits), { style: 'highlight' });
  const one = coordinatesOf(key, hits, current);
  ctx.highlight.set(CURRENT, one, { style: 'accent' });
  if (focus) ctx.focus.fromLoci(one, { origin: 'api' });
}

function Hit({
  codes,
  hits,
  i,
  on,
  onPick,
}: {
  codes: Uint8Array;
  hits: Hits;
  i: number;
  on: boolean;
  onPick: () => void;
}) {
  const [start, stop, score] = [hits.runs[i * 3], hits.runs[i * 3 + 1], hits.runs[i * 3 + 2]];
  const letters = String.fromCharCode(...codes.subarray(start, stop));
  return (
    <button type="button" className={q.hit} aria-current={on} onClick={onPick}>
      <span className={q.at}>{start + 1}</span>
      <span className={q.letters}>{letters}</span>
      {score > 0 && <span className={q.edits}>{score === 1 ? '1 edit' : `${score} edits`}</span>}
    </button>
  );
}

/** The search box and its hits beside the viewer on `ctx`; `set` the case's compiled set. */
function SearchPanel({ ctx, set, config }: { ctx: SeqViewerContext; set: SeqSet; config: CaseConfig }) {
  const [targetKey, setTargetKey] = useState(config.targets[0].key);
  const target = config.targets.find((t) => t.key === targetKey) ?? config.targets[0];
  const [mode, setMode] = useState<Mode>('exact');
  const [edits, setEdits] = useState(1);
  const [hitPolicy, setHitPolicy] = useState<HitPolicy>('best');
  const [tie, setTie] = useState<Tie>('last');
  const ends = endsOf(hitPolicy, tie);
  const [text, setText] = useState(target.examples.exact);
  const [current, setCurrent] = useState(0);
  const [focusNext, setFocusNext] = useState(false);

  const codes = useMemo(() => moleculeCodes(set, target.key), [set, target.key]);
  const matcher = useMemo(
    () => matcherFor(mode, text, target.alphabet, edits, ends),
    [mode, text, target.alphabet, edits, ends],
  );
  const hits = useMemo(() => search(codes, matcher), [codes, matcher]);
  const at = hits.count === 0 ? 0 : Math.min(current, hits.count - 1);

  // The hits light on every change of the search; scrolling to the current one only when the person steps or picks.
  // biome-ignore lint/correctness/useExhaustiveDependencies: `at` and `focusNext` are read for the focus, not re-lit by it
  useEffect(() => {
    light(ctx, target.key, hits, at, focusNext);
    setFocusNext(false);
  }, [ctx, target.key, hits, at]);
  useEffect(
    () => () => {
      ctx.highlight.clear(CHANNEL);
      ctx.highlight.clear(CURRENT);
    },
    [ctx],
  );

  const step = (dir: 1 | -1) => {
    if (hits.count === 0) return;
    setFocusNext(true);
    setCurrent((at + dir + hits.count) % hits.count);
  };
  const onKey = (e: KeyboardEvent<HTMLInputElement>) => {
    if (e.key === 'Enter') step(e.shiftKey ? -1 : 1);
  };
  const choose = (m: Mode) => {
    setMode(m);
    setText(target.examples[m]);
    setCurrent(0);
  };
  const chooseTarget = (key: string) => {
    const next = config.targets.find((t) => t.key === key) ?? target;
    setTargetKey(key);
    setText(next.examples[mode]);
    setCurrent(0);
  };
  const select = () => {
    if (hits.count > 0)
      void ctx.commands.dispatch(SeqCommands.Select, {
        target: { coordinates: coordinatesOf(target.key, hits) },
        origin: 'host:search',
      });
  };

  return (
    <div className={q.panel}>
      <div className={q.controls}>
        {config.targets.length > 1 && (
          <div className={q.group}>
            <span className={q.caption}>Search in</span>
            <TextTabs
              tabs={config.targets.map((t) => ({ value: t.key, label: t.label }))}
              value={target.key}
              onValueChange={chooseTarget}
            />
          </div>
        )}
        <div className={q.group}>
          <span className={q.caption}>Match</span>
          <TextTabs tabs={MODES as TextTab[]} value={mode} onValueChange={(v) => choose(v as Mode)} />
        </div>
        {mode === 'approximate' && (
          <label className={q.edits}>
            edits
            <input
              type="range"
              min={0}
              max={3}
              value={edits}
              onChange={(e) => setEdits(Number(e.target.value))}
              aria-label="Allowed edits"
            />
            <b>{edits}</b>
          </label>
        )}
        {mode === 'approximate' && (
          <div className={q.group}>
            <span className={q.caption}>Hits</span>
            <TextTabs
              tabs={ENDS as TextTab[]}
              value={hitPolicy}
              onValueChange={(v) => {
                setHitPolicy(v as HitPolicy);
                setCurrent(0);
              }}
            />
          </div>
        )}
        {mode === 'approximate' && hitPolicy === 'best' && (
          <div className={q.group}>
            <span className={q.caption}>Ties</span>
            <TextTabs
              tabs={TIES as TextTab[]}
              value={tie}
              onValueChange={(v) => {
                setTie(v as Tie);
                setCurrent(0);
              }}
            />
          </div>
        )}
      </div>
      <div className={q.box}>
        <input
          className={q.input}
          type="search"
          value={text}
          spellCheck={false}
          autoComplete="off"
          placeholder={`Search ${target.label}`}
          aria-label="Pattern"
          onChange={(e) => {
            setText(e.target.value);
            setCurrent(0);
          }}
          onKeyDown={onKey}
        />
        <span className={q.count} aria-live="polite">
          {matcher === null
            ? 'type a pattern'
            : hits.note
              ? hits.note
              : hits.count === 0
                ? 'no hits'
                : `${at + 1} of ${hits.count}`}
        </span>
        <TextAction onClick={() => step(-1)} disabled={hits.count === 0} title="Previous hit (Shift+Enter)">
          Previous
        </TextAction>
        <TextAction onClick={() => step(1)} disabled={hits.count === 0} title="Next hit (Enter)">
          Next
        </TextAction>
        <TextAction tone="quiet" onClick={select} disabled={hits.count === 0} title="Select every hit in the viewer">
          Select all
        </TextAction>
        <TextAction
          tone="quiet"
          onClick={() => ctx.selection.clear('host:search')}
          title="Clear the viewer's selection"
        >
          Clear selection
        </TextAction>
      </div>
      {codes && hits.count > 0 && (
        <ScrollArea className={q.hits} viewportClassName={q.hitsViewport}>
          <div className={q.hitList}>
            {Array.from({ length: Math.min(hits.count, 60) }, (_, i) => (
              <Hit
                key={`${hits.runs[i * 3]}-${hits.runs[i * 3 + 1]}`}
                codes={codes}
                hits={hits}
                i={i}
                on={i === at}
                onPick={() => {
                  setFocusNext(true);
                  setCurrent(i);
                }}
              />
            ))}
            {hits.count > 60 && <span className={q.more}>and {hits.count - 60} more</span>}
          </div>
        </ScrollArea>
      )}
    </div>
  );
}

/** One case: its viewer on the fixtures, the search panel above it. */
function CaseSearch({ config }: { config: CaseConfig }) {
  const c = use(compiledCase(config.name));
  const m = useCaseContext(c.doc, { parts: caseParts(c) });
  const ui = useStoryUi(m);
  return (
    <section className={f.card}>
      <h2 className={f.heading}>{config.tab}</h2>
      {m ? (
        <SeqViewerProvider context={m.ctx} ui={ui}>
          <SearchPanel ctx={m.ctx} set={c.set} config={config} />
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

/** The story: the application's search over the p53 case and the gene case. */
export function SearchStory(): ReactNode {
  const [name, setName] = useState(CASES[0].name);
  const config = CASES.find((c) => c.name === name) ?? CASES[0];
  return (
    <PrismPage
      title="A search the application composes"
      caption={
        <>
          <code>@molstar/seq-core/algorithms/search</code> finds the hits; this story decides what they mean: it lights
          them through the viewer's highlight channels, accents the current one and focuses it, and selects them on
          request. Neither the format nor the viewer has a search.
        </>
      }
    >
      <TextTabs tabs={CASES.map((c) => ({ value: c.name, label: c.tab }))} value={name} onValueChange={setName} />
      <CaseSearch key={name} config={config} />
    </PrismPage>
  );
}
