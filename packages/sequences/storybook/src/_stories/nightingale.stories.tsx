/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { resolveView } from '@molstar/seq-core/format';
import type { Issue } from '@molstar/seq-core/model';
import { extentLoci, materialize, NightingaleSupported } from '@molstar/seq-nightingale';
import { mount, type NightingaleEvent } from '@molstar/seq-nightingale/mount';
import type { Meta, StoryObj } from './support/seqstar-story.js';
import { readPromise } from './support/seqstar-read.js';
import { Suspense, useEffect, useMemo, useRef, useState } from 'react';
import { CASE_TITLES, type CompiledCase, caseSource, compiledCase, preview, reg } from './support/seqstar-cases.js';
import s from './support/seqstar-format.module.css';
import { IssueList } from './support/seqstar-issues.js';
import { PrismPage } from './support/seqstar-prism.js';
import { AnchoredTooltip, Fold, ScrollArea, TextAction } from './ui/index.js';

/**
 * The second viewer that checks the format: each acceptance
 * case compiled on the test fixtures, each view resolved with Nightingale's supported set,
 * materialized as Nightingale's elements and mounted. Below each: what Nightingale cannot draw (the resolution's
 * fallbacks, materialize's losses and progressive tiles that failed) and the element description. Hover a feature
 * or residue: the event mapped back to block, index, layer, item and the item's extent (every fragment), lit in
 * every section through the host's loci. Block names and range labels are mount's layer over the tracks, following
 * zoom; one cut short or with no room is the hover's tooltip. Groups are header rows: click one to collapse or expand it in place; a
 * track's ⓘ (its data's provenance) and a label cut short say the rest in a tooltip.
 */
const meta = { title: 'Seq*/Cases/Nightingale', parameters: { layout: 'fullscreen' } } satisfies Meta;
export default meta;
type Story = StoryObj<typeof meta>;

function ViewCard({ c, view }: { c: CompiledCase; view: { ref: string; label: string } }) {
  const r = useMemo(() => {
    const { view: resolved, issues } = resolveView(c.doc, c.set, view.ref, reg, NightingaleSupported);
    const render = (i: Issue) => i.code.startsWith('svs.render.');
    const other = issues.filter((i) => !render(i));
    const m = materialize(resolved, c.set);
    // The fallbacks the resolution took for Nightingale's supported set, then what materialize could not draw.
    const lost = [...issues.filter(render), ...m.issues];
    return { view: resolved, root: m.root, lost, other };
  }, [c, view.ref]);
  const host = useRef<HTMLDivElement>(null);
  const [event, setEvent] = useState('');
  const [open, setOpen] = useState(false);
  const [failed, setFailed] = useState<Issue[]>([]);
  const [hint, setHint] = useState<{ anchor: Element; text: string } | null>(null);

  useEffect(() => {
    const el = host.current;
    const root = r.root;
    if (!el || !root) return;
    const set = c.set;
    const label = (b: number) => (b >= 0 ? set.strings.get(r.view.blocks.label[b]) : '–');
    const mounted = mount(root, el, { view: r.view, set }, (e: NightingaleEvent) => {
      const layer = e.layer >= 0 ? set.strings.get(r.view.layers.ref[e.layer]) : '–';
      // Block indices, 0-based, inclusive: each run of the item's extent.
      const runs = e.extent.map((x) => `${label(x.block)} ${x.start}${x.stop - x.start > 1 ? `–${x.stop - 1}` : ''}`);
      const extent = runs.length ? `, extent ${runs.join(', ')}` : '';
      setEvent(`${e.kind}: block ${label(e.block)}, index ${e.index}, layer ${layer}, item ${e.item}${extent}`);
      // A pointer hover by position (no item) closes in what contains it: a nucleotide lights its residue's codon.
      const closure = e.item < 0 ? 'containing' : 'element';
      mounted.highlight(e.kind === 'leave' ? null : extentLoci(r.view, set, e.extent), { closure });
      // A feature's label not shown in full (cut short, or no room at this zoom): the whole of it at the feature.
      if (e.kind === 'hover') setHint(e.hint);
      if (e.kind === 'leave') setHint(null);
    });
    // The label column is mount's DOM: a track's hint mark or a label cut short, under the pointer or focus (a
    // feature's label comes with its hover event, above).
    const hintAt = (e: Event) => {
      const t = e.target as Element;
      const mark = t.closest?.('[data-svs-hint]');
      const cell = t.closest?.('.svs-label-text');
      const cut = cell && cell.scrollWidth > cell.clientWidth ? cell : null;
      const anchor = mark ?? cut;
      const text = mark?.getAttribute('data-svs-hint') ?? cut?.textContent ?? '';
      setHint((h) => (h?.anchor === anchor ? h : anchor ? { anchor, text } : null));
    };
    const clear = () => setHint(null);
    el.addEventListener('pointerover', hintAt);
    el.addEventListener('focusin', hintAt);
    el.addEventListener('pointerleave', clear);
    el.addEventListener('focusout', clear);
    let live = true;
    mounted.ready.then((issues) => live && setFailed(issues));
    return () => {
      live = false;
      setFailed([]);
      el.removeEventListener('pointerover', hintAt);
      el.removeEventListener('focusin', hintAt);
      el.removeEventListener('pointerleave', clear);
      el.removeEventListener('focusout', clear);
      mounted.dispose();
    };
  }, [r, c.set]);

  return (
    <section className={s.card}>
      <h2 className={s.heading}>
        {view.label} <span className={s.code}>{view.ref}</span>
      </h2>
      <ScrollArea orientation="horizontal">
        <div ref={host} className={s.nightingale} />
      </ScrollArea>
      <AnchoredTooltip anchor={hint?.anchor ?? null} content={hint?.text} />
      <div className={s.status}>{event || 'Hover a feature or a residue.'}</div>
      <h3 className={s.subheading}>What Nightingale cannot draw</h3>
      <IssueList issues={failed.length ? [...r.lost, ...failed] : r.lost} />
      <h3 className={s.subheading}>Other issues of the resolution</h3>
      <IssueList issues={r.other} />
      <div className={s.row}>
        <TextAction tone="quiet" onClick={() => setOpen((o) => !o)}>
          {open ? 'Hide the element description' : 'Show the element description'}
        </TextAction>
      </div>
      <Fold open={open}>
        <ScrollArea className={s.box} orientation="both">
          <pre className={s.pre}>{open ? preview(r.root) : ''}</pre>
        </ScrollArea>
      </Fold>
    </section>
  );
}

function Views({ name }: { name: string }) {
  const c = readPromise(compiledCase(name));
  return c.views.map((v) => <ViewCard key={v.ref} c={c} view={v} />);
}

/** The case's builder function, the source the story runs (cases.ts, by its region markers), behind a toggle. */
function BuilderCode({ name }: { name: string }) {
  const [open, setOpen] = useState(false);
  const { code, helpers, fixtures } = useMemo(() => caseSource(name), [name]);
  const uses = [
    helpers.length > 0 && `from cases.ts: ${helpers.join(', ')}`,
    fixtures.length > 0 && `from the fixtures (load-fixtures.ts): ${fixtures.join(', ')}`,
  ].filter(Boolean);
  return (
    <div className={s.stack}>
      <div className={s.row}>
        <TextAction tone="quiet" onClick={() => setOpen((o) => !o)}>
          {open ? 'Hide the builder code' : 'Builder code'}
        </TextAction>
      </div>
      <Fold open={open}>
        <div className={s.stack}>
          <div className={s.code}>
            {uses.length > 0 ? `Uses ${uses.join('; ')}.` : 'Uses only the builder (@molstar/seq-core/builder).'}
          </div>
          <ScrollArea className={s.box} orientation="both">
            <pre className={s.pre}>{open ? code : ''}</pre>
          </ScrollArea>
        </div>
      </Fold>
    </div>
  );
}

function Case({ name }: { name: string }) {
  return (
    <PrismPage
      title={CASE_TITLES[name]}
      caption="Nightingale 5 (^5, from npm), through @molstar/seq-nightingale; fixtures, no network."
    >
      <BuilderCode name={name} />
      <Suspense fallback={<div className={s.status}>Compiling…</div>}>
        <Views name={name} />
      </Suspense>
    </PrismPage>
  );
}

export const Case1: Story = { name: 'p53: text and tracks', render: () => <Case name="case1" /> };
export const Case2: Story = { name: 'Antibody Fv, Kabat and IMGT', render: () => <Case name="case2" /> };
export const Case3: Story = { name: 'SH2 family', render: () => <Case name="case3" /> };
export const Case4: Story = { name: 'tRNA family, base pairs', render: () => <Case name="case4" /> };
export const Case5: Story = { name: 'DNA construct', render: () => <Case name="case5" /> };
export const Case6: Story = { name: 'Imatinib in ABL', render: () => <Case name="case6" /> };
export const Case7: Story = { name: 'Barnase–barstar', render: () => <Case name="case7" /> };
export const Case8: Story = { name: 'AlphaFold DB complex', render: () => <Case name="case8" /> };
export const Case9: Story = { name: 'PF00042 by view', render: () => <Case name="case9" /> };
export const Case10: Story = { name: 'The Stories drafts', render: () => <Case name="case10" /> };
export const Case12: Story = { name: 'Designed dimer', render: () => <Case name="case12" /> };
export const Case13: Story = { name: 'p53 features, UniProt-style', render: () => <Case name="case13" /> };
export const Case14: Story = { name: 'Protein over its coding DNA', render: () => <Case name="case14" /> };
export const Case15: Story = { name: 'Indels in protein designs', render: () => <Case name="case15" /> };
export const Case16: Story = { name: 'A plasmid map', render: () => <Case name="case16" /> };
export const Case17: Story = { name: 'A genotype matrix', render: () => <Case name="case17" /> };
export const Case18: Story = { name: 'Finding things in p53', render: () => <Case name="case18" /> };
export const Case19: Story = { name: 'A chromosome by region', render: () => <Case name="case19" /> };
export const Case20: Story = { name: 'A protein on its genome', render: () => <Case name="case20" /> };
export const Case21: Story = { name: 'A gene on its reverse strand', render: () => <Case name="case21" /> };
export const Case22: Story = { name: 'A family triaged', render: () => <Case name="case22" /> };
