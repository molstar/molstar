/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `Seq*/Cases/Mol*` demo: a Mol* canvas (a plugin, no Mol* UI) beside the Seq*
// viewer over the p53 document (p53-demo.ts). What 3D shows is state here, `Shown`: a structure and the track coloring
// it (null: its own default, pLDDT for the AlphaFold model, none for a PDB entry). Each change composes an MVS state
// from the viewer's state at that moment (`composeStructureView`) and loads it; the structure is then linked to the
// viewer: hover and selection both ways, a selection focusing the other view. "Show in 3D" comes from the hover card
// and the track menu, "Apply to 3D" from the track menu and the chooser over the canvas. SeqQL over the linked
// structure (QM1): `molql` runs a PyMOL selection in Mol*, Mol*'s selection is `(input structure)`, a result lights the
// sequence and goes back to 3D. Offline: the fixtures are the loaders' bytes and Mol*'s assets at the services' URLs.

import { compileQuery, createQueryEnv, SeqSelection } from '@molstar/seq-core/model/query';
import { bindQuery, readQuery, renderIssue } from '@molstar/seq-core/query';
import {
  MolstarSvsEntry,
  molstarQueryExtension,
  StructureLink,
  type StructureRef,
  selectionLoci,
  structureInput,
  structureReach,
  structureUrl,
} from '@molstar/seq-molstar-poc';
import {
  mountPlugin,
  type PluginContext,
  pluginView3D,
  selectedLoci,
  setBackground,
  showMvs,
} from '@molstar/seq-molstar-poc/plugin';
import {
  colorByItemAction,
  colorTracks,
  composeStructureView,
  linkViews,
  type Shown,
  type StructureView3D,
  showInStructureAction,
  structureTrackActions,
} from '@molstar/seq-molstar-poc/viewer';
import {
  createSeqViewer,
  DefaultRegistry,
  DefaultSeqViewerSpec,
  SeqViewerBehaviors,
  SeqViewerConfig,
  type SeqViewerContext,
  SeqViewerSpec,
} from '@molstar/seq-viewer';
import { windowScheduler } from '@molstar/seq-viewer/canvas';
import { SeqViewerLayout, SeqViewerProvider } from '@molstar/seq-viewer/react';
import { useEffect, useMemo, useRef, useState } from 'react';
import { P53 } from '../../../../core/_test/format/load-fixtures.js';
import features from '../../../../molstar-poc/_test/fixtures/p53/P04637.json?raw';
import ptms from '../../../../molstar-poc/_test/fixtures/p53/P04637-ptms.tsv?raw';
import variants from '../../../../molstar-poc/_test/fixtures/p53/variation-P04637.json?raw';
import c1c26 from '../../../../molstar-poc/_test/fixtures/structures/1c26.bcif?url';
import c1tup from '../../../../molstar-poc/_test/fixtures/structures/1tup.bcif?url';
import c2ocj from '../../../../molstar-poc/_test/fixtures/structures/2ocj.bcif?url';
import c2xwr from '../../../../molstar-poc/_test/fixtures/structures/2xwr.bcif?url';
import caf from '../../../../molstar-poc/_test/fixtures/structures/AF-P04637-F1-model_v6.bcif?url';
import {
  coverageTrack,
  P53_STRUCTURES,
  type P53Structure,
  p53DemoDoc,
  p53Loaders,
  p53ShownPatch,
} from '../../../../molstar-poc/_test/p53-demo.js';
import { TextAction } from '../ui/index.js';
import f from './seqstar-format.module.css';
import m from './seqstar-molstar.module.css';
import { storyPatch } from './seqstar-patch.js';
import { PrismCanvas } from './seqstar-prism.js';
import q from './seqstar-query-runtime.module.css';
import { partSlots, STORY_CHROME, useStoryUi } from './seqstar-viewer.js';

/** The p53 demo's `headerHover` limit, in its viewer spec: a header of more items lights its row band alone. */
const HEADER_MAX_ITEMS = 500;
const FIXTURES: Record<string, string> = {
  '1tup': c1tup,
  '2ocj': c2ocj,
  '2xwr': c2xwr,
  '1c26': c1c26,
  'AF-P04637-F1': caf,
};
const refOf = (s: P53Structure): StructureRef => ({
  source: s.source,
  ...('assembly' in s ? { assembly: s.assembly } : {}),
});
const sourceId = (s: StructureRef | null | undefined) => (s && 'id' in s.source ? s.source.id : null);
const fasta = (seq: string) => `>sp|P04637|P53_HUMAN Cellular tumor antigen p53\n${seq.replace(/(.{60})/g, '$1\n')}\n`;

/** Each structure's bytes, fetched once from the bundled fixtures, by the URL its service resolves to. */
async function structureBytes(): Promise<Map<string, Uint8Array>> {
  const out = new Map<string, Uint8Array>();
  for (const s of P53_STRUCTURES) {
    const at = structureUrl({ source: s.source });
    const res = await fetch(FIXTURES[s.source.id]);
    if (at) out.set(at.url, new Uint8Array(await res.arrayBuffer()));
  }
  return out;
}

/** The Mol* canvas on `host` once the fixtures are in (`mountPlugin`); null until mounted. */
function usePlugin(host: React.RefObject<HTMLDivElement | null>, bytes: Map<string, Uint8Array> | null) {
  const [plugin, setPlugin] = useState<PluginContext | null>(null);
  useEffect(() => {
    const el = host.current;
    if (!el || !bytes) return;
    let live = true;
    let made: PluginContext | null = null;
    void mountPlugin(el, { files: bytes }).then((p) => {
      made = p;
      if (live) setPlugin(p);
      else p.dispose();
    });
    return () => {
      live = false;
      setPlugin(null);
      made?.dispose();
    };
  }, [host, bytes]);
  return plugin;
}

/** The viewer: the document, "Show in 3D" in its hover card and its track menu, "Apply to 3D". */
function useViewer(
  bytes: Map<string, Uint8Array> | null,
  shown: { readonly current: Shown | null },
  show: (s: Shown) => void,
) {
  const ctx = useMemo(() => {
    if (!bytes) return null;
    const item = showInStructureAction(show, () => shown.current);
    const color = colorByItemAction(() => shown.current, show);
    const tracks = structureTrackActions({
      structures: P53_STRUCTURES.map((s) => ({ label: s.label, structure: refOf(s) })),
      shown: () => shown.current,
      show,
    });
    const builtIn = SeqViewerConfig.UI.TrackActions.defaultValue;
    // The viewer's reactions are its spec's behaviors (seqspec:viewer#behaviors): here a header's hover lights at most
    // `HEADER_MAX_ITEMS` items (the default is 5,000), past it the row band alone.
    const spec = DefaultSeqViewerSpec();
    const behaviors = spec.behaviors.map((b) =>
      b.behavior === SeqViewerBehaviors.HeaderHover
        ? SeqViewerSpec.Behavior(SeqViewerBehaviors.HeaderHover, { maxItems: HEADER_MAX_ITEMS })
        : b,
    );
    const c = createSeqViewer(
      {
        ...spec,
        behaviors,
        registry: [DefaultRegistry, { svs: MolstarSvsEntry, actions: [item, color] }],
        config: [
          [SeqViewerConfig.UI.Chrome, STORY_CHROME],
          [SeqViewerConfig.UI.ItemActions, [item.name, color.name]],
          [SeqViewerConfig.UI.TrackActions, [...(builtIn === 'none' ? [] : builtIn), ...tracks]],
        ],
      },
      {
        loaders: p53Loaders({ structures: bytes, fasta: fasta(P53.sequence), features, variants, ptms }),
        scheduler: windowScheduler(),
      },
    );
    void c.setDocument(p53DemoDoc(P53.sequence));
    return c;
  }, [bytes]);
  useEffect(() => () => ctx?.dispose(), [ctx]);
  return ctx;
}

/**
 * The PDB coverage entry of the structure shown in 3D (`show`): marked "· in 3D" (a patch) and opened, with its
 * group; the entry the story opened before is closed again unless the person toggled it since.
 */
function shownEntry(ctx: SeqViewerContext) {
  const opened = new Set<string>();
  let ours = false;
  let before = ctx.ui.expanded.value;
  // A toggle not ours: the person's, so that entry stays as they left it.
  const sub = ctx.ui.expanded.subscribe((now) => {
    if (!ours) for (const r of opened) if (before.has(r) !== now.has(r)) opened.delete(r);
    before = now;
  });
  const toggle = (ref: string) => {
    ours = true;
    ctx.ui.toggleGroup(ref);
    ours = false;
  };
  return {
    async show(s: P53Structure | null) {
      await storyPatch(ctx, p53ShownPatch(s), 'shown in 3D');
      const want = s ? ['g-pdb', coverageTrack(s)] : [];
      for (const r of [...opened]) {
        if (want.includes(r)) continue;
        if (ctx.ui.expanded.value.has(r)) toggle(r);
        opened.delete(r);
      }
      for (const r of want) {
        if (ctx.ui.expanded.value.has(r)) continue;
        toggle(r);
        opened.add(r);
      }
    },
    dispose: () => sub.unsubscribe(),
  };
}

/** "Colored by" over the canvas: the tracks that reach the shown structure, its own default first. */
function ColorChooser({ ctx, shown, show }: { ctx: SeqViewerContext; shown: Shown; show: (s: Shown) => void }) {
  const options = useMemo(() => colorTracks(ctx, shown.structure), [ctx, shown.structure]);
  return (
    <label className={m.chooser}>
      Colored by{' '}
      <select
        value={shown.colorBy ?? ''}
        onChange={(e) => show({ ...shown, colorBy: e.target.value === '' ? null : e.target.value, where: null })}
      >
        <option value="">Its default</option>
        {options.map((o) => (
          <option key={o.track} value={o.track}>
            {o.label}
          </option>
        ))}
      </select>
    </label>
  );
}

/** The structure shown in 3D, linked: what `molql` and `(input structure)` read, and where a result goes back. */
interface Linked {
  readonly plugin: PluginContext;
  readonly link: StructureLink;
  readonly view3d: StructureView3D;
}

const MOLQL_DNA = '(molql "polymer.protein and (polymer.nucleic around 4)")';
// Results stay on the chains (instance space) unless a query projects them onto P04637, which joins the chains.
const PRESETS = [
  { label: 'Mol* selection', text: '(input structure)', says: 'select residues in 3D first' },
  { label: 'Zinc site', text: '(molql "polymer and (resn ZN around 3)")', says: 'within 3 Å of a zinc' },
  { label: 'DNA contacts', text: MOLQL_DNA, says: 'within 4 Å of DNA (1TUP)' },
  { label: 'DNA contacts on P04637', text: `(onto ${MOLQL_DNA} p53)`, says: "any chain's, on the sequence" },
  {
    label: 'Variants at the DNA',
    text: `(pick (items :of variants :where (= (length) 1)) :where (overlaps (onto ${MOLQL_DNA} p53)))`,
    says: 'single-residue variants on them',
  },
];

/** SeqQL over the linked structure: presets and a box; the result lights the sequence view, "Select in 3D" sends it. */
function MolQueries({ ctx, linked }: { ctx: SeqViewerContext; linked: { readonly current: Linked | null } }) {
  const env = useMemo(() => createQueryEnv([molstarQueryExtension(() => linked.current?.link ?? null)]), [linked]);
  const [text, setText] = useState(PRESETS[1].text);
  const [ran, setRan] = useState<{ text: string; sel: SeqSelection | null; said: string[]; reach: string } | null>(
    null,
  );
  const lit = useRef(0);
  const run = (t: string) => {
    setText(t);
    const set = ctx.graph.assembled.value?.set;
    if (!set) return;
    const read = readQuery(t);
    const bound = read.node ? bindQuery(read.node, env.language) : null;
    const issues = [...read.issues, ...(bound?.issues ?? [])];
    let sel: SeqSelection | null = null;
    if (bound?.expr) {
      const l = linked.current?.link.set === set ? linked.current : null;
      const structure = l ? structureInput(l.link, selectedLoci(l.plugin, l.link.structure)) : SeqSelection.empty(set);
      const r = compileQuery(bound.expr, env, bound.source).run({ set, inputs: { structure } });
      issues.push(...r.issues);
      if (!r.issues.some((i) => i.severity === 'error')) sel = r.selection;
    }
    for (let i = 0; i < lit.current; i++) ctx.highlight.clear(`molql:${i}`);
    const loci = sel ? SeqSelection.toLoci(sel) : [];
    for (const [i, l] of loci.entries()) ctx.highlight.set(`molql:${i}`, l, { style: 'highlight' });
    lit.current = loci.length;
    const l = linked.current;
    const r = l && sel?.groupCount ? structureReach(l.link, sel) : null;
    const n = r?.chains.length ?? 0;
    const which = `chain${n === 1 ? '' : 's'} ${r?.chains.join(', ')}`;
    const reach = !r || n === 0 ? '' : r.exact ? `in 3D: ${which}` : `on the sequence: any of ${which} in 3D`;
    setRan({ text: t, sel, said: issues.map((i) => renderIssue(i).split('\n')[0]), reach });
  };
  const clear = () => {
    for (let i = 0; i < lit.current; i++) ctx.highlight.clear(`molql:${i}`);
    lit.current = 0;
    setRan(null);
  };
  const to3d = () => {
    const l = linked.current;
    if (!l || !ran?.sel) return;
    const loci = selectionLoci(l.link, ran.sel);
    l.view3d.select(loci);
    if (loci.elements.length) l.view3d.focus(loci);
  };
  return (
    <section className={f.card}>
      <div className={m.switch} role="toolbar" aria-label="SeqQL over Mol*">
        {PRESETS.map((p) => (
          <button
            key={p.label}
            type="button"
            className={m.choice}
            aria-pressed={ran?.text === p.text}
            onClick={() => (ran?.text === p.text ? clear() : run(p.text))}
          >
            {p.label}
            <span>{p.says}</span>
          </button>
        ))}
      </div>
      <textarea
        className={q.box}
        value={text}
        rows={2}
        onChange={(e) => setText(e.target.value)}
        onKeyDown={(e) => {
          if (e.key === 'Enter' && !e.shiftKey) {
            e.preventDefault();
            run(text);
          } else if (e.key === 'Escape') clear();
        }}
        spellCheck={false}
        aria-label="SeqQL query"
      />
      <div className={f.row}>
        <TextAction onClick={() => run(text)}>Run</TextAction>
        <TextAction tone="quiet" disabled={!ran?.sel?.groupCount} onClick={to3d}>
          Select in 3D
        </TextAction>
        <TextAction tone="quiet" disabled={!ran} onClick={clear}>
          Clear
        </TextAction>
        {ran?.sel && <span className={q.line}>{SeqSelection.summary(ran.sel, { maxParts: 6 })}</span>}
        {ran?.reach && <span className={q.note}>{ran.reach}</span>}
      </div>
      {ran && ran.said.length > 0 && <pre className={q.issues}>{ran.said.join('\n')}</pre>}
    </section>
  );
}

/** The 3D canvas's color per toolbar theme: Prism's base of that theme (its theme.css), as the sequence canvas has. */
const CANVAS = { light: '#eeece7', dark: '#111318', blue: '#2d3e50' } as const;

export function MolstarDemo({ theme }: { theme: keyof typeof CANVAS }) {
  const [bytes, setBytes] = useState<Map<string, Uint8Array> | null>(null);
  useEffect(() => void structureBytes().then(setBytes), []);
  const host = useRef<HTMLDivElement>(null);
  const plugin = usePlugin(host, bytes);
  const [shown, setShown] = useState<Shown | null>(null);
  const shownRef = useRef<Shown | null>(null);
  shownRef.current = shown;
  const show = useRef((s: Shown) => setShown(s)).current;
  const ctx = useViewer(bytes, shownRef, show);
  const linked = useRef<Linked | null>(null);
  // The track menu has "Show in 3D", "Color 3D by this track" and "Apply to 3D": `clean` leaves it out, the story needs it.
  const ui = useStoryUi(null, { components: partSlots(['TrackActions']) });
  const [status, setStatus] = useState('Loading…');
  const background = useRef('#ffffff');
  background.current = CANVAS[theme] ?? CANVAS.light;

  useEffect(() => {
    if (plugin) setBackground(plugin, Number.parseInt((CANVAS[theme] ?? CANVAS.light).slice(1), 16));
  }, [plugin, theme]);

  // 1TUP in 3D at first, unless a structure was picked meanwhile.
  useEffect(() => {
    if (!ctx) return;
    const sub = ctx.graph.assembled.subscribe((a) => {
      if (a?.view) setShown((was) => was ?? { structure: refOf(P53_STRUCTURES[0]), colorBy: null });
    });
    return () => sub.unsubscribe();
  }, [ctx]);

  // The coverage entry of the structure in 3D: marked and open (its chain rows carry the per-row marks).
  const entry = useMemo(() => (ctx ? shownEntry(ctx) : null), [ctx]);
  useEffect(() => () => entry?.dispose(), [entry]);

  // What is shown: its coverage entry marked and opened first, then the 3D state composed now and loaded (the camera
  // kept when only its coloring changed); then linked.
  const last = useRef<string | null>(null);
  useEffect(() => {
    if (!plugin || !ctx || !shown) return;
    let unlink: (() => void) | null = null;
    let link: StructureLink | null = null;
    let live = true;
    const id = sourceId(shown.structure) ?? '';
    const same = last.current === id;
    last.current = id;
    void (async () => {
      if (!same) await entry?.show(P53_STRUCTURES.find((x) => x.source.id === id) ?? null);
      const r = live
        ? composeStructureView(ctx, shown.structure, {
            ...(shown.colorBy ? { tracks: [shown.colorBy], where: shown.where ?? null } : {}),
            defaults: shown.defaults ?? { representation: 'cartoon', base: '#c8c8c8' },
            background: background.current,
          })
        : null;
      if (!r) return;
      let structure: Awaited<ReturnType<typeof showMvs>> = null;
      try {
        structure = await showMvs(plugin, r.state, { keepCamera: same });
      } catch (e) {
        setStatus(`Mol* could not load it: ${e instanceof Error ? e.message : e}`);
        return;
      }
      const set = ctx.graph.assembled.value?.set;
      if (!live || !structure || !set) return;
      link = StructureLink.attach(ctx.systems, '3d', structure, set, shown.structure.complex ?? id.toLowerCase());
      const view3d = pluginView3D(plugin, structure);
      // The link's reactions, composed in its options: a Seq* hover of more than 1,000 items lights its whole chains.
      unlink = linkViews(ctx, view3d, link, { hover: { maxItems: 1_000, over: 'chains' } });
      linked.current = { plugin, link, view3d };
      setStatus(`${id.toUpperCase()} in 3D, linked: hover or select in either view.`);
    })();
    return () => {
      live = false;
      linked.current = null;
      unlink?.();
      link?.detach();
    };
  }, [plugin, ctx, shown, entry]);

  const current = sourceId(shown?.structure);
  return (
    <div className={f.stack}>
      <div className={m.switch}>
        {P53_STRUCTURES.map((s) => (
          <button
            key={s.complex}
            type="button"
            className={m.choice}
            aria-pressed={current === s.source.id}
            onClick={() => setShown({ structure: refOf(s), colorBy: null })}
          >
            {s.label}
            <span>{s.covers}</span>
          </button>
        ))}
      </div>
      <div className={f.status}>{status}</div>
      {ctx && <MolQueries ctx={ctx} linked={linked} />}
      <div className={m.split}>
        <div className={m.stage}>
          {ctx && shown && <ColorChooser ctx={ctx} shown={shown} show={show} />}
          <div ref={host} className={m.canvas} />
        </div>
        <section className={f.card}>
          {ctx ? (
            <SeqViewerProvider context={ctx} ui={ui}>
              <PrismCanvas>
                <SeqViewerLayout />
              </PrismCanvas>
            </SeqViewerProvider>
          ) : (
            <div className={f.status}>Loading…</div>
          )}
        </section>
      </div>
    </div>
  );
}
