/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The genome-scale story: the synthetic genome of genome-fixture.ts on the real
// viewer, timed. The renderers are the kernel's, each prepare and draw wrapped by a timer (the viewer's
// _test/canvas/genome-probe.ts: allocation-free, sums in a typed array), read back as `__genome` in the console; load phases
// shown under the viewer.

import {
  compile,
  createRegistries,
  KernelRegistry,
  type Loaders,
  register,
  type SeqViewSpecDoc,
  validate,
} from '@molstar/seq-core/format';
import {
  createSeqViewer,
  DefaultRegistry,
  DefaultSeqViewerSpec,
  SeqViewerConfig,
  type SeqViewerConfigItem,
} from '@molstar/seq-viewer';
import { windowScheduler } from '@molstar/seq-viewer/canvas';
import { SeqViewerLayout, SeqViewerProvider } from '@molstar/seq-viewer/react';
import { type CSSProperties, useEffect, useMemo, useState } from 'react';
import { genomeDoc, synthGenome } from '../../../../core/_test/format/genome-fixture.js';
import SOURCE from '../../../../core/_test/format/genome-fixture.ts?raw';
import { PROBE, TIMED } from '../../../../viewer/_test/canvas/genome-probe.js';
import { PrismPage } from '../support/seqstar-prism.js';
import { Choices, TextAction } from '../ui/index.js';
import f from './seqstar-format.module.css';
import { PrismCanvas, prism } from './seqstar-prism.js';
import {
  BuilderCode,
  type Mounted,
  type NavigatorChoice,
  presets,
  STORY_CHROME,
  type UiChoice,
  useStoryUi,
} from './seqstar-viewer.js';

const SECTION = 'genome.main';
(globalThis as { __genome?: typeof PROBE }).__genome = PROBE;
/** Everything is by value: a load is a bug, never a fetch. */
const NO_LOADS: Loaders = { load: () => Promise.reject(new Error('the genome story loads nothing')) };

export interface GenomeChoice {
  readonly proteins: number;
  readonly gcTable: 'one' | 'each';
}
interface Built {
  readonly doc: SeqViewSpecDoc;
  readonly synthMs: number;
  readonly builderMs: number;
}
const BUILT = new Map<string, Built>();
function built(c: GenomeChoice): Built {
  const key = `${c.proteins}/${c.gcTable}`;
  let b = BUILT.get(key);
  if (!b) {
    const t0 = performance.now();
    const g = synthGenome({ proteins: c.proteins });
    const t1 = performance.now();
    const doc = genomeDoc(g, { gcTable: c.gcTable });
    b = { doc, synthMs: t1 - t0, builderMs: performance.now() - t1 };
    BUILT.set(key, b);
  }
  PROBE.doc = b.doc;
  return b;
}

interface Load {
  graphMs: number;
  firstFrameMs: number;
  /** The first frame's main-thread ms (the surface's `lastMs`), and a full redraw's with its raster. */
  frameMs: number;
  redrawMs: number;
  /** JS heap before `setDocument` and after the first frame (collected first when `gc` is exposed: Playwright). */
  heapBeforeMb: number;
  heapMb: number;
  prepare: Record<string, number>;
  issues: number;
}
const heapMb = () => {
  (globalThis as { gc?: () => void }).gc?.();
  const m = (performance as { memory?: { usedJSHeapSize: number } }).memory;
  return m ? Math.round(m.usedJSHeapSize / 1e6) : Number.NaN;
};

/** A context over the genome with the timed renderers; its load timed from `setDocument` to the first whole frame. */
function useGenome(doc: SeqViewSpecDoc, navigator: NavigatorChoice) {
  const [m, setM] = useState<Mounted | null>(null);
  const [load, setLoad] = useState<Load | null>(null);
  useEffect(() => {
    const spec = DefaultSeqViewerSpec();
    const registry = [{ ...DefaultRegistry, presets: presets({ navigator }), renderers: TIMED }];
    const config: [SeqViewerConfigItem, unknown][] = [[SeqViewerConfig.UI.Chrome, STORY_CHROME]];
    const ctx = createSeqViewer({ ...spec, registry, config }, { loaders: NO_LOADS, scheduler: windowScheduler() });
    const next: Mounted = { ctx, surfaces: new Map() };
    const all = (globalThis as { __seqViewers?: Set<Mounted> }).__seqViewers;
    all?.add(next);
    for (const k of Object.keys(PROBE.prepare)) delete PROBE.prepare[k];
    setLoad(null);
    const heapBeforeMb = heapMb();
    const t0 = performance.now();
    let raf = 0;
    void ctx.setDocument(doc).then((r) => {
      const graphMs = performance.now() - t0;
      // The first whole frame: the surface mounted and drawn (a layer is never presented half drawn).
      const wait = () => {
        const surface = next.surfaces.get(SECTION);
        if (ctx.disposed) return;
        if (!surface || surface.stats.frames === 0) {
          raf = requestAnimationFrame(wait);
          return;
        }
        const firstFrameMs = performance.now() - t0;
        const frameMs = surface.stats.lastMs;
        // Then one full redraw timed, as the browser test does (the deferred raster flushed by a one-pixel read).
        const t1 = performance.now();
        surface.composer.redraw();
        surface.invalidate(7);
        surface.renderNow();
        surface.element.querySelector('canvas')?.getContext('2d')?.getImageData(0, 0, 1, 1);
        setLoad({
          graphMs,
          firstFrameMs,
          frameMs,
          redrawMs: performance.now() - t1,
          heapBeforeMb,
          heapMb: heapMb(),
          prepare: { ...PROBE.prepare },
          issues: r.issues.length,
        });
      };
      wait();
    });
    setM(next);
    return () => {
      cancelAnimationFrame(raf);
      all?.delete(next);
      setM(null);
      ctx.dispose();
    };
  }, [doc, navigator]);
  return { m, load };
}

/** validate and compile on their own (the context runs them per node inside its graph): the format's share. */
async function phases(doc: SeqViewSpecDoc) {
  const reg = createRegistries();
  register(reg, KernelRegistry);
  const t0 = performance.now();
  const issues = validate(doc, reg);
  const t1 = performance.now();
  const c = await compile(doc, reg, NO_LOADS);
  const t2 = performance.now();
  return { validateMs: t1 - t0, compileMs: t2 - t1, issues: issues.length + c.issues.length };
}

const ms = (x: number) => `${Math.round(x).toLocaleString('en-US')} ms`;
const PROTEINS = [0, 30, 300] as const;

export function GenomeView({ navigator, ui }: { navigator: NavigatorChoice; ui: UiChoice }) {
  const [choice, setChoice] = useState<GenomeChoice>({ proteins: 300, gcTable: 'each' });
  const b = useMemo(() => built(choice), [choice]);
  const { m, load } = useGenome(b.doc, navigator);
  const spec = useStoryUi(m);
  const [timed, setTimed] = useState<string>('');
  const style = { '--seq-section-max-height': '720px' } as CSSProperties;
  const dpr = typeof window === 'undefined' ? 1 : window.devicePixelRatio;
  return (
    <PrismPage
      title="A genome at scale: 3,000,000 nucleotides on one axis"
      caption={
        <>
          Synthetic and seeded, made in the page: 24 chromosomes in the human proportions as blocks, cytobands, 3,000
          genes as transcripts (exons as fragments of one item, UTR and CDS), {choice.proteins} proteins translated from
          their CDS with codon mappings, GC content per position (1 kb windows), 20,000 variants, 10,000 repeats.
        </>
      }
    >
      <BuilderCode
        source={{
          code: SOURCE.slice(SOURCE.indexOf('// #region genome\n') + 18, SOURCE.indexOf('// #endregion')).trimEnd(),
          uses: 'The data from synthGenome (core/_test/format/genome-fixture.ts): a seeded PRNG, no network',
        }}
      />
      <div className={prism.controls}>
        <div className={prism.narrow}>
          <Choices
            label="Proteins"
            choices={PROTEINS.map((n) => ({ id: String(n), label: String(n) }))}
            active={String(choice.proteins)}
            onPick={(id) =>
              setChoice({ ...choice, proteins: PROTEINS.find((n) => String(n) === id) ?? choice.proteins })
            }
          />
        </div>
        <div className={prism.wide}>
          <Choices
            label="GC"
            choices={[
              { id: 'each', label: 'a table per chromosome' },
              { id: 'one', label: 'one table, a target per row' },
            ]}
            active={choice.gcTable}
            onPick={(id) => setChoice({ ...choice, gcTable: id as 'each' | 'one' })}
          />
        </div>
      </div>
      <section className={f.card} style={style}>
        <h2 className={f.heading}>Genome</h2>
        {m ? (
          <SeqViewerProvider context={m.ctx} ui={spec}>
            <PrismCanvas>
              <SeqViewerLayout />
            </PrismCanvas>
          </SeqViewerProvider>
        ) : (
          <div className={f.status}>Loading…</div>
        )}
        {ui === 'DOM' && <div className={f.code}>This story draws the React UI only.</div>}
      </section>
      <div className={f.code}>
        {`Window ${window.innerWidth} × ${window.innerHeight} CSS px, DPR ${dpr}. Synthesis ${ms(b.synthMs)}, builder ${ms(b.builderMs)}. `}
        {load
          ? `Graph (load, resolve, prepare) ${ms(load.graphMs)}; first whole frame ${ms(load.firstFrameMs)} after setDocument (that frame ${ms(load.frameMs)}, a full redraw ${ms(load.redrawMs)}); JS heap ${load.heapBeforeMb} MB before, ${load.heapMb} MB after; ${load.issues} issues. Prepare: ${Object.entries(
              load.prepare,
            )
              .map(([k, v]) => `${k} ${ms(v)}`)
              .join(', ')}.`
          : 'Loading…'}
      </div>
      <div className={f.row}>
        <TextAction tone="quiet" onClick={() => void phases(b.doc).then((r) => setTimed(JSON.stringify(r)))}>
          Time validate and compile alone
        </TextAction>
        <span className={f.code}>{timed}</span>
      </div>
    </PrismPage>
  );
}
