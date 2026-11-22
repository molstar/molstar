/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import type { SeqRegistryEntry } from '@molstar/seqx-core/context';
import {
  delayedSource,
  type FormatRuntime,
  openSeq,
  type OpenedSource,
  type RegionRequest,
  type SeqData,
  SourceOpenError,
  type StoreClock,
} from '@molstar/seqx-core/data';
import { assembly, createSeqSpec, genome, type ViewBuilder } from '@molstar/seqx-core/lang';
import type { ExampleView } from '../views.js';
import { BinsPitch, Grape, grapePlum, Peach, peachGrape, peachSelf, Plum, type ToyGenome } from './genomes.js';

// The views across sections (Features/Across sections; across-sections, Defaults): genomes stacked, a bridge between
// neighbours drawing an alignment's blocks as ribbons. Each alignment a source generated at open and answered after a
// delay on the store's clock, a lower genome's chromosomes twice the upper's, so a band fills a side at a time. A source
// names the format `seq-across` and its alignment by `uri` (`generated:peach-grape`).

/** The format the views' sources name. */
export const AcrossFormat = 'seq-across';

const Alignments: Readonly<Record<string, () => SeqData>> = {
  'generated:peach-grape': () => peachGrape(),
  'generated:peach-grape-dense': () => peachGrape(BinsPitch),
  'generated:grape-plum': grapePlum,
  'generated:peach-self': peachSelf,
};

/** The chromosomes read late: a lower genome's (grape's below peach, plum's below grape). */
const Late: ReadonlySet<string> = new Set([...Grape.chromosomes, ...Plum.chromosomes].map((c) => c.name));

export interface AcrossOptions {
  /** The clock the answers wait on: real timers by default. */
  readonly clock?: StoreClock;
  /** A read's delay, ms (a lower genome's twice it: a band fills from its upper side); by request when a function. */
  readonly delayMs?: number | ((request: RegionRequest) => number);
}

/** A lower genome's chromosomes wait twice the upper's delay, so the tiles arrive a side at a time. */
const bySide = (delayMs: number) => (r: RegionRequest) => (Late.has(r.system) ? 2 * delayMs : delayMs);

/** The runtime of `seq-across` sources: `uri` names an alignment, opened as a `seq` source answered after a delay. */
export function acrossRuntime(options: AcrossOptions = {}): FormatRuntime {
  return {
    name: 'seq-across',
    formats: [AcrossFormat],
    async open(source): Promise<OpenedSource> {
      const make = source.uri === undefined ? undefined : Alignments[source.uri];
      if (!make) throw new SourceOpenError([]);
      const delay = options.delayMs ?? 300;
      return delayedSource(openSeq(make()), {
        ...(options.clock ? { clock: options.clock } : {}),
        delayMs: typeof delay === 'number' ? bySide(delay) : delay,
      });
    },
  };
}

/** What a view needs besides the default registry: the alignments' runtime. */
export function registry(options: { readonly delayMs?: number } = {}): SeqRegistryEntry[] {
  return [{ kind: 'format-runtime', runtime: acrossRuntime(options) }];
}

/** The sections' addresses (what the tests read). */
export const PeachSection = '/views/main/peach';
export const GrapeSection = '/views/main/grape';
export const PlumSection = '/views/main/plum';
/** The lower peach of one genome against itself. */
export const PeachBelowSection = '/views/main/peach-below';

type Builder = ReturnType<typeof createSeqSpec>;

const asm = (b: Builder, id: string, g: ToyGenome) => b.data.assembly(id, assembly(g.name, g.chromosomes));
const alignment = (b: Builder, id: string, uri: string) => b.data.source(id, { format: AcrossFormat, uri } as never);
/** A bridge between two sections drawing a pairs node (its `legend` none: the sides' hues are the chromosomes'). */
const bridge = (
  v: ViewBuilder,
  id: string,
  a: string,
  b: string,
  data: unknown,
  params: Record<string, unknown> = {},
  layer: Record<string, unknown> = {},
) =>
  v.bridge(id, { a, b, ...params } as never, (br) =>
    br.layer(id, { type: 'pairs', data, legend: 'none', ...layer } as never),
  );

interface TwoOptions {
  readonly title: string;
  readonly uri?: string;
  /** Grape's `chr4` cut to its first 300 kb (a section of blocks, one with a `range`): the stub's view. */
  readonly cut?: boolean;
  readonly hidden?: boolean;
}

/** Peach above, grape below, the bridge between drawing their alignment. */
function peachOverGrape(o: TwoOptions): Node {
  const b = createSeqSpec({ title: o.title });
  const peach = asm(b, 'peach', Peach);
  const grape = asm(b, 'grape', Grape);
  const src = alignment(b, 'alignment', o.uri ?? 'generated:peach-grape');
  const aln = b.data.pairs('aln', { a: peach, b: grape, source: src });
  b.views.view('main', (v) => {
    genome(v, 'peach', peach);
    bridge(v, 'peach-grape', PeachSection, GrapeSection, aln, o.hidden ? { hidden: true } : {});
    if (!o.cut) genome(v, 'grape', grape);
    else
      v.section('grape', {}, (s) => {
        s.track('ruler').layer('ruler', { type: 'ruler', units: 'bases' } as never);
        for (const c of Grape.chromosomes)
          s.block(c.name, {
            of: `/data/grape/(sequence :name ${c.name})`,
            ...(c.name === 'chr4' ? { range: [0, CutChr4] } : {}),
          } as never);
      });
  });
  return b.build();
}

/**
 * One chromosome each, read a side at a time (each section its own sequence: a read a side, grape's after peach's):
 * the pending hatch along the edge of the side not loaded, then the ribbons.
 */
function pendingView(): Node {
  const b = createSeqSpec({ title: 'Peach and grape, arriving' });
  const src = alignment(b, 'alignment', 'generated:peach-grape');
  const pp = b.data.sequence('Pp02', { alphabet: 'dna', source: src, name: 'Pp02' });
  const g2 = b.data.sequence('chr2', { alphabet: 'dna', source: src, name: 'chr2' });
  const aln = b.data.pairs('aln', { a: pp, b: g2, source: src } as never);
  b.views.view('main', (v) => {
    v.section('peach', { of: pp })
      .track('ruler')
      .layer('ruler', { type: 'ruler', units: 'bases' } as never);
    bridge(v, 'peach-grape', PeachSection, GrapeSection, aln);
    v.section('grape', { of: g2 })
      .track('ruler')
      .layer('ruler', { type: 'ruler', units: 'bases' } as never);
  });
  return b.build();
}

/** Peach against itself: two sections of one assembly, the bridge drawing its duplications both ways (`mirror`). */
function peachSelfView(): Node {
  const b = createSeqSpec({ title: 'Peach against itself' });
  const peach = asm(b, 'peach', Peach);
  const dups = b.data.pairs('dups', { a: peach, b: peach, source: alignment(b, 'alignment', 'generated:peach-self') });
  b.views.view('main', (v) => {
    genome(v, 'peach', peach);
    bridge(v, 'self', PeachSection, PeachBelowSection, dups, {}, { mirror: true });
    // equal block lists: the lower section follows the upper (one genome, moved together)
    v.section('peach-below', { of: peach, follow: PeachSection } as never, (s) =>
      s.track('ruler').layer('ruler', { type: 'ruler', units: 'bases' } as never),
    );
  });
  return b.build();
}

/** Three genomes chained: peach, grape, plum, a bridge between each two. */
function chainView(): Node {
  const b = createSeqSpec({ title: 'Three genomes chained' });
  const peach = asm(b, 'peach', Peach);
  const grape = asm(b, 'grape', Grape);
  const plum = asm(b, 'plum', Plum);
  const pg = b.data.pairs('pg', { a: peach, b: grape, source: alignment(b, 'pg-alignment', 'generated:peach-grape') });
  const gp = b.data.pairs('gp', { a: grape, b: plum, source: alignment(b, 'gp-alignment', 'generated:grape-plum') });
  b.views.view('main', (v) => {
    genome(v, 'peach', peach);
    bridge(v, 'peach-grape', PeachSection, GrapeSection, pg);
    genome(v, 'grape', grape);
    bridge(v, 'grape-plum', GrapeSection, PlumSection, gp);
    genome(v, 'plum', plum);
  });
  return b.build();
}

const view = (id: string, title: string, summary: string, build: () => Node): ExampleView => ({
  dataset: 'across',
  id,
  title,
  summary,
  build,
  registry,
});

/** Where the stub view cuts grape's `chr4`. */
export const CutChr4 = 300_000;

/** Where the zoomed-in view opens: peach's `Pp02` from 1.2 to 1.5 Mb (Defaults, 3). */
export const ZoomedLocus = 'Pp02:1,200,001-1,500,000';

export const AcrossViews: readonly ExampleView[] = [
  view(
    'across-fit',
    'Two genomes at the fit',
    'Peach above, grape below, a 96 px band between: each alignment block a ribbon in its peach chromosome’s hue; Pp02’s inversion crossed, Pp03’s tail carried to grape’s chr4 (a translocation), Pp01’s lost stretch no ribbon (a gap).',
    () => peachOverGrape({ title: 'Peach and grape' }),
  ),
  view(
    'across-bins',
    'Zoomed out: bins',
    'The same genomes aligned in about 6,000 blocks: past the 4,096 a band draws, bin pairs of each side’s width, their opacity by the bases aligned; zoomed in past the gate, the blocks replace them.',
    () => peachOverGrape({ title: 'Peach and grape, dense', uri: 'generated:peach-grape-dense' }),
  ),
  view(
    'across-zoomed',
    'One side zoomed in',
    'Peach zoomed into Pp02 from 1.2 to 1.5 Mb, grape at the fit: the zoomed side’s blocks fan to their true x on grape, the ribbons narrowing to grape’s 300 kb; ends past a side’s edge leave the band sideways.',
    () => peachOverGrape({ title: 'Peach zoomed in' }),
  ),
  view(
    'across-stub',
    'A stub',
    'Grape’s chr4 cut to its first 300 kb (a block’s `range`): the ribbons of peach Pp03’s tail landing past it have no place below, each a quarter-height taper from its peach end.',
    () => peachOverGrape({ title: 'Grape’s chr4 cut', cut: true }),
  ),
  view(
    'across-pending',
    'Pending, a side at a time',
    'One chromosome each, grape’s reads after peach’s: the pending hatch along the edge of the side not loaded, then the ribbons.',
    pendingView,
  ),
  view(
    'across-hidden',
    'A reader without bridges',
    'The bridge hidden (`hidden`, kept in the spec): the two sections stacked with the sections’ gap, the band gone.',
    () => peachOverGrape({ title: 'Peach and grape, no band', hidden: true }),
  ),
  view(
    'across-self',
    'One genome against itself',
    'Peach above and below, its duplications drawn both ways (`mirror`): one within Pp02, one onto Pp03, one reversed onto Pp04; the lower section following the upper, so they move together.',
    peachSelfView,
  ),
  view(
    'across-chain',
    'Three genomes chained',
    'Peach, grape, plum: a bridge between each two; grape’s chr1 tail carried reversed onto plum’s Pd4, its chr3 and chr4 fused in plum’s Pd3.',
    chainView,
  ),
];
