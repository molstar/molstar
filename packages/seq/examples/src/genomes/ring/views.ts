/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import type { SeqRegistryEntry } from '@molstar/seqx-core/context';
import { delayedRuntime, type StoreClock } from '@molstar/seqx-core/data';
import { assembly, createSeqSpec } from '@molstar/seqx-core/lang';
import { BedFormats, BedRuntime } from '@molstar/seqx-core/runtimes/bed';
import { rondo, rondoRing } from '@molstar/seqx-ext-rondo';
import { withRecordedPages } from '../../recorded-pages.js';
import type { ExampleView } from '../../views.js';
import { CytobandsFile, RecordedFiles } from '../files.js';
import { Grch38Chromosomes, RingInclude, ringBedpe, ringLinks, ringRows, syntheticBands } from './data.js';

// The genome ring's views (the showcase of `@molstar/seqx-ext-rondo`, Features/Genome ring): GRCh38's 23 chromosomes,
// 12,000 generated links between them as BEDPE (or held), UCSC's recorded cytobands (or synthetic ones) and the hue
// ring, written by `rondoRing`; one chromosome and a pair by the section's `include`, the square beside.

/** The ring section's address (what the tests and stories read). */
export const RingSection = '/views/main/ring';

export interface RingSpecOptions {
  /** The section's `include` (a chromosome, a pair); the 23 when absent. */
  readonly include?: readonly string[];
  /** The links held inline (rows over the assembly) rather than read from BEDPE. */
  readonly held?: boolean;
  /** Synthetic cytobands (no recording read). */
  readonly synthetic?: boolean;
  /** The links (12,000 by default). */
  readonly links?: number;
  /** The 24 hues (not colour-blind safe) rather than the safe categories. */
  readonly hues?: boolean;
  readonly square?: boolean;
  /** The links read from a file that is not there (the failed default: the placeholder and its issue). */
  readonly failed?: boolean;
  /** Written without the extension listed (a reader of it draws the same; what the fallback story shows). */
  readonly bare?: boolean;
  readonly title?: string;
}

/** A BEDPE the views name that is not committed (the failed default's read). */
export const MissingLinksFile = 'fixtures/genomes/ring-links-missing.bedpe';

/** The genome ring over GRCh38 (`rondoRing`), its links read from BEDPE or held. */
export function ringSpec(o: RingSpecOptions = {}): Node {
  const b = createSeqSpec({ title: o.title ?? 'The genome ring', ...(o.bare ? {} : { extensions: [rondo] }) });
  const g = b.data.assembly('grch38', assembly('GRCh38', Grch38Chromosomes, { include: [...RingInclude] } as never));
  const links = ringLinks(o.links ?? 12_000);
  const pairs = o.held
    ? b.data.pairs('links', { a: g, rows: ringRows(links) } as never)
    : b.data.pairs('links', {
        a: g,
        source: b.data.source(
          'links-bedpe',
          o.failed ? { format: 'bedpe', uri: MissingLinksFile } : { format: 'bedpe', value: ringBedpe(links) as never },
        ),
      } as never);
  const bandsSource = o.synthetic
    ? b.data.source('bands-bed', {
        format: 'bed',
        value: syntheticBands(Grch38Chromosomes) as never,
        fields: { stain: '4' },
      } as never)
    : b.data.source('bands-bed', {
        format: 'bed',
        uri: CytobandsFile,
        fields: { stain: '4' },
        attribution: 'UCSC Genome Browser, hg38 cytoBandIdeo',
        licence: 'LicenseRef-UCSC-freely-usable',
      } as never);
  const bands = b.data.ranges('bands', { target: g, source: bandsSource } as never);
  b.views.view('main', (v) => {
    rondoRing(v, {
      assembly: g,
      links: pairs,
      bands,
      opens: true,
      ...(o.hues !== undefined ? { hues: o.hues } : { hues: true }),
      ...(o.square !== undefined ? { square: o.square } : {}),
    });
  });
  let root = b.build();
  if (o.include) root = withInclude(root, o.include);
  return root;
}

/** The spec with the ring section's `include` set (as the `section.include` command patches it). */
export function withInclude(root: Node, include: readonly string[]): Node {
  const map = (n: Node, path: string): Node => {
    const here = n.id ? `${path}/${n.id}` : path;
    if (n.kind === 'section' && here === RingSection) return { ...n, params: { ...n.params, include: [...include] } };
    return n.children ? { ...n, children: n.children.map((c) => map(c, here)) } : n;
  };
  return map(root, '');
}

/**
 * What the views need besides the default registry: the `bed` runtime (BEDPE, the cytobands), the cytobands read from
 * their recording, each answer after `delayMs` (150) on `clock`.
 */
export function registry(options: { readonly delayMs?: number; readonly clock?: StoreClock } = {}): SeqRegistryEntry[] {
  const runtime = delayedRuntime(withRecordedPages(BedRuntime, RecordedFiles), {
    name: 'bed-ring',
    formats: BedFormats,
    delayMs: options.delayMs ?? 150,
    ...(options.clock ? { clock: options.clock } : {}),
  });
  return [{ kind: 'format-runtime', runtime }];
}

/** The genome ring's views (Features/Genome ring): their own list, beside the genome formats'. */
export const RingViews: readonly ExampleView[] = [
  {
    dataset: 'genomes',
    id: 'ring',
    title: 'The genome ring',
    summary:
      'GRCh38’s 23 chromosomes on a ring with UCSC’s cytobands and a hue ring, 12,000 generated links (BEDPE) bundled into ribbons between stems, wide by count; the square beside picks a chromosome or a pair. Hover a stem for its links; double-click a ribbon to open its pair.',
    build: () => ringSpec(),
    registry,
  },
  {
    dataset: 'genomes',
    id: 'ring-chromosome',
    title: 'One chromosome',
    summary:
      'chr4 alone on the ring (`include`): its own links as ribbons, items once they fit; the rest as stubs inward.',
    build: () => ringSpec({ include: ['chr4'], title: 'chr4' }),
    registry,
  },
  {
    dataset: 'genomes',
    id: 'ring-pair',
    title: 'A pair',
    summary: 'chr4 and chr7: their links between them as ribbons, the rest stubs in their partners’ hues.',
    build: () => ringSpec({ include: ['chr4', 'chr7'], title: 'chr4 × chr7' }),
    registry,
  },
];
