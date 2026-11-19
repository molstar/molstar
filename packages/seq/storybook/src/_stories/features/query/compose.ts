/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type Issue, parseAddress, printAddress } from '@molstar/core/repr/tree';
import type { SeqPorts } from '@molstar/seqx-core/context';
import { positionCount, type QueryResultType } from '@molstar/seqx-core/query';
import type { ExampleView } from '@molstar/seqx-examples';
import { DefaultSeqViewerSpec } from '@molstar/seqx-viewer/default-spec';
import { mountToolbar } from '@molstar/seqx-viewer/mount-toolbar';
import { createSeqViewer, type SeqViewer } from '@molstar/seqx-viewer/viewer';

// The query stories' composition, DOM-free: the story and its headless twin call these. A typed address is parsed,
// checked and resolved by the library's own paths (`parseAddress`, `printAddress`, the viewer's `AddressModel`), and
// its result is selected or highlighted through the address model, never around it.

export interface ComposeOptions {
  /** Where the canvas mounts: none in the headless twin. */
  readonly target?: HTMLElement;
  /** The ports: the twin's reads fixtures from memory. */
  readonly ports?: SeqPorts;
}

/** A viewer of an example view, composed from the default spec with the view's registry. */
export function composeViewer(view: ExampleView, options: ComposeOptions = {}): SeqViewer {
  const base = DefaultSeqViewerSpec();
  const extra = view.registry?.() ?? [];
  return createSeqViewer(
    {
      ...base,
      ui: options.target ? { ...base.ui, toolbar: mountToolbar } : undefined,
      registry: [...(base.registry ?? []), ...extra],
      ...(options.ports ? { ports: options.ports } : {}),
    },
    options.target,
  );
}

/** What the result does in the view: selected (session state), or highlighted on a channel of its own. */
export type Emphasis = 'select' | 'highlight';

/** The highlight channel the query box uses, apart from the links' and the hover's. */
export const QueryChannel = 'query';

export type QueryStatus = 'empty' | 'resolved' | 'waiting' | 'missing' | 'failed';

/** What a typed query became: its canonical form, what it resolved to, and the issues with their hints. */
export interface QueryReport {
  readonly text: string;
  readonly status: QueryStatus;
  /** The canonical (printed) form: absent when the text does not parse. */
  readonly printed?: string;
  /** The result's kind: nodes, elements or a location. */
  readonly kind?: QueryResultType;
  /** What the result holds, as a phrase (`3 positions on 1 system`). */
  readonly size?: string;
  /** The result as positions (elements by their extents), per system. */
  readonly systems: readonly SystemRanges[];
  /** A Location's rects, each its two sides (rows by columns): none for other results. */
  readonly rects: readonly { readonly a: SystemRanges; readonly b: SystemRanges }[];
  readonly issues: readonly Issue[];
}

/** Half-open ranges of positions on a system. */
export interface SystemRanges {
  readonly system: string;
  readonly ranges: readonly (readonly [number, number])[];
}

const pairs = (flat: Int32Array): (readonly [number, number])[] =>
  Array.from({ length: flat.length / 2 }, (_, i) => [flat[2 * i], flat[2 * i + 1]] as const);

const plural = (n: number, one: string, many = `${one}s`) => `${n} ${n === 1 ? one : many}`;

/**
 * Parses, checks and resolves a query as typed, and selects or highlights its result in the view (what it named
 * before is cleared first, so the view never shows an older query). The library's own paths: `parseAddress`,
 * `printAddress`, `AddressModel.locate` and `AddressModel.select` or `highlightAddress`.
 */
export function runQuery(viewer: SeqViewer, text: string, emphasis: Emphasis): QueryReport {
  clearQuery(viewer);
  const trimmed = text.trim();
  if (!trimmed) return { text, status: 'empty', systems: [], rects: [], issues: [] };
  const parsed = parseAddress(trimmed);
  const printed = parsed.pattern ? printAddress(parsed.pattern.segments, parsed.pattern.subtree) : undefined;
  const { resolution, location } = viewer.addresses.locate(trimmed);
  const r = resolution.result;
  if (!r) {
    return {
      text,
      status: resolution.status === 'resolved' ? 'failed' : resolution.status,
      ...(printed !== undefined ? { printed } : {}),
      systems: [],
      rects: [],
      issues: resolution.issues,
    };
  }
  const systems = location.parts.map((p) => ({ system: p.system, ranges: pairs(p.ranges) }));
  const rects = (r.type === 'location' ? r.location.rects : []).map((x) => ({
    a: { system: x.a.system, ranges: pairs(x.a.ranges) },
    b: { system: x.b.system, ranges: pairs(x.b.ranges) },
  }));
  const size =
    r.type === 'nodes'
      ? plural(r.nodes.length, 'node')
      : r.type === 'elements'
        ? `${plural(r.elements.keys.length, 'element')} (${r.elements.keys.map((k) => k.key).join(', ')})`
        : [
            r.location.parts.length
              ? `${plural(positionCount(r.location), 'position')} on ${plural(r.location.parts.length, 'system')}`
              : '',
            r.location.rects.length ? plural(r.location.rects.length, 'rect') : '',
          ]
            .filter(Boolean)
            .join(', ');
  // the canonical form is what the view holds: the same query typed another way is the same address
  const address = printed ?? trimmed;
  if (emphasis === 'select') viewer.addresses.select(address);
  else viewer.addresses.highlightAddress(address, QueryChannel);
  return {
    text,
    status: 'resolved',
    ...(printed !== undefined ? { printed } : {}),
    kind: r.type,
    size,
    systems,
    rects,
    issues: resolution.issues,
  };
}

/** Clears what the query box set: its selection and its highlight channel. */
export function clearQuery(viewer: SeqViewer): void {
  viewer.selection.clear();
  viewer.addresses.clearHighlight(QueryChannel);
}

/** One worked example: its explanation, the queries it offers (the first is the one it opens with) and how it shows. */
export interface QueryExample {
  readonly id: string;
  readonly title: string;
  /** The view the example's data comes from (`FeatureViews` id). */
  readonly view: string;
  /** What it shows, a sentence or two. */
  readonly about: string;
  readonly queries: readonly { readonly label: string; readonly query: string }[];
  readonly emphasis: Emphasis;
}

export const QueryExamples: readonly QueryExample[] = [
  {
    id: 'category',
    title: 'Items by category',
    view: 'query',
    about:
      'An item form picks the items of an annotation by a field: its id, its label or its category. A list of values matches any of them, and the result is Elements (keys of what the data makes), drawn by their extents.',
    emphasis: 'select',
    queries: [
      { label: 'DNA-binding domain', query: '/data/domains/(item :label DNA-binding)' },
      { label: 'Two categories', query: '/data/domains/(item :category [binding oligomer])' },
      { label: 'By id', query: '/data/domains/(item :id [TAD CTD])' },
    ],
  },
  {
    id: 'ranges',
    title: 'Ranges and slices',
    view: 'query',
    about:
      'A range of positions on a system: labels are 1-based unless a numbering says otherwise, and `index` makes them 0-based slices. A single residue is `at`. The result is a Location.',
    emphasis: 'select',
    queries: [
      { label: 'Residues 94 to 292', query: '/data/p53/(range :from 94 :to 292)' },
      { label: 'Residue 175', query: '/data/p53/(at :label 175)' },
      { label: 'A 0-based slice', query: '/data/p53/(range :from 0 :numbering index :to 9)' },
      { label: 'On the gene', query: '/data/tp53/(range :from 101 :to 400)' },
    ],
  },
  {
    id: 'on-gene',
    title: 'A result on the gene',
    view: 'query',
    about:
      'A result on the protein is drawn on the gene that reads along it: the selection is projected through the codon map, three nucleotides a residue, split at the introns. The result itself stays on the protein.',
    emphasis: 'select',
    queries: [
      { label: 'The DNA-binding domain', query: '/data/domains/(item :id DBD)' },
      { label: 'Residues 94 to 292', query: '/data/p53/(range :from 94 :to 292)' },
      { label: 'Three hotspots', query: '/data/p53/(or :of [(at :label 175) (at :label 248) (at :label 273)])' },
    ],
  },
  {
    id: 'combined',
    title: 'And, or, not',
    view: 'query',
    about:
      'Combinators work on one result type: positions here. `and` intersects, `or` joins and `not` takes the complement within the system. Joining elements with positions is an issue, not a guess: convert with `extent` first.',
    emphasis: 'highlight',
    queries: [
      {
        label: 'The domain less a stretch',
        query: '/data/p53/(and :of [(range :from 94 :to 292) (not :of (range :from 150 :to 200))])',
      },
      { label: 'Two ranges', query: '/data/p53/(or :of [(range :from 1 :to 40) (range :from 320 :to 360)])' },
      { label: 'Where two overlap', query: '/data/p53/(and :of [(range :from 94 :to 292) (range :from 250 :to 360)])' },
      { label: 'Mixed types (an issue)', query: '/data/domains/(or :of [(item :id TET) (extent :of (item :id TET))])' },
    ],
  },
  {
    id: 'alignment',
    title: 'An alignment’s rows and columns',
    view: 'query-alignment',
    about:
      'An alignment has two axes: a row form names rows by their labels, a range names columns. A rect pairs them, rows by columns; the selection is a rectangle of cells.',
    emphasis: 'select',
    queries: [
      { label: 'A row', query: '/data/globins/(row :key "HBB_HUMAN/1-147")' },
      { label: 'Some columns', query: '/data/globins/(range :from 20 :to 60)' },
      {
        label: 'Rows by columns',
        query: '/data/globins/(rect :a (row :key ["HBB_HUMAN/1-147" "GLB1_GEN/1-137"]) :b (range :from 20 :to 60))',
      },
    ],
  },
];
