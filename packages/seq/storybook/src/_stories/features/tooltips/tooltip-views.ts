/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import { assembly, createSeqSpec, genome } from '@molstar/seqx-core/lang';
import { type ExampleView, GenomeFormatsViews, P53Views } from '@molstar/seqx-examples';
import { toyKaryotype } from '../genome/views.js';

// The Tooltips stories' views, DOM-free (the stories and their headless twin read them): the spec's words on rendered
// marks. By data, a GFF's transcripts speak their own attributes; by range, a curator's notes on p53 (inline, and from
// a source); over an assembly, notes by their sequences' names; Markdown in the template with an in-view link; and an
// export whose items are titled by their cards.

/** A spec with nodes added: `data` children and `views` children (the tooltips). */
function withNodes(spec: Node, add: { readonly data?: readonly Node[]; readonly views?: readonly Node[] }): Node {
  return {
    ...spec,
    children: (spec.children ?? []).map((c) =>
      c.kind === 'data' && add.data
        ? { ...c, children: [...(c.children ?? []), ...add.data] }
        : c.kind === 'views' && add.views
          ? { ...c, children: [...(c.children ?? []), ...add.views] }
          : c,
    ),
  };
}

const tooltip = (id: string, params: Record<string, unknown>): Node => ({
  kind: 'tooltip',
  id,
  params: params as never,
});

const p53Tracks = P53Views.find((v) => v.id === 'tracks')!;
const gff = GenomeFormatsViews.find((v) => v.id === 'gff-tp53')!;

/** p53's curated notes: a hotspot residue each, and the core domain's own (the protein's positions, 0-based). */
export const Notes = {
  start: [174, 247, 272, 93],
  end: [175, 248, 273, 292],
  tooltip: ['R175H: conformational', 'R248Q: contact', 'R273H: contact', 'the DNA-binding core'],
} as const;

/** Where each story's twin rests: a section's track and a position on its system, and the row it expects there. */
export interface TooltipProbe {
  /** Its track: by index, or by label. */
  readonly track: number | string;
  readonly position: number;
  /** The item hovered (by data), its owner and key. */
  readonly element?: { readonly owner: string; readonly key: string | number } | 'first';
  readonly label: string;
  readonly value: RegExp;
}

export interface TooltipView extends ExampleView {
  readonly hint: string;
  readonly probe: TooltipProbe;
  /** The story shows the SVG export beside the canvas (its items titled by their cards). */
  readonly svg?: boolean;
}

export const TooltipViews: readonly TooltipView[] = [
  {
    ...gff,
    dataset: 'tooltips',
    id: 'by-data',
    title: 'By data: a GFF’s attributes',
    summary: 'TP53’s transcripts from Ensembl’s GFF3, each speaking its own attributes (`on` the transcripts).',
    hint: 'Rest the pointer on a transcript: its card ends with a Transcript row, its name and its gene’s, from the file’s attributes (a tooltip `on` the data node). An empty lane says nothing.',
    build: () =>
      withNodes(gff.build(), {
        views: [
          tooltip('transcripts', { on: '/data/transcripts', label: 'Transcript', text: '{label} of *{gene_name}*' }),
        ],
      }),
    probe: { track: 2, position: 7_676_000, element: 'first', label: 'Transcript', value: / of TP53$/ },
  },
  {
    ...p53Tracks,
    dataset: 'tooltips',
    id: 'by-range-inline',
    title: 'By range: a curator’s notes',
    summary:
      'A table of ranges over p53 inline in the tooltip: each note on every mark it covers (a letter, a domain, a bar).',
    hint: 'Rest on residue 175, 248 or 273 (the letters, the DNA-binding domain or the conservation bars): a Note row says the hotspot. Anywhere in 94–292 the core’s own note shows.',
    build: () =>
      withNodes(p53Tracks.build(), {
        views: [tooltip('notes', { target: '/data/p53', rows: Notes, label: 'Note' })],
      }),
    probe: { track: 0, position: 174, label: 'Note', value: /^R175H/ },
  },
  {
    ...p53Tracks,
    dataset: 'tooltips',
    id: 'by-range-source',
    title: 'By range: notes from a source',
    summary: 'The same notes as JSON records read through a `source` under `data`, the store’s like any data.',
    hint: 'The notes from a source (JSON records, read through the store and wanted beside the view’s marks): rest on residue 248.',
    build: () =>
      withNodes(p53Tracks.build(), {
        data: [
          {
            kind: 'source',
            id: 'notes-file',
            params: {
              format: 'json',
              system: 'p53',
              value: Notes.start.map((s, i) => ({ start: s, end: Notes.end[i], note: Notes.tooltip[i] })) as never,
            },
          },
        ],
        views: [tooltip('notes', { target: '/data/p53', source: '/data/notes-file', field: 'note', label: 'Note' })],
      }),
    probe: { track: 0, position: 247, label: 'Note', value: /^R248Q/ },
  },
  {
    dataset: 'tooltips',
    id: 'over-assembly',
    title: 'Over an assembly by name',
    summary: 'Notes over a toy genome, a row naming its chromosome (`name`), on the features they cover.',
    hint: 'Rest on a feature of chr1 or chr2: a Region row names the note its range holds.',
    build: overAssembly,
    probe: { track: 'Features', position: 1_500_000, element: 'first', label: 'Region', value: /^first arm/ },
  },
  {
    ...p53Tracks,
    dataset: 'tooltips',
    id: 'markdown-link',
    title: 'Markdown with an in-view link',
    summary: 'The template’s Markdown (strong, code, a `!seq-focus` link to a ref); the values text.',
    hint: 'Rest on residue 175: the card’s link is entered across the gap and focuses the hotspots (a link acts through commands only).',
    build: () => {
      const spec = withNodes(p53Tracks.build(), {
        views: [
          tooltip('notes', {
            target: '/data/p53',
            rows: Notes,
            label: 'Note',
            text: '**{}** · [all hotspots](!seq-focus=/refs/hotspots)',
          }),
        ],
      });
      return {
        ...spec,
        children: [
          ...(spec.children ?? []),
          {
            kind: 'refs',
            id: 'refs',
            children: [{ kind: 'ref', id: 'hotspots', params: { to: '/data/p53/(range :from 175 :to 273)' } }],
          },
        ],
      };
    },
    probe: { track: 0, position: 174, label: 'Note', value: /^R175H: conformational · all hotspots$/ },
  },
  {
    ...p53Tracks,
    dataset: 'tooltips',
    id: 'svg-titles',
    title: 'The SVG’s titles',
    summary:
      'An export whose items are titled by their cards (`export.titles: card`): the browser’s tooltip on the SVG.',
    hint: 'Below the canvas, the SVG export: rest on a domain there, the browser’s tooltip is the domain’s card (its name, the Pfam row), the export’s `<title>`.',
    svg: true,
    build: () => {
      const spec = withNodes(p53Tracks.build(), {
        views: [tooltip('domains', { on: '/data/domains', label: 'Domain', text: '{label} ({id})' })],
      });
      return {
        ...spec,
        children: (spec.children ?? []).map((c) =>
          c.kind === 'views'
            ? {
                ...c,
                children: (c.children ?? []).map((v) =>
                  v.kind === 'view' ? { ...v, params: { ...v.params, export: { titles: 'card' } } } : v,
                ),
              }
            : c,
        ),
      };
    },
    probe: {
      track: 1,
      position: 150,
      element: { owner: '/data/domains', key: 'DBD' },
      label: 'Domain',
      value: /^DNA-binding \(DBD\)$/,
    },
  },
];

/** A toy genome's first chromosomes with features over it by name, and notes by name over them. */
function overAssembly(): Node {
  const b = createSeqSpec({ title: 'Notes over an assembly' });
  const g = b.data.assembly('toy', assembly('Toy', toyKaryotype().slice(0, 3)));
  const features = b.data.ranges('features', {
    target: g,
    rows: {
      target: ['chr1', 'chr1', 'chr2'],
      start: [1_000_000, 40_000_000, 2_000_000],
      end: [3_000_000, 42_000_000, 9_000_000],
      label: ['f1', 'f2', 'f3'],
    } as never,
  });
  b.views.tooltip('regions', {
    target: g,
    rows: {
      target: ['chr1', 'chr2'],
      start: [0, 1_000_000],
      end: [20_000_000, 5_000_000],
      tooltip: ['first arm of chr1', 'a stretch of chr2'],
    } as never,
    label: 'Region',
  });
  b.views.view('main', (v) =>
    genome(v, 'genome', g, {}, (s) =>
      s.track('features', { label: 'Features' }).layer('features', { type: 'ranges', data: features }),
    ),
  );
  return b.build();
}
