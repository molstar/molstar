/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import { type DynamicNode, type SeqData, seqToValue } from '@molstar/seqx-core/data';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { itemAddress, nodeAddress } from '@molstar/seqx-core/resolve';

// A gene read by region whose answers bring transcripts as dynamic nodes: each a coding ranges node (its exons, minus
// strand) and a nested source of its domains, which lists a ranges node on the transcript's residues up front and
// serves its rows by region. Two regions of 750: T1 in the first, T2 in the second.

export const Genes = '/data/genes';
export const GeneLength = 1500;

const Transcripts: Record<string, { exons: [number, number][]; domains: [number, number, string][] }> = {
  T1: {
    exons: [
      [100, 400],
      [480, 900],
    ],
    domains: [
      [0, 60, 'A'],
      [90, 200, 'B'],
    ],
  },
  T2: {
    exons: [
      [801, 900],
      [960, 1422],
    ],
    domains: [[10, 150, 'C']],
  },
};

/** The nested source of a transcript's domains: on the transcript's own system, named by its id. */
/** With `byRegion`, the domains node comes with the region's answer, not up front. */
function domainsSource(id: string, byRegion = false): SeqData {
  const t = Transcripts[id];
  const length = Math.floor(t.exons.reduce((n, [a, b]) => n + b - a, 0) / 3);
  const nested = nodeAddress(Genes, `${id}-domains`);
  const domains: Node = { kind: 'ranges', params: { target: itemAddress(nodeAddress(Genes, id), id), source: nested } };
  return {
    coverage: {
      systems: [{ name: id, length }],
      shapes: [{ system: id, shape: 'ranges' }],
      levels: [1],
      ...(byRegion ? {} : { nodes: [{ key: 'domains', node: domains }] }),
    },
    pairs: [
      {
        request: { system: id, start: 0, end: length, level: 1, shapes: ['ranges'] },
        answer: {
          ok: true,
          tables: {
            ranges: {
              form: 'list',
              start: Int32Array.from(t.domains, (d) => d[0]),
              end: Int32Array.from(t.domains, (d) => d[1]),
              row: Int32Array.from(t.domains, (_, i) => i),
            },
          },
          ...(byRegion ? { nodes: [{ key: 'domains', node: domains }] } : {}),
        },
      },
    ],
  };
}

function transcriptNodes(id: string, format: string, byRegion: boolean): DynamicNode[] {
  const t = Transcripts[id];
  return [
    {
      key: id,
      node: {
        kind: 'ranges',
        params: {
          target: '/data/gene',
          map: 'codon',
          rows: {
            start: t.exons.map((e) => e[0]),
            end: t.exons.map((e) => e[1]),
            id: t.exons.map(() => id),
            strand: t.exons.map(() => '-'),
          },
        },
      },
    },
    {
      key: `${id}-domains`,
      node: { kind: 'source', params: { format, value: seqToValue(domainsSource(id, byRegion)) } },
    },
  ];
}

/** The gene source: its letters by region, each region's answer with the transcripts that start in it. */
export function genesData(format = 'seq', nestedByRegion = false): SeqData {
  const letters = 'ACGT'.repeat(GeneLength / 4);
  const pair = (start: number, end: number, ids: string[]) => ({
    request: { system: 'tp53', start, end, level: 1, shapes: ['letters' as const] },
    answer: {
      ok: true as const,
      tables: { letters: { form: 'letters' as const, start, letters: letters.slice(start, end) } },
      nodes: ids.flatMap((id) => transcriptNodes(id, format, nestedByRegion)),
    },
  });
  return {
    coverage: {
      systems: [{ name: 'tp53', length: GeneLength }],
      shapes: [{ system: 'tp53', shape: 'letters' }],
      levels: [1],
    },
    pairs: [pair(0, 750, ['T1']), pair(750, GeneLength, ['T2'])],
  };
}

/** The gene shown with a track a transcript and, with `domains`, a track for each transcript's domains. */
export function transcriptsSpec(
  options: { domains?: boolean; format?: string; nestedByRegion?: boolean; templates?: boolean } = {},
): Node {
  const b = createSeqSpec({ title: 'transcripts' });
  const genes = b.data.source('genes', {
    format: options.format ?? 'seq',
    value: seqToValue(genesData('seq', options.nestedByRegion)),
  });
  const gene = b.data.sequence('gene', { alphabet: 'dna', source: genes, name: 'tp53' });
  b.views.view('main', (v) => {
    const s = v.section('gene', { of: gene });
    s.track('letters').layer({ type: 'letters' });
    if (options.templates === false) return;
    s.track('transcripts', {
      label: 'Transcript',
      for_each: { source: genes, node_kind: 'ranges', filter: { map: 'codon' } },
    }).layer({ type: 'ranges' });
    if (options.domains)
      s.track('domains', { label: 'Domains', for_each: { source: genes, node_kind: 'ranges', depth: 1 } }).layer({
        type: 'ranges',
      });
  });
  return b.build();
}
