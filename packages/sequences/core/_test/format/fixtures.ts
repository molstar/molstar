/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Documents for the V1 tests: one valid document with every node kind, and helpers to reach its nodes.

import { createRegistries, KernelRegistry, MvsRegistry, register, type SvsRegistries } from '../../src/format/index.js';

// biome-ignore lint/suspicious/noExplicitAny: test documents are mutated freely as JSON.
export type J = any;

export function kernel(): SvsRegistries {
  const reg = createRegistries();
  register(reg, [KernelRegistry, MvsRegistry]);
  return reg;
}

const layer = (ref: string, params: J, children?: J[]) => ({
  kind: 'layer',
  ref,
  params,
  ...(children ? { children } : {}),
});
const track = (ref: string, ...layers: J[]) => ({ kind: 'track', ref, children: layers });

/** Every node kind, valid. Paths in the tests index into it: keep the order. */
export function fullDoc(): J {
  return {
    format: 'seq_view_spec',
    version: '1.0.0',
    revision: 'r1',
    metadata: { title: 'Every kind' },
    root: {
      kind: 'root',
      params: { default_view: 'v' },
      children: [
        {
          kind: 'data',
          children: [
            /* 0 */ { kind: 'source', ref: 'fasta', params: { value: '>a\nMKV\n', format: 'fasta' } },
            /* 1 */ {
              kind: 'sequence',
              ref: 'p53',
              params: { value: 'MEEPQSDPSVEPPLSQ', alphabet: 'protein' },
              children: [
                {
                  kind: 'numbering',
                  ref: 'p53.kabat',
                  params: { name: 'kabat', method: 'labels', labels: [1, '1A', 2, null] },
                },
              ],
            },
            /* 2 */ { kind: 'sequence', ref: 'vh', params: { source_ref: 'fasta', record: 'a' } },
            /* 3 */ { kind: 'ligand', ref: 'sti', params: { value: ['C', 'N', 'O'], chem_comp_ids: ['STI'] } },
            /* 4 */ {
              kind: 'unit',
              ref: 'u1',
              params: { sequence: 'p53', range: [0, 10], origin: 'designed' },
              custom: { mvs: { label_asym_id: 'A' }, 'org.example': { anything: 1 } },
            },
            /* 5 */ {
              kind: 'complex',
              ref: 'dimer',
              params: { label: 'Dimer' },
              children: [
                { kind: 'unit', ref: 'dimer.A', params: { sequence: 'p53', label: 'A' } },
                { kind: 'unit', ref: 'dimer.L', params: { ligand: 'sti' } },
              ],
            },
            /* 6 */ {
              kind: 'alignment',
              ref: 'aln',
              params: { value: ['MK-V', 'MKLV'], labels: ['a', 'b'] },
              children: [{ kind: 'tree', ref: 'aln.tree', params: { value: '(a,b);' } }],
            },
            /* 7 */ { kind: 'mapping', ref: 'map', params: { from: 'aln/0', to: 'p53/unit', method: 'identity' } },
            /* 8 */ {
              kind: 'annotation',
              ref: 'regions',
              params: { geometry: 'ranges', target: 'p53', value: [{ start: 0, stop: 4, label: 'TA' }] },
            },
            /* 9 */ {
              kind: 'annotation',
              ref: 'sites',
              params: { geometry: 'ranges', target: 'p53', numbering: 'kabat', value: { beg: [1, '1A'] } },
            },
            /* 10 */ {
              kind: 'computed',
              ref: 'cons',
              params: { provider: 'column_stats', alignment: 'aln', measure: 'conservation' },
            },
            /* 11 */ {
              kind: 'annotation',
              ref: 'contacts',
              params: {
                geometry: 'pairs',
                value: { target: ['p53'], start: [0], b_target: ['sti/unit'], b_start: [1], category: ['hbond'] },
              },
            },
          ],
        },
        {
          kind: 'view',
          ref: 'v',
          params: { label: 'Chains' },
          children: [
            {
              kind: 'section',
              ref: 'v.main',
              params: { numbering: 'kabat', positions: 'regions' },
              children: [
                /* 0 */ { kind: 'block', ref: 'v.b1', params: { sequence: 'p53' } },
                /* 1 */ { kind: 'block', ref: 'v.b2', params: { complex: 'dimer' } },
                /* 2 */ track('v.letters', layer('v.letters.l', { type: 'letters' })),
                /* 3 */ {
                  kind: 'group',
                  ref: 'v.g',
                  params: { label: 'Features', collapsed: true },
                  children: [
                    track(
                      'v.ranges',
                      layer('v.ranges.l', {
                        type: 'ranges',
                        data: ['regions', 'sites'],
                        color: { theme: 'field_palette', params: { palette: { kind: 'categorical', colors: 'Set2' } } },
                        glyph: { theme: 'uniform', params: { glyph: { char: '★', label: 'site' } } },
                      }),
                    ),
                  ],
                },
                /* 4 */ track(
                  'v.own',
                  layer('v.own.l', { type: 'values', domain: [0, 1] }, [
                    {
                      kind: 'annotation',
                      ref: 'v.own.data',
                      params: { geometry: 'values', target: 'p53', value: { value: [0.1, 0.2, null] } },
                    },
                  ]),
                ),
                /* 5 */ { kind: 'focus', ref: 'v.focus', params: { target: 'p53', start: 0, stop: 5 } },
              ],
            },
          ],
        },
        {
          kind: 'view',
          ref: 'msa',
          params: { label: 'Alignment', layout: 'tracks' },
          children: [
            {
              kind: 'section',
              ref: 'msa.main',
              children: [
                { kind: 'block', ref: 'msa.b', params: { alignment: 'aln' } },
                track(
                  'msa.rows',
                  layer('msa.rows.l', {
                    type: 'rows',
                    color: { theme: 'conservation_shading', params: { data: 'cons', threshold: 0.4 } },
                  }),
                ),
              ],
            },
          ],
        },
      ],
    },
    snapshots: [{ ref: 's1', title: 'One', ops: [{ op: 'set', ref: 'v.main', params: { numbering: null } }] }],
  };
}

/** The smallest valid document: a sequence and a view of its letters. */
export function minimalDoc(): J {
  return {
    format: 'seq_view_spec',
    version: '1.0.0',
    root: {
      kind: 'root',
      children: [
        { kind: 'data', children: [{ kind: 'sequence', ref: 'p53', params: { value: 'MEEPQSDPSV' } }] },
        {
          kind: 'view',
          ref: 'v',
          params: { label: 'V' },
          children: [{ kind: 'section', ref: 'v.main', children: [track('v.t', layer('v.t.l', { type: 'letters' }))] }],
        },
      ],
    },
  };
}

export const dataOf = (d: J): J[] => d.root.children[0].children;
export const sectionOf = (d: J): J[] => d.root.children[1].children[0].children;
export { layer, track };
