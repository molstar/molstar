/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { SeqSets } from '@molstar/seq-core/model';
import { describe, expect, it } from 'vitest';
import {
  float,
  KernelRegistry,
  type LayerTypeEntry,
  OptionalField,
  type ProviderEntry,
  RequiredField,
  ref,
  register,
  SimpleParamsSchema,
  SvsTreeSchema,
  validate,
} from '../../src/format/index.js';
import { dataOf, fullDoc, type J, kernel, layer, minimalDoc, sectionOf, track } from './fixtures.js';

const errors = (d: unknown, reg = kernel()) => validate(d, reg).filter((i) => i.severity === 'error');

describe('valid documents', () => {
  it('the document with every kind validates with no issues', () => {
    expect(validate(fullDoc(), kernel())).toEqual([]);
  });

  /** A minimal document per kind: the smallest addition that holds that kind. */
  const perKind: Record<string, (d: J) => void> = {
    root: () => {},
    data: () => {},
    view: () => {},
    section: () => {},
    track: () => {},
    layer: () => {},
    sequence: () => {},
    source: (d) =>
      dataOf(d).push({ kind: 'source', ref: 's', params: { service: 'uniprot', id: 'P04637', resource: 'sequence' } }),
    ligand: (d) =>
      dataOf(d).push({ kind: 'ligand', ref: 'hem', params: { value: ['FE', 'C'], chem_comp_ids: ['HEM'] } }),
    numbering: (d) => {
      dataOf(d)[0].children = [
        { kind: 'numbering', ref: 'n', params: { name: 'uniprot:P04637', method: 'offset', first: 1 } },
      ];
    },
    unit: (d) => dataOf(d).push({ kind: 'unit', ref: 'u', params: { sequence: 'p53', observed: [[0, 5]] } }),
    complex: (d) => {
      dataOf(d).push({ kind: 'source', ref: 'pdb', params: { uri: 'https://example.org/1brs.cif', format: 'mmcif' } });
      dataOf(d).push({ kind: 'complex', ref: '1brs', params: { provider: 'molstar:structure', source_ref: 'pdb' } });
      d.root.children[1].children[0].children.unshift({ kind: 'block', ref: 'v.b', params: { unit: '1brs/A' } });
    },
    alignment: (d) =>
      dataOf(d).push({ kind: 'alignment', ref: 'a', params: { value: ['MK', 'M-'], flavor: 'aligned' } }),
    tree: (d) =>
      dataOf(d).push({
        kind: 'alignment',
        ref: 'a',
        params: { value: ['MK', 'M-'] },
        children: [{ kind: 'tree', ref: 'a.t', params: { value: '(0,1);', leaves: 'order' } }],
      }),
    mapping: (d) => {
      dataOf(d).push({ kind: 'sequence', ref: 'q', params: { value: 'MEEPQ' } });
      dataOf(d).push({ kind: 'mapping', ref: 'm', params: { from: 'q', to: 'p53', method: 'pairwise', cigar: '5M' } });
    },
    annotation: (d) =>
      dataOf(d).push({
        kind: 'annotation',
        ref: 'an',
        params: { geometry: 'matrix', target: 'p53', states: ['A', 'B'], value: { index: [0], value: [[0.1, null]] } },
      }),
    computed: (d) => dataOf(d).push({ kind: 'computed', ref: 'cov', params: { provider: 'coverage', of: 'p53/unit' } }),
    block: (d) => sectionOf(d).unshift({ kind: 'block', ref: 'v.b', params: { sequence: 'p53', label: 'p53' } }),
    group: (d) =>
      sectionOf(d).push({
        kind: 'group',
        ref: 'v.g',
        params: { label: 'G' },
        children: [track('v.g.t', layer('v.g.l', { type: 'ruler' }))],
      }),
    focus: (d) =>
      sectionOf(d).push({ kind: 'focus', ref: 'v.f', params: { target: 'p53', beg: 2, end: 4, numbering: 'ordinal' } }),
  };

  it("the table of kinds is the schema's", () => {
    expect(Object.keys(perKind).sort()).toEqual(Object.keys(SvsTreeSchema.nodes).sort());
  });

  for (const [kind, add] of Object.entries(perKind)) {
    it(`a document with a ${kind} validates`, () => {
      const d = minimalDoc();
      add(d);
      expect(validate(d, kernel())).toEqual([]);
      expect(JSON.stringify(d)).toContain(`"kind":"${kind}"`);
    });
  }

  it('a newer minor reads unknowns as warnings; strict reads them as errors', () => {
    const d = fullDoc();
    d.version = '1.3.0';
    d.root.children[1].params.zoom = 2;
    dataOf(d).push({ kind: 'widget', ref: 'w' });
    const issues = validate(d, kernel());
    expect(issues.map((i) => [i.code, i.severity, i.path])).toEqual([
      ['svs.version.newer_minor', 'warning', '/version'],
      ['svs.version.newer_minor', 'warning', '/root/children/0/children/12/kind'],
      ['svs.version.newer_minor', 'warning', '/root/children/1/params/zoom'],
    ]);
    expect(errors(d).length).toBe(0);
    expect(
      validate(d, kernel(), { strict: true })
        .filter((i) => i.severity === 'error')
        .map((i) => i.code),
    ).toEqual(['svs.node.kind', 'svs.params.unknown']);
  });
});

const D = '/root/children/0/children';
const S = '/root/children/1/children/0/children';
const RANGES = `${S}/3/children/0/children/0/params`;

/** [code, path, mutation of the full document]: one per static code. */
const QUERY = { seqql: '1.0', expr: 'p53', text: 'p53' };
const SLICE = { head: 'seq.modifier.slice', args: { '0': 'p53', '1': [9] } };
const SEARCH = { head: 'seq.generator.search', args: { '0': 'GXGXXG' } };
const INPUT = { head: 'seq.generator.input', args: { '0': 'selection' } };
const pushQuery = (d: J, query: unknown) =>
  dataOf(d).push({ kind: 'computed', ref: 'q', params: { provider: 'query', query } });
const mutations: [string, string, (d: J) => void][] = [
  ['svs.format.invalid', '/format', (d) => (d.format = 'mvs')],
  ['svs.version.missing', '/version', (d) => delete d.version],
  ['svs.version.unsupported_major', '/version', (d) => (d.version = '2.0.0')],
  ['svs.node.kind', `${D}/12/kind`, (d) => dataOf(d).push({ kind: 'widget', ref: 'w' })],
  ['svs.node.parent', `${D}/12`, (d) => dataOf(d).push({ kind: 'block', ref: 'b', params: { sequence: 'p53' } })],
  ['svs.node.count', '/root/children/2', (d) => (d.root.children[2].children = [])],
  [
    'svs.node.count',
    '/root/children/0/children/6',
    (d) => dataOf(d)[6].children.push({ kind: 'tree', ref: 't2', params: { value: '(a,b);' } }),
  ],
  ['svs.params.missing', '/root/children/1/params/label', (d) => delete d.root.children[1].params.label],
  ['svs.params.missing', `${D}/8/params/value`, (d) => (dataOf(d)[8].params.value = [{ label: 'x' }])],
  ['svs.params.unknown', `${D}/1/params/colour`, (d) => (dataOf(d)[1].params.colour = 'red')],
  ['svs.params.type', '/root/children/1/params/layout', (d) => (d.root.children[1].params.layout = 'grid')],
  ['svs.params.type', `${D}/7/params/method`, (d) => (dataOf(d)[7].params.method = 'magic')],
  ['svs.params.type', `${D}/1/children/0/params/name`, (d) => (dataOf(d)[1].children[0].params.name = 'ordinal')],
  ['svs.params.exclusive', `${D}/1/params`, (d) => (dataOf(d)[1].params.source_ref = 'fasta')],
  ['svs.params.exclusive', `${D}/5`, (d) => (dataOf(d)[5].params.provider = 'molstar:structure')],
  ['svs.ref.syntax', `${D}/2/ref`, (d) => (dataOf(d)[2].ref = 'v h')],
  ['svs.ref.syntax', `${D}/8/params/target`, (d) => (dataOf(d)[8].params.target = 'p53//A')],
  ['svs.ref.missing', `${D}/2`, (d) => delete dataOf(d)[2].ref],
  ['svs.ref.duplicate', `${D}/2/ref`, (d) => (dataOf(d)[2].ref = 'p53')],
  ['svs.ref.duplicate', '/snapshots/1/ref', (d) => d.snapshots.push({ ref: 's1' })],
  ['svs.ref.unresolved', `${D}/8/params/target`, (d) => (dataOf(d)[8].params.target = 'nope')],
  [
    'svs.ref.unresolved',
    `${D}/11/params/value/b_target/0`,
    (d) => (dataOf(d)[11].params.value.b_target = ['nope/unit']),
  ],
  ['svs.ref.wrong_kind', `${S}/0/params/sequence`, (d) => (sectionOf(d)[0].params.sequence = 'aln/0')],
  ['svs.ref.wrong_kind', `${S}/1/params/complex`, (d) => (sectionOf(d)[1].params.complex = 'p53')],
  ['svs.custom.namespace', `${D}/1/custom/svs`, (d) => (dataOf(d)[1].custom = { svs: {} })],
  ['svs.custom.invalid', `${D}/4/custom/mvs/label_asym_id`, (d) => (dataOf(d)[4].custom.mvs.label_asym_id = 5)],
  ['svs.custom.invalid', `${D}/1/custom/mvs`, (d) => (dataOf(d)[1].custom = { mvs: { label_asym_id: 'A' } })],
  ['svs.requires.unknown', '/requires/0', (d) => (d.requires = ['org.example'])],
  ['svs.position.syntax', `${D}/9/params/value/beg/1`, (d) => (dataOf(d)[9].params.value.beg = [1, '1 A'])],
  [
    'svs.position.syntax',
    `${D}/1/children/0/params/labels/1`,
    (d) => (dataOf(d)[1].children[0].params.labels[1] = 'N13'),
  ],
  ['svs.position.form', `${D}/9/params/value`, (d) => (dataOf(d)[9].params.value.start = [0, 1])],
  ['svs.position.form', `${D}/9/params/value`, (d) => delete dataOf(d)[9].params.numbering],
  ['svs.position.form', `${S}/5/params/numbering`, (d) => (sectionOf(d)[5].params = { target: 'p53', beg: 1, end: 3 })],
  ['svs.columns.length', `${D}/9/params/value/label`, (d) => (dataOf(d)[9].params.value.label = ['one'])],
  ['svs.columns.unknown', `${D}/8/params/value`, (d) => (dataOf(d)[8].params.value[0].row = 1)],
  ['svs.columns.type', `${D}/8/params/value/0/start`, (d) => (dataOf(d)[8].params.value[0].start = -1)],
  [
    'svs.columns.type',
    `${D}/11/params/value/category/codes/0`,
    (d) => (dataOf(d)[11].params.value.category = { dict: ['hb'], codes: [3] }),
  ],
  [
    'svs.color.invalid',
    `${RANGES}/color/params/palette/colors`,
    (d) => {
      sectionOf(d)[3].children[0].children[0].params.color.params.palette.colors = 'Set2x';
    },
  ],
  [
    'svs.color.invalid',
    `${RANGES}/color/params/palette/colors/1/0`,
    (d) => {
      sectionOf(d)[3].children[0].children[0].params.color.params.palette = {
        kind: 'continuous',
        colors: [
          ['red', 0],
          ['#ggg', 1],
        ],
      };
    },
  ],
  [
    'svs.color.invalid',
    `${RANGES}/color/theme`,
    (d) => (sectionOf(d)[3].children[0].children[0].params.color = { theme: 'rainbow' }),
  ],
  [
    'svs.scheme.unknown',
    `${RANGES}/color/params/name`,
    (d) => (sectionOf(d)[3].children[0].children[0].params.color = { theme: 'scheme', params: { name: 'rainbow' } }),
  ],
  [
    'svs.scheme.letter',
    `${RANGES}/color/params/custom/groups/0/letters`,
    (d) => {
      const custom = { kind: 'letters', alphabets: ['dna'], groups: [{ letters: 'AJ', color: 'red' }] };
      sectionOf(d)[3].children[0].children[0].params.color = { theme: 'scheme', params: { custom } };
    },
  ],
  [
    'svs.glyph.invalid',
    `${RANGES}/glyph/params/glyph/char`,
    (d) => {
      sectionOf(d)[3].children[0].children[0].params.glyph.params.glyph = { char: 'ab' };
    },
  ],
  [
    'svs.source.uri',
    `${D}/0/params/uri`,
    (d) => (dataOf(d)[0].params = { uri: 'http://example.org/a.fasta', format: 'fasta' }),
  ],
  ['svs.layer.data', `${S}/2/children/0/params/data`, (d) => (sectionOf(d)[2].children[0].params.data = 'regions')],
  [
    'svs.layer.data',
    `${S}/4/children/0/params/data`,
    (d) => (sectionOf(d)[4].children[0].children[0].params.geometry = 'pairs'),
  ],
  [
    'svs.layer.fallback',
    `${S}/2/children/0/params/fallback`,
    (d) => (sectionOf(d)[2].children[0].params.type = 'org.example:fancy'),
  ],
  [
    'svs.layer.fallback',
    `${RANGES}/color/fallback`,
    (d) => {
      sectionOf(d)[3].children[0].children[0].params.color = { theme: 'org.example:heat' };
    },
  ],
  [
    'svs.layer.annotation',
    `${S}/4/children/0/params/data`,
    (d) => (sectionOf(d)[4].children[0].params.data = 'regions'),
  ],
  // Split rows: by a column an inline table lacks ("sites" has no label), a second split layer, one beside rows.
  [
    'svs.layer.split',
    `${RANGES}/split`,
    (d) => (sectionOf(d)[3].children[0].children[0].params.split = { by: 'label' }),
  ],
  [
    'svs.params.type',
    `${RANGES}/split/empty`,
    (d) => (sectionOf(d)[3].children[0].children[0].params.split = { by: 'label', empty: 'none' }),
  ],
  [
    'svs.layer.split',
    `${S}/3/children/0`,
    (d) => {
      const split = { type: 'points', data: 'regions', split: { by: 'label' } };
      sectionOf(d)[3].children[0].children.push(layer('v.p1', split), layer('v.p2', split));
    },
  ],
  [
    'svs.layer.split',
    '/root/children/2/children/0/children/1',
    (d) =>
      d.root.children[2].children[0].children[1].children.push(layer('v.p', { type: 'points', split: { by: 'x' } })),
  ],
  [
    'svs.section.axis',
    '/root/children/2/children/0',
    (d) => {
      d.root.children[2].children[0].children.unshift({ kind: 'block', ref: 'msa.p53', params: { sequence: 'p53' } });
    },
  ],
  [
    'svs.section.positions',
    '/root/children/1/children/0/params/positions/1',
    (d) => {
      dataOf(d).push({
        kind: 'annotation',
        ref: 'meta',
        params: { geometry: 'rows', target: 'aln', value: { row: [0, 1] } },
      });
      d.root.children[1].children[0].params.positions = ['regions', 'meta'];
    },
  ],
  [
    'svs.tree.streamed',
    `${D}/6/children/0/params/leaves`,
    (d) => {
      dataOf(d).push({
        kind: 'source',
        ref: 'tiles',
        params: { uri: 'https://example.org/pf', protocol: 'seqstar-tiles/1' },
      });
      dataOf(d)[6].params = { source_ref: 'tiles', row_count: 2, column_count: 4 };
    },
  ],
  // QD1: the query provider's (seqspec:seqql#in-documents)
  ['svs.query.invalid', `${D}/12/params/query/expr/args/1/0`, (d) => pushQuery(d, { ...QUERY, expr: SLICE })],
  ['svs.query.version', `${D}/12/params/query/seqql`, (d) => pushQuery(d, { ...QUERY, seqql: '2.0' })],
  ['svs.query.search', `${D}/12/params/query/expr`, (d) => pushQuery(d, { ...QUERY, expr: SEARCH })],
  ['svs.query.input', `${D}/12/params/query/expr`, (d) => pushQuery(d, { ...QUERY, expr: INPUT })],
];

describe('one mutation per static code', () => {
  it.each(mutations)('%s at %s', (code, path, mutate) => {
    const d = fullDoc();
    mutate(d);
    const issues = errors(d);
    expect(issues.map((i) => `${i.code} ${i.path}`)).toContain(`${code} ${path}`);
    for (const i of issues) expect(i.message.length).toBeGreaterThan(10);
  });

  it('covers every static code the validator reports', () => {
    const covered = new Set(mutations.map(([code]) => code));
    const statics = (KernelRegistry.issueCodes ?? [])
      .filter((c) => c.severity === 'error' && !LATER.has(c.code))
      .map((c) => c.code);
    expect(statics.filter((c) => !covered.has(c) && c !== 'svs.ref.cycle')).toEqual([]);
  });

  it('refs that form a cycle are svs.ref.cycle', () => {
    const reg = kernel();
    const echo: ProviderEntry = {
      name: 'org.example:echo',
      options: SimpleParamsSchema({ of: OptionalField(ref('computed'), '', '') }),
      geometry: 'values',
      compute: () => Promise.reject(new Error('not in V1')),
    };
    register(reg, { providers: [echo] });
    const d = fullDoc();
    dataOf(d).push({ kind: 'computed', ref: 'c1', params: { provider: 'org.example:echo', options: { of: 'c2' } } });
    dataOf(d).push({ kind: 'computed', ref: 'c2', params: { provider: 'org.example:echo', options: { of: 'c1' } } });
    expect(errors(d, reg).map((i) => `${i.code} ${i.path}`)).toEqual([`svs.ref.cycle ${D}/12`]);
  });

  it('reports every issue of a document in one pass', () => {
    const d = fullDoc();
    for (const [, , mutate] of mutations.slice(3, 12)) mutate(d);
    expect(new Set(errors(d).map((i) => i.code)).size).toBeGreaterThanOrEqual(5);
  });
});

/** Codes of later phases (load, resolve, render; the theme channel and streamed-tree checks) and of patches and snapshots (V3). */
const LATER = new Set(
  `svs.version.newer_minor svs.ref.derived_unknown svs.position.range svs.section.scroll_group svs.tree.leaves
  svs.patch.base svs.patch.ref svs.patch.anchor svs.patch.cycle svs.snapshot.ref svs.source.load svs.source.format
  svs.source.record svs.source.category svs.source.size_bound svs.sequence.alphabet svs.ligand.unsupported
  svs.numbering.length svs.numbering.unknown svs.numbering.fallback svs.alignment.ragged svs.alignment.shape
  svs.mapping.cigar svs.mapping.range svs.mapping.codon svs.mapping.ambiguous svs.focus.range
  svs.numbering.invalid svs.mapping.align svs.mapping.endpoint svs.section.follow
  svs.theme.channel svs.tree.streamed`.split(/\s+/),
);

describe('edges', () => {
  const msaLayer = (d: J) => d.root.children[2].children[0].children[1].children[0].params;

  it('themes nest at most 4 deep', () => {
    const nest = (depth: number): J =>
      depth === 0 ? { theme: 'letter' } : { theme: 'identity_to_row', params: { base: nest(depth - 1) } };
    const d = fullDoc();
    msaLayer(d).color = nest(3);
    expect(errors(d)).toEqual([]);
    msaLayer(d).color = nest(4);
    expect(errors(d).map((i) => i.code)).toEqual(['svs.color.invalid']);
  });

  it('an extension theme with a kernel fallback validates; its params pass through', () => {
    const d = fullDoc();
    msaLayer(d).color = { theme: 'org.example:heat', params: { k: 1 }, fallback: { theme: 'letter' } };
    msaLayer(d).glyph = { theme: 'field_map', params: { map: { hbond: 'diamond', ionic: { char: '+' } } } };
    expect(errors(d)).toEqual([]);
  });

  it("refs may name a set's keys (data given without a document)", () => {
    const { set } = SeqSets.build({ molecules: [{ key: 'ext', residues: 'MKV' }] });
    const d = minimalDoc();
    d.root.children[1].children[0].children.unshift({ kind: 'block', ref: 'v.b', params: { sequence: 'ext' } });
    expect(errors(d).map((i) => i.code)).toEqual(['svs.ref.unresolved']);
    expect(validate(d, kernel(), { set })).toEqual([]);
  });

  it("an extension layer type's options are checked by its registered table", () => {
    const reg = kernel();
    const type: LayerTypeEntry = {
      ...(KernelRegistry.layerTypes ?? [])[4],
      name: 'org.example:bars',
      options: SimpleParamsSchema({ height: RequiredField(float, '') }),
    };
    register(reg, { layerTypes: [type] });
    const d = fullDoc();
    sectionOf(d)[2].children[0].params = { type: 'org.example:bars', fallback: 'ranges', data: 'regions', options: {} };
    expect(errors(d, reg).map((i) => `${i.code} ${i.path}`)).toEqual([
      `svs.params.missing ${S}/2/children/0/params/options/height`,
    ]);
  });
});
