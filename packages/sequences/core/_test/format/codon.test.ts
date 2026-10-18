/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// C1: a protein over its coding DNA (seqspec:format#mappings, seqspec:format#placement, seqspec:format#layer-types):
// the `codon` mapping loaded from the document, protein data placed on the gene's axis as whole codons, `letters` with
// `of`, and the issues.

import { createSeqViewSpec, GENETIC_CODES, translate } from '@molstar/seq-core/builder';
import { Coords, GC_PRT, type Issue, MappingStatus, Paths, type SeqSet } from '@molstar/seq-core/model';
import { describe, expect, it } from 'vitest';
import {
  compile,
  type ResolvedView,
  resolveView,
  type SeqViewSpecDoc,
  supportedAll,
  validate,
} from '../../src/format/index.js';
import { case14 } from './cases.js';
import { fixtureLoaders, registries } from './load-fixtures.js';

const reg = registries();
const codes = (issues: Issue[]) => issues.map((i) => `${i.code} ${i.path}`);
const errors = (issues: Issue[]) => issues.filter((i) => i.severity === 'error');

/** The placements of a layer: [block, annotation, path's span, its runs over rows' ranges]. */
function placed(v: ResolvedView, set: SeqSet, layer: string) {
  const l = Array.from(v.layers.ref).indexOf(set.strings.find(layer));
  const out: unknown[] = [];
  for (let k = v.layers.placementStart[l]; k < v.layers.placementStart[l] + v.layers.placementCount[l]; k++) {
    const path = v.placements.path[k] >= 0 ? v.paths[v.placements.path[k]] : null;
    const ann = set.annotations[v.placements.annotation[k]];
    const runs: number[][] = [];
    for (let r = 0; ann && path && r < ann.count; r++)
      Coords.forEachTo(path.map, ann.start[r], ann.stop[r], (s, t) => runs.push([s, t]), null);
    out.push([v.placements.block[k], v.placements.annotation[k], path?.span ?? 0, runs]);
  }
  return out;
}

/** Case 14's data with a mapping or a layer changed. */
function variant(codon: Record<string, unknown> = {}, letters: Record<string, unknown> = {}): SeqViewSpecDoc {
  const b = createSeqViewSpec();
  const gene = b.data.sequence({ ref: 'gene', alphabet: 'dna', value: 'ATGGCTTCCTAA' });
  b.data.sequence({ ref: 'protein', alphabet: 'protein', value: 'MAS' }); // M and S are also DNA ambiguity codes
  b.data.sequence({ ref: 'other', alphabet: 'dna', value: 'ACGT' });
  b.data.mapping({ ref: 'm', from: 'protein', to: 'gene', method: 'codon', ...codon } as never);
  b.view({ ref: 'v', label: 'V' })
    .block({ sequence: gene })
    .track({ label: 'P' })
    .letters({ of: 'protein', ...letters });
  return b.getDocument();
}

describe('case 14: a protein over its coding DNA', () => {
  it('compiles the translated protein and its codon mapping, every status exact', async () => {
    const { set, issues } = await compile(case14(), reg, fixtureLoaders());
    expect(issues).toEqual([]);
    const protein = set.molecules.find((m) => m && set.strings.get(m.key) === 'protein');
    expect([protein?.alphabet, protein && String.fromCharCode(...protein.code)]).toEqual([
      'protein',
      'MASKELTPRDGHWYFNSQVE',
    ]);
    const m = set.mappings.find((x) => x && set.strings.get(x.key) === 'protein.codon');
    expect([m?.method, m?.span, set.strings.get(m?.provenance.method ?? -1)]).toEqual(['codon', 3, 'translation']);
    expect(Array.from(m?.status ?? []).every((s) => s === MappingStatus.exact)).toBe(true);
  });

  it("places the protein's letters and data on the gene's axis as whole codons; the ruler in nucleotides", async () => {
    const doc = case14();
    const { set } = await compile(doc, reg, fixtureLoaders());
    const { view, issues } = resolveView(doc, set, 'gene-view', reg, supportedAll(reg));
    expect(issues).toEqual([]);
    const ann = (key: string) => set.annotations.findIndex((a) => a && set.strings.get(a.key) === key);
    expect(placed(view, set, 'gene-view.letters-2')).toEqual([[0, -1, 3, []]]);
    expect(placed(view, set, 'gene-view.ranges')).toEqual([[0, ann('exons'), 0, []]]);
    expect(placed(view, set, 'gene-view.ranges-2')).toEqual([
      [
        0,
        ann('domain'),
        3,
        [
          [24, 46],
          [76, 84],
        ],
      ],
    ]);
    expect(placed(view, set, 'gene-view.points')).toEqual([[0, ann('site'), 3, [[90, 93]]]]);
    // His12 split by the intron: 45, then [76, 78).
    const path = view.paths[view.placements.path[view.layers.placementStart[2]]];
    expect([0, 1, 2].map((k) => Paths.spanAt(path, 11, k))).toEqual([45, 76, 77]);
    expect([set.strings.get(view.blocks.numbering[0]?.name ?? -1), view.blocks.numberingPath[0]]).toEqual([
      'ordinal',
      -1,
    ]);
  });

  it("does not carry a protein's numbering across the codon step", async () => {
    const doc = case14();
    const [data, view] = doc.root.children ?? [];
    const section = view.children?.[0];
    const protein = data.children?.find((n) => n.ref === 'protein');
    if (!section || !protein) throw new Error('case 14 changed');
    section.params = { ...section.params, numbering: 'uniprot:X' };
    protein.children = [{ kind: 'numbering', ref: 'protein.x', params: { name: 'uniprot:X', method: 'offset' } }];
    const { set } = await compile(doc, reg, fixtureLoaders());
    const { issues } = resolveView(doc, set, 'gene-view', reg, supportedAll(reg));
    expect(codes(issues)).toEqual(['svs.numbering.fallback /root/children/1/children/0/children/0']);
  });
});

describe('codon mappings and letters with of: issues', () => {
  const load = async (doc: SeqViewSpecDoc) => {
    expect(validate(doc, reg)).toEqual([]);
    return compile(doc, reg, fixtureLoaders());
  };

  it('loads strand, genetic code and a partial construct', async () => {
    const { set, issues } = await load(variant({ cds: [[0, 9]], genetic_code: 11, strand: '+' }));
    expect(issues).toEqual([]);
    expect(set.mappings[0]?.codon?.geneticCode).toBe(11);
  });

  it('reports svs.mapping.codon and svs.mapping.range at the mapping', async () => {
    const at = '/root/children/0/children/3';
    for (const codon of [
      { from: 'other', to: 'gene' },
      {
        cds: [
          [0, 6],
          [3, 9],
        ],
      },
      {
        cds: [
          [3, 9],
          [0, 3],
        ],
      },
      { cds: [[0, 8]] },
      { genetic_code: 99 },
    ]) {
      expect(codes(errors((await load(variant(codon))).issues))).toEqual([`svs.mapping.codon ${at}`]);
    }
    expect(codes(errors((await load(variant({ cds: [[0, 15]] }))).issues))).toEqual([`svs.mapping.range ${at}`]);
  });

  it('keeps of and row apart, and reports letters of a polymer that reaches no block', async () => {
    const both = variant({}, { row: 0 });
    expect(codes(validate(both, reg))).toEqual([
      'svs.params.exclusive /root/children/1/children/0/children/1/children/0/params/of',
    ]);
    const doc = variant({}, { of: 'other' });
    const { set } = await load(doc);
    expect(codes(resolveView(doc, set, 'v', reg, supportedAll(reg)).issues)).toEqual([
      'svs.layer.unplaced /root/children/1/children/0/children/1/children/0/params/of',
    ]);
  });

  it('places `of` as a list or `mapped`: a placement per protein a codon mapping brings to the block', async () => {
    const b = createSeqViewSpec();
    const gene = b.data.sequence({ ref: 'gene', alphabet: 'dna', value: 'ATGGCTTCCTAA' });
    b.data.sequence({ ref: 'p1', alphabet: 'protein', value: 'MAS' });
    b.data.sequence({ ref: 'p2', alphabet: 'protein', value: 'AS' });
    b.data.sequence({ ref: 'other', alphabet: 'dna', value: 'ACGTAC' });
    b.data.mapping({ ref: 'm1', from: 'p1', to: 'gene', method: 'codon', cds: [[0, 9]] } as never);
    b.data.mapping({ ref: 'm2', from: 'p2', to: 'gene', method: 'codon', cds: [[3, 9]] } as never);
    const v = b.view({ ref: 'v', label: 'V' }).block({ sequence: gene }).block({ sequence: 'other' });
    v.track({ label: 'List' }).letters({ ref: 'lst', of: ['p1', 'p2'] });
    v.track({ label: 'Mapped' }).letters({ ref: 'mapped', of: 'mapped' });
    const doc = b.getDocument();
    expect(validate(doc, reg)).toEqual([]);
    const { set } = await load(doc);
    const r = resolveView(doc, set, 'v', reg, supportedAll(reg));
    expect(r.issues).toEqual([]);
    const name = (k: number) => set.strings.get(set.molecules[r.view.paths[k].from.id]?.key ?? -1);
    for (const layer of ['lst', 'mapped']) {
      const l = Array.from(r.view.layers.ref).indexOf(set.strings.find(layer));
      const at = Array.from({ length: r.view.layers.placementCount[l] }, (_, k) => r.view.layers.placementStart[l] + k);
      // Both proteins on the gene's block, each through a codon path from its own space; none on the other block.
      expect(
        at.map((k) => [
          r.view.placements.block[k],
          name(r.view.placements.path[k]),
          r.view.paths[r.view.placements.path[k]].span,
        ]),
      ).toEqual([
        [0, 'p1', 3],
        [0, 'p2', 3],
      ]);
    }
    const bad = variant({}, { of: 42 });
    expect(codes(validate(bad, reg))).toEqual([
      'svs.params.type /root/children/1/children/0/children/1/children/0/params/of',
    ]);
  });

  it("keeps the builder's genetic codes equal to the model's", () => {
    expect(GENETIC_CODES).toEqual(GC_PRT);
  });

  it('translates by the builder what the model reads back as exact', async () => {
    const b = createSeqViewSpec();
    const g = b.data.sequence({ ref: 'g', alphabet: 'rna', value: 'GUGAAAUGAUAG' });
    translate(b.data, { ref: 'p', from: g, cds: [[0, 12]], genetic_code: 11 });
    b.view({ label: 'V' }).block({ sequence: g });
    const { set, issues } = await load(b.getDocument());
    expect(issues).toEqual([]);
    expect(Array.from(set.mappings[0]?.status ?? [])).toEqual([0, 0, 0]);
  });
});
