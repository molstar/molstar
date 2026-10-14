/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import { type Issue, KeyKind, type KeyTarget, type SeqSet } from '@molstar/seq-core/model';
import { describe, expect, it } from 'vitest';
import {
  ColorLists,
  compile,
  fromJSON,
  makePaletteFunction,
  normalize,
  paletteValue,
  type SeqViewSpecDoc,
  toJSON,
  validate,
} from '../../src/format/index.js';
import { UNNAMED_ITEM } from '../../src/format/palettes/fields.js';
import { CASES, case6, case7, case9, case10, case12 } from './cases.js';
import { fixtureLoaders, IEP, registries } from './load-fixtures.js';

const errors = (issues: Issue[]) => issues.filter((i) => i.severity === 'error');
const codes = (issues: Issue[]) => [...new Set(issues.map((i) => i.code))].sort();
const t: KeyTarget = { kind: -1, id: -1, row: -1 };
const find = (set: SeqSet, key: string) => (set.byKey(set.strings.find(key), t) ? { ...t } : null);
const annotation = (set: SeqSet, key: string) => {
  const k = find(set, key);
  return k?.kind === KeyKind.annotation ? set.annotations[k.id] : null;
};
const nodeRefs = (doc: SeqViewSpecDoc) => {
  const out: string[] = [];
  const walk = (n: SeqViewSpecDoc['root']) => {
    if (n.ref) out.push(n.ref);
    n.children?.forEach(walk);
  };
  walk(doc.root);
  return out;
};

/** Warnings each case expects on fixtures; everything else is an error-free compile. */
const EXPECTED_WARNINGS: Record<string, string[]> = {
  case3: ['svs.alignment.inserts'],
  case4: ['svs.numbering.descending'],
  case6: ['svs.provider.unavailable'],
  case7: ['svs.provider.unavailable'],
  case9: [],
  case10: ['svs.provider.unavailable'],
  case12: [],
  case15: ['svs.alignment.inserts'],
};

describe.each(Object.entries(CASES))('%s', (name, make) => {
  const reg = registries();
  it('builds, validates and round-trips', () => {
    const doc = make();
    expect(validate(doc, reg)).toEqual([]);
    expect(fromJSON(toJSON(doc))).toEqual({ doc, issues: [] });
    const n = normalize(doc, reg);
    expect(normalize(n, reg)).toEqual(n);
    expect(validate(n, reg)).toEqual([]);
  });

  it('compiles on fixtures', async () => {
    const { set, issues } = await compile(make(), reg, fixtureLoaders());
    expect(errors(issues)).toEqual([]);
    for (const code of codes(issues)) expect([...(EXPECTED_WARNINGS[name] ?? [])]).toContain(code);
    expect(set.molecules.length + set.alignments.length).toBeGreaterThan(0);
  });
});

describe('cases on fixtures', () => {
  const reg = registries();
  const load = async (doc: SeqViewSpecDoc, r = reg) => compile(doc, r, fixtureLoaders());

  it('1: a sequence from uniprot with its numbering; labels and indices; a CSV read where it is used', async () => {
    const { set, issues } = await load(CASES.case1());
    const p53 = set.molecules[find(set, 'p53')?.id ?? -1];
    expect(p53?.length).toBe(393);
    expect(p53?.numberings.byName(set.strings.find('uniprot:P04637'))).not.toBeNull();
    const zinc = annotation(set, 'zinc');
    expect([...(zinc?.start ?? [])]).toEqual([175, 178, 237, 241]);
    const iupred = set.annotations.find((a) => a?.geometry === 'values');
    expect(iupred?.count).toBe(393);
    // The tracks view's data: UniProt's motifs and the variants load, the disorder profile is shaped like p53's.
    expect(issues).toEqual([]);
    expect(annotation(set, 'p53-variants')?.count).toBe(10);
    const types = new Set(Array.from(annotation(set, 'p53-features')?.category ?? [], (c) => set.strings.get(c)));
    expect(types.has('Motif')).toBe(true);
    const d = iupred?.sparse;
    const mean = (a: number, b: number) => Array.from(d?.slice(a, b) ?? []).reduce((x, y) => x + y, 0) / (b - a);
    expect([mean(0, 85) > 0.6, mean(110, 280) < 0.25, mean(300, 320) > 0.6, mean(330, 350) < 0.5]).toEqual([
      true,
      true,
      true,
      true,
    ]);
    expect(find(set, 'p53/unit')?.kind).toBe(KeyKind.unit);
  });

  it('2: one FASTA read by two sequences; IMGT ranges and a Kabat end the numbering lacks', async () => {
    const { set } = await load(CASES.case2());
    expect(set.molecules[find(set, 'vh')?.id ?? -1]?.length).toBe(120);
    const cdrs = annotation(set, 'cdrs');
    expect(cdrs?.count).toBe(6);
    const h1 = annotation(set, 'cdr-h1-kabat');
    // Kabat 31..35B on a domain without 35A or 35B ends at 35: indices 30..34.
    expect([h1?.start[0], h1?.stop[0]]).toEqual([30, 35]);
  });

  it("3: Stockholm with a tree; sites in each row's numbering read from a TSV; species per row", async () => {
    const { set } = await load(CASES.case3());
    const a = set.alignments[find(set, 'sh2')?.id ?? -1];
    expect(a?.rowCount).toBe(4);
    expect(a?.tree?.leafOrder.length).toBe(4);
    const sites = annotation(set, 'sites');
    expect(sites?.count).toBe(3);
    expect([...(sites?.target.row ?? [])]).toEqual([0, 2, 1]);
    expect(sites?.links.length).toBe(1);
    expect(annotation(set, 'species')?.count).toBe(4);
    expect(annotation(set, 'cons-a')?.geometry).toBe('values');
    expect(annotation(set, 'profile')?.matrix).not.toBeNull();
  });

  it("4: SS_cons from the alignment's own file and RNAfold's pairs on a row, both base pairs", async () => {
    const { set, issues } = await load(CASES.case4());
    const ss = annotation(set, 'ss-cons');
    expect(ss?.count).toBe(4);
    expect(ss?.pairs?.bStart[0]).toBe(11);
    expect(set.strings.get(ss?.category?.[0] ?? -1)).toBe('base_pair');
    const fold = annotation(set, 'ss-fold');
    expect(fold?.count).toBe(5);
    expect(fold?.target.row[0]).toBe(3);
    expect(issues.filter((i) => i.code === 'svs.numbering.descending')).toHaveLength(1);
  });

  it('5: DNA with features from a TSV (strand a field) and boundaries', async () => {
    const { set } = await load(CASES.case5());
    expect(set.molecules[find(set, 'insert')?.id ?? -1]?.alphabet).toBe('dna');
    const feats = annotation(set, 'features');
    expect(feats?.fields.map((f) => set.strings.get(f.name))).toEqual(['strand']);
    const cuts = annotation(set, 'cuts');
    expect([...(cuts?.start ?? [])]).toEqual([...(cuts?.stop ?? [])]);
    // The default palette value (seqspec:format#themes): the kind of feature, so each kind its own color; no `category`
    // column at all (the cuts): the item, so no track is left uncolored.
    const by = (a: typeof feats) => Array.from({ length: a?.count ?? 0 }, (_, r) => a && paletteValue(a, set, r, null));
    expect(by(feats)).toEqual(['start_codon', 'restriction_site', 'CDS', 'restriction_site']);
    expect(cuts?.category).toBeNull();
    const items = by(cuts);
    expect(new Set(items).size).toBe(items.length);
    expect(feats && paletteValue(feats, set, 0, 'label')).toBe('start codon');
  });

  it("the default palette value of items with no category, id or label: one value, the palette's first color", async () => {
    const b = createSeqViewSpec({ title: 'Unnamed' });
    const a = b.data.sequence({ ref: 'a', value: 'MKTAYIAKQRQISFVKSHFSRQ' });
    const unnamed = b.data.annotation({ ref: 'unnamed', geometry: 'ranges', target: a, value: { start: [2, 5, 9] } });
    const named = b.data.annotation({
      ref: 'named',
      geometry: 'ranges',
      target: a,
      value: { start: [2, 5], label: ['Catalytic H', 'Ligand contact'] },
    });
    const v = b.view({ ref: 'v', label: 'Sites' });
    v.block({ sequence: a });
    v.track({ label: 'Unnamed' }).points({ data: unnamed });
    v.track({ label: 'Named' }).points({ data: named });
    const { set, issues } = await load(b.getDocument());
    expect(errors(issues)).toEqual([]);
    const by = (key: string) => {
      const ann = annotation(set, key);
      return Array.from({ length: ann?.count ?? 0 }, (_, r) => ann && paletteValue(ann, set, r, null));
    };
    const values = by('unnamed');
    expect(values).toEqual([UNNAMED_ITEM, UNNAMED_ITEM, UNNAMED_ITEM]);
    const f = makePaletteFunction({ kind: 'categorical' }, values);
    expect(values.map(f)).toEqual(values.map(() => ColorLists.Category10[0]));
    // Names still key the palette (a legend lists them).
    expect(by('named')).toEqual(['Catalytic H', 'Ligand contact']);
  });

  it('6: chains and SDF poses from a provider; isoforms; the pocket; the extension provider unavailable', async () => {
    const { set, issues } = await load(case6());
    const poses = [0, 1, 2, 3].map((k) => find(set, `poses/${k}/ligand`)?.id ?? -1);
    const groups = new Set(poses.map((u) => set.units.group[u]));
    expect(groups.size).toBe(1);
    expect(find(set, 'abl1a-1iep')?.kind).toBe(KeyKind.mapping);
    const contacts = annotation(set, 'contacts');
    expect(contacts?.count).toBe(IEP_CONTACTS);
    expect(annotation(set, 'pocket')?.count).toBeGreaterThan(0);
    expect(issues.filter((i) => i.code === 'svs.provider.unavailable').map((i) => i.path)).toHaveLength(1);
  });

  it("7: copies, a precursor's segment, alignment, UniProt features, secondary structure", async () => {
    const { set } = await load(case7());
    const A = find(set, '1brs/A')?.id ?? -1;
    const B = find(set, '1brs/B')?.id ?? -1;
    expect(set.units.group[A]).toBe(set.units.group[B]);
    expect(set.mappings.filter(Boolean)).toHaveLength(2);
    const align = set.mappings.find((m) => set.strings.get(m?.key ?? -1) === 'mapping-2');
    expect(set.strings.get(align?.provenance.source ?? -1)).toBe('biomsa');
    expect(annotation(set, 'feat-P00648')?.count).toBe(18);
    expect(annotation(set, 'ss')?.count).toBeGreaterThan(0);
    expect(annotation(set, 'contacts')?.count).toBe(43);
    expect(annotation(set, 'effects')?.matrix?.states.length).toBe(20);
  });

  it('7 without the provider: 1brs and ss unavailable, nothing cascades', async () => {
    const { set, issues } = await load(case7(), registries(false));
    expect(errors(issues)).toEqual([]);
    const unavailable = issues.filter((i) => i.code === 'svs.provider.unavailable').map((i) => i.message.split('"')[1]);
    expect(unavailable.sort()).toEqual(['1brs', 'hyd', 'pdbe', 'ss']);
    expect(find(set, '1brs')).toBeNull();
    expect(find(set, 'barnase/unit')?.kind).toBe(KeyKind.unit);
    expect(find(set, 'barstar/unit')?.kind).toBe(KeyKind.unit);
    expect(annotation(set, 'contacts')?.count).toBe(0);
    expect(annotation(set, 'feat-P11540')?.count).toBeGreaterThan(0);
  });

  it("8: a CIF category from the structure's file, chains as suffixes", async () => {
    const { set, issues } = await load(CASES.case8());
    expect(issues).toEqual([]);
    const plddt = annotation(set, 'plddt');
    expect(plddt?.count).toBe(110 + 47 + 90);
    expect(plddt?.slots.unit.length).toBe(2);
    // The synthetic scores are smooth, high in the folded cores and low at the termini and barnase's signal peptide.
    const v = Array.from(plddt?.sparse ?? []);
    const mean = (a: number, b: number) => v.slice(a, b).reduce((x, y) => x + y, 0) / (b - a);
    expect(Math.max(...v.slice(1).map((x, i) => Math.abs(x - v[i])))).toBeLessThan(12);
    expect([mean(0, 47) < 70, mean(100, 105) > 90, mean(157 + 60, 157 + 64) > 85, v[246] < 75]).toEqual([
      true,
      true,
      true,
      true,
    ]);
    // The model's helices and strands, and the interface the contacts between the chains give.
    const ss = annotation(set, 'ss');
    expect(set.strings.get(ss?.category?.[0] ?? -1)).toBe('helix');
    expect(ss?.count).toBe(10 + 9);
    expect(annotation(set, 'contacts')?.count).toBe(43);
    expect(annotation(set, 'iface')?.count).toBeGreaterThan(0);
  });

  it('9: without a progressive loader the alignment is unavailable, and what reads it with it', async () => {
    const { set, issues } = await compile(case9(), reg, { load: fixtureLoaders().load });
    expect(codes(issues)).toEqual(['svs.provider.unavailable']);
    expect(find(set, 'pf00042')).toBeNull();
    expect(find(set, '1a3n/A')?.kind).toBe(KeyKind.unit);
    expect(annotation(set, 'ss')?.count).toBe(2);
  });

  it('10: the drafts: a unit over part of a sequence, an MSA by value, snapshots', async () => {
    const { set } = await load(case10());
    const kd = find(set, 'abl1-kd')?.id ?? -1;
    expect([set.units.rangeStart[kd], set.units.rangeStop[kd]]).toEqual([241, 493]);
    expect(set.alignments[find(set, 'gatekeeper-msa')?.id ?? -1]?.parent[1]).toBe(0);
    expect(annotation(set, 'kd-sites')?.count).toBe(4);
  });

  it('12: a complex by value with no extension; a pasted Stockholm; one binding reaching both chains', async () => {
    const r = registries(false);
    const { set, issues } = await load(case12(), r);
    expect(issues).toEqual([]);
    const A = find(set, 'dimer.A')?.id ?? -1;
    expect(set.units.group[A]).toBe(set.units.group[find(set, 'dimer.B')?.id ?? -1]);
    expect(set.complexes[find(set, 'dimer')?.id ?? -1]?.units.length).toBe(3);
    const map = set.mappings[find(set, 'homologs.row-0')?.id ?? -1];
    expect(map?.method).toBe('pairwise');
    // Dense values are no rows: one array per target (seqspec:model#annotations), a value per column.
    const cons = annotation(set, 'cons');
    expect([cons?.count, cons?.dense?.[0].filter((x) => !Number.isNaN(x)).length]).toEqual([0, 55]);
  });

  it("13: UniProt's features and the Proteins API's variants on one sequence, by type in its categories", async () => {
    const { set, issues } = await load(CASES.case13());
    expect(issues).toEqual([]);
    const feats = annotation(set, 'p53-features');
    const types = new Set(Array.from(feats?.category ?? [], (c) => set.strings.get(c)));
    expect([...types]).toEqual(expect.arrayContaining(['Chain', 'Domain', 'Binding site', 'Compositional bias']));
    const variants = annotation(set, 'p53-variants');
    expect(variants?.count).toBe(10);
    expect(variants?.start[2]).toBe(174); // R175H
  });
});

const IEP_CONTACTS = (IEP.poses as { contacts: unknown[] }[]).reduce((n, p) => n + p.contacts.length, 0);

describe('generated refs', () => {
  it("match the cases'", () => {
    const six = nodeRefs(case6());
    expect(six).toEqual(
      expect.arrayContaining(['chains.main', 'chains.block', 'chains.block-2', '1iep.src', 'poses.src']),
    );
    const seven = nodeRefs(case7());
    expect(seven).toEqual(
      expect.arrayContaining(['complex-one.residues', 'complex-each.residues', 'mapping', 'mapping-2']),
    );
    expect(seven).toEqual(
      expect.arrayContaining(['complex-one.structure-1brs', 'complex-one.ruler', 'feat-P00648.src']),
    );
    expect(nodeRefs(case9())).toEqual(expect.arrayContaining(['pf00042-hba', 'msa.main', 'msa.chains']));
    expect(nodeRefs(case12())).toEqual(expect.arrayContaining(['homologs.row-0', 'homologs.src']));
    const all = [six, seven].map((r) => new Set(r).size === r.length);
    expect(all).toEqual([true, true]);
  });
});

describe('loading', () => {
  const reg = registries(false);
  const seqDoc = (
    more: (
      b: ReturnType<typeof createSeqViewSpec>,
      p: ReturnType<ReturnType<typeof createSeqViewSpec>['data']['sequence']>,
    ) => void,
  ) => {
    const b = createSeqViewSpec();
    const p = b.data.sequence({ ref: 'p53', source: { service: 'uniprot', resource: 'sequence', id: 'P04637' } });
    more(b, p);
    b.view({ ref: 'v', label: 'V' }).track({ label: 'L' }).letters();
    return b.getDocument();
  };

  it('index and label forms of one range load equal; mixing them is svs.position.form', async () => {
    const doc = seqDoc((b, p) => {
      b.data.annotation({ ref: 'i', geometry: 'ranges', target: p, value: { start: [99, 174], stop: [120, 175] } });
      b.data.annotation({
        ref: 'l',
        geometry: 'ranges',
        target: p,
        numbering: 'uniprot:P04637',
        value: { beg: [100, 175], end: [120, 175] },
      });
    });
    const { set, issues } = await compile(doc, reg, fixtureLoaders());
    expect(issues).toEqual([]);
    const [i, l] = [annotation(set, 'i'), annotation(set, 'l')];
    expect([...(l?.start ?? [])]).toEqual([...(i?.start ?? [])]);
    expect([...(l?.stop ?? [])]).toEqual([...(i?.stop ?? [])]);
    const mixed = seqDoc((b, p) => {
      b.data.annotation({
        ref: 'm',
        geometry: 'ranges',
        target: p,
        numbering: 'ordinal',
        value: { start: [1], beg: [3] },
      });
    });
    const r = await compile(mixed, reg, fixtureLoaders());
    expect(codes(errors(r.issues))).toEqual(['svs.position.form']);
  });

  it('mixing forms in a file is svs.position.form at load', async () => {
    const doc = seqDoc((b, p) => {
      b.data.annotation({
        ref: 'f',
        geometry: 'ranges',
        target: p,
        numbering: 'ordinal',
        source: { value: 'start\tbeg\n1\t\n\t5\n', format: 'tsv' },
      });
    });
    const { issues } = await compile(doc, reg, fixtureLoaders());
    expect(codes(errors(issues))).toEqual(['svs.position.form']);
    expect(issues[0].path).toBe('/root/children/0/children/3');
  });

  it('a source given by value: FASTA, read by record', async () => {
    const b = createSeqViewSpec();
    const src = b.data.source({ ref: 'pasted', value: '>a first\nMKV\n>b\nMKVLA\n', format: 'fasta' });
    b.data.sequence({ ref: 'a', source_ref: src, record: 'b' });
    b.data.sequence({ ref: 'c', source_ref: src, record: 0 });
    b.data.sequence({ ref: 'd', source_ref: src });
    b.view({ ref: 'v', label: 'V' }).track({ label: 'L' }).letters();
    const loaders = fixtureLoaders();
    const { set, issues } = await compile(b.getDocument(), reg, loaders);
    expect(set.molecules[find(set, 'a')?.id ?? -1]?.length).toBe(5);
    expect(set.molecules[find(set, 'c')?.id ?? -1]?.length).toBe(3);
    expect(issues.map((i) => [i.code, i.path])).toEqual([
      ['svs.source.record', '/root/children/0/children/3/params/record'],
    ]);
    expect(loaders.requests).toEqual([]);
  });

  it('a source over its bound, an unregistered format, unreadable bytes, a missing file', async () => {
    const b = createSeqViewSpec();
    b.data.sequence({ ref: 'a', source: { value: '>a\nMKV\n', format: 'fasta', max_bytes: 4 } });
    b.data.sequence({ ref: 'b', source: { value: '>a\nMKV\n', format: 'genbank' } });
    b.data.sequence({ ref: 'c', source: { host_file: 'nothing.fasta', format: 'fasta' } });
    b.data.sequence({ ref: 'd', source: { value: 'garbage', format: 'fasta' } });
    b.view({ ref: 'v', label: 'V' }).track({ label: 'L' }).letters();
    const { set, issues } = await compile(b.getDocument(), reg, fixtureLoaders());
    expect(codes(issues)).toEqual([
      'svs.provider.unavailable',
      'svs.source.format',
      'svs.source.load',
      'svs.source.size_bound',
    ]);
    // Bytes the format cannot read are one issue: the source is unreadable, its reader unavailable quietly.
    expect(find(set, 'd')).toBeNull();
    // A format no package registered: the source is unavailable (a warning), its reader quietly so
    // (seqspec:sources#data-by-reference).
    const unregistered = issues.filter((i) => i.message.includes('genbank'));
    expect(unregistered.map((i) => [i.code, i.severity])).toEqual([['svs.provider.unavailable', 'warning']]);
    expect(find(set, 'b')).toBeNull();
    // A registered format that cannot read the bytes is svs.source.format, an error.
    expect(issues.filter((i) => i.code === 'svs.source.format').map((i) => i.severity)).toContain('error');
  });

  it('an align mapping per aligner; a ligand where a polymer is needed is svs.ligand.unsupported', async () => {
    const b = createSeqViewSpec();
    const a = b.data.sequence({ ref: 'a', value: 'MKVLAAGIVGLLLAQ' });
    const c = b.data.sequence({ ref: 'c', value: 'MKVLGIVGLLLAQ' });
    const lig = b.data.ligand({ ref: 'lig', value: ['C', 'N', 'O'] });
    for (const aligner of reg.list('aligners')) {
      b.data.mapping({ ref: `m-${aligner.name}`, from: a, to: c, method: 'align', aligner: aligner.name });
    }
    b.data.mapping({ ref: 'bad', from: lig, to: c, method: 'align' });
    b.view({ ref: 'v', label: 'V' }).track({ label: 'L' }).letters();
    const { set, issues } = await compile(b.getDocument(), reg, fixtureLoaders());
    expect(reg.list('aligners').map((x) => x.name)).toEqual(['biomsa', 'needleman_wunsch', 'myers']);
    for (const aligner of reg.list('aligners')) {
      const m = set.mappings[find(set, `m-${aligner.name}`)?.id ?? -1];
      expect(m?.forward.filter((x) => x >= 0).length).toBe(13);
      expect(m?.provenance.parameters).toMatchObject({ aligner: aligner.name });
    }
    expect(issues.map((i) => [i.code, i.path])).toEqual([['svs.ligand.unsupported', '/root/children/0/children/6']]);
  });

  it("an aligner's options are checked by its table, and a missing aligner leaves the mapping unavailable", () => {
    const b = createSeqViewSpec();
    const a = b.data.sequence({ ref: 'a', value: 'MKV' });
    b.data.mapping({
      ref: 'm',
      from: a,
      to: a,
      method: 'align',
      aligner: 'myers',
      options: { equalities: 'N', mode: 'x' },
    });
    b.view({ ref: 'v', label: 'V' }).track({ label: 'L' }).letters();
    expect(validate(b.getDocument(), reg).map((i) => [i.code, i.path.split('/').pop()])).toEqual([
      ['svs.params.type', 'mode'],
      ['svs.params.type', 'equalities'],
    ]);
  });

  it('Stockholm tags: GR values at gap columns dropped, GS as row fields', async () => {
    const sto = '# STOCKHOLM 1.0\n#=GS r1 OS Mus\nr0 AC-GU\nr1 ACCGU\n#=GR r0 PP 12345\n#=GR r1 PP 12345\n//\n';
    const b = createSeqViewSpec();
    const aln = b.data.alignment({ ref: 'aln', source: { value: sto, format: 'stockholm' } });
    b.data.annotation({
      ref: 'pp',
      geometry: 'values',
      target: aln,
      from_source: aln,
      tag: 'GR:PP',
      fields: { index: 'start' },
    });
    b.data.annotation({ ref: 'os', geometry: 'rows', target: aln, from_source: aln, tag: 'GS:OS' });
    b.view({ ref: 'v', label: 'V' }).block({ alignment: aln }).track({ label: 'L' }).rows();
    const { set, issues } = await compile(b.getDocument(), reg, fixtureLoaders());
    expect(issues).toEqual([]);
    const pp = annotation(set, 'pp');
    // Dense: per row its values at its residues (gap columns dropped), no rows.
    expect([pp?.count, pp?.dense?.map((d) => [...d])]).toEqual([
      0,
      [
        [1, 2, 4, 5],
        [1, 2, 3, 4, 5],
      ],
    ]);
    expect(annotation(set, 'os')?.fields.map((f) => set.strings.get(f.name))).toEqual(['OS']);
  });

  it('CSV with quotes and JSON tables in both forms', async () => {
    const b = createSeqViewSpec();
    const p = b.data.sequence({ ref: 'p', value: 'MKVLAAG' });
    b.data.annotation({
      ref: 'csv',
      geometry: 'ranges',
      target: p,
      source: { value: 'start,label\n0,"a, ""b"""\n3,c\n', format: 'csv' },
    });
    b.data.annotation({
      ref: 'rows',
      geometry: 'ranges',
      target: p,
      source: { value: '[{"start":1},{"start":2}]', format: 'json' },
    });
    b.data.annotation({
      ref: 'cols',
      geometry: 'ranges',
      target: p,
      source: { value: '{"start":[1,2],"category":{"dict":["x"],"codes":[0,null]}}', format: 'json' },
    });
    b.view({ ref: 'v', label: 'V' }).track({ label: 'L' }).letters();
    const { set, issues } = await compile(b.getDocument(), reg, fixtureLoaders());
    expect(issues).toEqual([]);
    const csv = annotation(set, 'csv');
    expect(set.strings.get(csv?.rowLabel?.[0] ?? -1)).toBe('a, "b"');
    expect(annotation(set, 'rows')?.count).toBe(2);
    expect(set.strings.get(annotation(set, 'cols')?.category?.[0] ?? -1)).toBe('x');
  });

  it("a host's shared parse is used instead of the bytes", async () => {
    const loaders = {
      ...fixtureLoaders(),
      async parsed() {
        return { kind: 'molecules' as const, molecules: [{ key: 'x', residues: 'MKVLA' }] };
      },
    };
    const b = createSeqViewSpec();
    b.data.sequence({ ref: 'a', source: { host_file: 'anything', format: 'fasta' } });
    b.view({ ref: 'v', label: 'V' }).track({ label: 'L' }).letters();
    const { set } = await compile(b.getDocument(), reg, loaders);
    expect(set.molecules[find(set, 'a')?.id ?? -1]?.length).toBe(5);
  });
});
