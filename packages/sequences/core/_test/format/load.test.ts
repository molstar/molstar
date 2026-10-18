/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// VW3's tests (V2, V4, V7): unavailability decided bottom-up and no entry able to reject `compile` (the load half),
// the readers, the transport, load caching, and align mappings with stubbed aligners (the format's half).

import { readFileSync } from 'node:fs';
import { type Aligner, AlignerOptionsError, AlignerSizeError } from '@molstar/seq-core/algorithms';
import { createSeqViewSpec } from '@molstar/seq-core/builder';
import { Fasta, type Issue, type KeyTarget, type RowIndex, type SeqSet } from '@molstar/seq-core/model';
import { tilesStream } from '@molstar/seq-core/stream/tiles';
import { tileServer } from '@molstar/seq-core/stream/tiles/server';
import { describe, expect, it, vi } from 'vitest';
import {
  type AlignerEntry,
  assembleSet,
  compile,
  createRegistries,
  type DocNode,
  indexDoc,
  KernelRegistry,
  type Loaders,
  loadNode,
  MvsRegistry,
  register,
  SimpleParamsSchema,
  type SourceResult,
  type StageContext,
  snapshotDoc,
  sourceParse,
  validate,
} from '../../src/format/index.js';
import { alignerEntry } from '../../src/format/load/entries.js';
import { KernelFormats, UniprotService } from '../../src/format/load/formats.js';
import { CASES } from './cases.js';
import { fixtureLoaders, pf00042Server, registries } from './load-fixtures.js';
import { done, loaders } from './load-support.js';

const errors = (issues: Issue[]) => issues.filter((i) => i.severity === 'error');
const codes = (issues: Issue[]) => issues.map((i) => i.code);
const t: KeyTarget = { kind: -1, id: -1, row: -1 };
const find = (set: SeqSet, key: string) => (set.byKey(set.strings.find(key), t) ? { ...t } : null);
const annotation = (set: SeqSet, key: string) => (find(set, key) ? set.annotations[t.id] : null);
const text = (set: SeqSet, ids: ArrayLike<number> | null | undefined) =>
  Array.from(ids ?? [], (id) => (id < 0 ? null : set.strings.get(id)));
// --- V2: availability, bottom-up ------------------------------------------------------------------------------------

describe('every case and snapshot with no extension registered', () => {
  const reg = registries(false);
  // The data's own warnings (the case convention, a descending numbering) are not the extensions'.
  const allowed = new Set(['svs.provider.unavailable', 'svs.alignment.inserts', 'svs.numbering.descending']);
  const docs = Object.entries(CASES).flatMap(([name, make]) => {
    const doc = make();
    const refs = (doc.snapshots ?? []).map((s) => s.ref);
    return [[name, doc] as const, ...refs.map((ref) => [`${name} ${ref}`, snapshotDoc(doc, ref, reg).doc] as const)];
  });
  it.each(docs)('%s: no errors; only unavailability, once per node', async (_, doc) => {
    const { issues } = await compile(doc, reg, fixtureLoaders());
    expect(errors(issues)).toEqual([]);
    expect(issues.filter((i) => !allowed.has(i.code))).toEqual([]);
    const unavailable = issues.filter((i) => i.code === 'svs.provider.unavailable').map((i) => i.path);
    expect(new Set(unavailable).size).toBe(unavailable.length);
  });
});

describe('a source that does not load or read', () => {
  /** Two sources, each read by a node and by annotations over it: one issue each, whatever reads them. */
  function doc(options: { fasta?: string; sto?: string } = {}) {
    const b = createSeqViewSpec();
    const s1 = b.data.source({
      ref: 's1',
      host_file: 'down.fasta',
      format: 'fasta',
      ...(options.fasta ? { host_file: options.fasta } : {}),
    });
    const p = b.data.sequence({ ref: 'p', source_ref: s1 });
    const s2 = b.data.source({ ref: 's2', host_file: options.sto ?? 'down.sto', format: 'stockholm' });
    const aln = b.data.alignment({ ref: 'aln', source_ref: s2 });
    const start = Array.from({ length: 30 }, (_, i) => i);
    b.data.annotation({ ref: 'rows-on-p', geometry: 'ranges', target: p, value: { start } });
    b.data.annotation({ ref: 'rows-on-aln', geometry: 'ranges', target: aln, value: { start } });
    b.data.mapping({ ref: 'm', from: p, to: p, method: 'identity' });
    b.data.computed({ ref: 'cons', provider: 'column_stats', alignment: aln, measure: 'conservation' });
    b.view({ ref: 'v', label: 'V' }).block({ sequence: p }).track({ label: 'L' }).letters();
    return b.getDocument();
  }

  it('a failed load is one issue per source; the nodes over it are left out quietly', async () => {
    const { set, issues } = await compile(doc(), registries(false), loaders({}));
    expect(issues.map((i) => [i.code, i.path])).toEqual([
      ['svs.source.load', '/root/children/0/children/0'],
      ['svs.source.load', '/root/children/0/children/2'],
    ]);
    for (const key of ['p', 'aln', 'rows-on-p', 'rows-on-aln', 'm', 'cons']) expect(find(set, key)).toBeNull();
  });

  it('a record the source lacks is one issue; its readers are quiet', async () => {
    const b = createSeqViewSpec();
    const src = b.data.source({ ref: 'src', value: '>a\nMKV\n>b\nMKVLA\n', format: 'fasta' });
    const p = b.data.sequence({ ref: 'p', source_ref: src, record: 'nope' });
    b.data.annotation({ ref: 'on-p', geometry: 'ranges', target: p, value: { start: [0, 1, 2] } });
    const { set, issues } = await compile(done(b), registries(false), loaders({}));
    expect(issues.map((i) => i.code)).toEqual(['svs.source.record']);
    expect([find(set, 'p'), find(set, 'on-p')]).toEqual([null, null]);
  });

  it('a source whose bytes the format cannot read is one issue, whoever reads it; its readers are unavailable', async () => {
    const b = createSeqViewSpec();
    const p = b.data.sequence({ ref: 'p', value: 'MKVLAAG' });
    const src = b.data.source({ ref: 'src', host_file: 'bad.json', format: 'json' });
    for (const ref of ['one', 'two']) b.data.annotation({ ref, geometry: 'ranges', target: p, source_ref: src });
    const { set, issues } = await compile(done(b), registries(false), loaders({ 'bad.json': '{ nope' }));
    expect(issues.map((i) => [i.code, i.severity])).toEqual([['svs.source.format', 'error']]);
    expect(issues[0].path).toBe('/root/children/0/children/1');
    // A read that reported an error is unreadable: its readers are unavailable, not empty.
    expect([find(set, 'one'), find(set, 'two')]).toEqual([null, null]);
  });
});

describe("a complex by value leaves out its unavailable units, whatever the document's order", () => {
  /** `a` is unavailable (a format nothing reads); `b` is not. Units are written before or after their sequences. */
  function doc(unitsFirst: boolean, onlyUnavailable = false) {
    const b = createSeqViewSpec();
    const complex = b.data.complex({ ref: 'c' });
    const sequences = () => {
      const ext = b.data.source({ ref: 'ext', value: '>x\nMKV\n', format: 'genbank' });
      b.data.sequence({ ref: 'a', source_ref: ext });
      b.data.sequence({ ref: 'b', value: 'MKVLA' });
    };
    if (!unitsFirst) sequences();
    complex.unit({ ref: 'ua', sequence: 'a' });
    if (!onlyUnavailable) complex.unit({ ref: 'ub', sequence: 'b' });
    if (unitsFirst) sequences();
    b.view({ ref: 'v', label: 'V' }).block({ complex }).track({ label: 'L' }).letters();
    return b.getDocument();
  }

  it.each([false, true])('units first: %s', async (unitsFirst) => {
    const { set, issues } = await compile(doc(unitsFirst), registries(false), fixtureLoaders());
    expect(issues.map((i) => [i.code, i.severity])).toEqual([['svs.provider.unavailable', 'warning']]);
    expect(set.complexes[find(set, 'c')?.id ?? -1]?.units.length).toBe(1);
    expect([find(set, 'ua'), find(set, 'ub')?.kind]).toEqual([null, expect.anything()]);
  });

  it('with no unit left, the complex is unavailable, and nothing else is reported', async () => {
    for (const first of [false, true]) {
      const { set, issues } = await compile(doc(first, true), registries(false), fixtureLoaders());
      expect(codes(issues)).toEqual(['svs.provider.unavailable']);
      expect(find(set, 'c')).toBeNull();
    }
  });

  it('a unit left out for naming a ligand as a sequence is left out of its complex too, quietly', async () => {
    const reg = createRegistries();
    register(reg, [
      KernelRegistry,
      MvsRegistry,
      {
        complexProviders: [
          {
            name: 'org.example:ligands',
            options: SimpleParamsSchema({}),
            formats: ['json'],
            build: async () => ({ molecules: [{ key: 'LIG', kind: 'ligand', elements: ['C'] }], units: [] }),
          },
        ],
      },
    ]);
    const b = createSeqViewSpec();
    b.data.complex({ ref: 'prov', provider: 'org.example:ligands', source: { value: '{}', format: 'json' } } as never);
    const c = b.data.complex({ ref: 'c' });
    c.unit({ ref: 'bad', sequence: 'prov/LIG' } as never);
    c.unit({ ref: 'good', sequence: b.data.sequence({ ref: 's', value: 'MKV' }) });
    const { set, issues } = await compile(done(b), reg, fixtureLoaders());
    expect(issues.map((i) => [i.code, i.path])).toEqual([
      ['svs.ref.wrong_kind', expect.stringMatching(/params\/sequence$/)],
    ]);
    expect(set.complexes[find(set, 'c')?.id ?? -1]?.units.length).toBe(1);
  });
});

describe('entries that cannot run never reject compile', () => {
  const structure = (build: () => never | Promise<never>) => {
    const reg = createRegistries();
    register(reg, [
      KernelRegistry,
      MvsRegistry,
      {
        complexProviders: [
          {
            name: 'org.example:structure',
            options: SimpleParamsSchema({}),
            formats: ['mmcif'],
            build: vi.fn(async () => build()),
          },
        ],
        providers: [
          {
            name: 'org.example:stats',
            options: SimpleParamsSchema({}),
            geometry: 'values',
            compute: vi.fn(async () => build()),
          },
        ],
      },
    ]);
    return reg;
  };
  const complexDoc = (source: Record<string, unknown>) => {
    const b = createSeqViewSpec();
    const c = b.data.complex({ ref: 'c', provider: 'org.example:structure', source: source as never });
    b.data.sequence({ ref: 's', value: 'MKV' });
    b.view({ ref: 'v', label: 'V' }).block({ sequence: 's' }).track({ label: 'L' }).letters();
    return { doc: b.getDocument(), c };
  };

  it("a complex provider never runs on a source its format does not read; a service's format counts", async () => {
    const reg = structure(() => {
      throw new Error('must not run');
    });
    for (const source of [
      { value: '>a\nMKV\n', format: 'fasta' },
      { service: 'uniprot', id: 'P04637', resource: 'sequence' },
    ]) {
      const { doc } = complexDoc(source);
      const { set, issues } = await compile(doc, reg, fixtureLoaders());
      expect(issues.map((i) => [i.code, i.severity])).toEqual([['svs.source.format', 'error']]);
      expect(find(set, 'c')).toBeNull();
    }
    const entry = reg.get('complexProviders', 'org.example:structure');
    expect(entry?.build).not.toHaveBeenCalled();
  });

  it('a provider that throws makes its node unavailable with its issue', async () => {
    const reg = structure(() => {
      throw new Error('boom');
    });
    reg.get('complexProviders', 'org.example:structure');
    register(reg, { formats: [{ ...KernelFormats[0], name: 'mmcif' }] });
    const { doc } = complexDoc({ value: '>a\nMKV\n', format: 'mmcif' });
    const { set, issues } = await compile(doc, reg, fixtureLoaders());
    expect(issues.map((i) => [i.code, i.severity, i.message])).toEqual([
      ['svs.provider.unavailable', 'warning', expect.stringContaining('boom')],
    ]);
    expect(find(set, 'c')).toBeNull();
    expect(find(set, 's')?.kind).toBeDefined();
    // The same for an annotation provider.
    const b = createSeqViewSpec();
    b.data.sequence({ ref: 's', value: 'MKV' });
    b.data.computed({ ref: 'x', provider: 'org.example:stats', options: {} } as never);
    const r = await compile(done(b), reg, fixtureLoaders());
    expect(r.issues.map((i) => [i.code, i.severity])).toEqual([['svs.provider.unavailable', 'warning']]);
  });

  it('an abort goes on: it is not a failed provider', async () => {
    const control = new AbortController();
    const reg = structure(() => {
      control.abort();
      throw new Error('not the abort');
    });
    register(reg, { formats: [{ ...KernelFormats[0], name: 'mmcif' }] });
    const { doc } = complexDoc({ value: '>a\nMKV\n', format: 'mmcif' });
    await expect(compile(doc, reg, fixtureLoaders(), control.signal)).rejects.toThrow();
  });
});

// --- A1's line: align mappings ---------------------------------------------------------------------------------------

describe('align mappings', () => {
  const stub = (name: string, align: Aligner['align']): Aligner => ({ name, description: '', params: {}, align });
  const registered = (...aligners: Aligner[]) => {
    const reg = registries(false);
    register(reg, { aligners: aligners.map(alignerEntry) });
    return reg;
  };
  const mapDoc = (aligner: string, provenance?: Record<string, unknown>) => {
    const b = createSeqViewSpec();
    const a = b.data.sequence({ ref: 'a', value: 'MKVLAAGIVGLLLAQ' });
    const c = b.data.sequence({ ref: 'c', value: 'MKVLGIVGLLLAQ' });
    b.data.mapping({ ref: 'm', from: a, to: c, method: 'align', aligner, provenance } as never);
    b.view({ ref: 'v', label: 'V' }).track({ label: 'L' }).letters();
    return b.getDocument();
  };
  const report = async (aligner: Aligner) => {
    const reg = registered(aligner);
    const doc = mapDoc(aligner.name);
    expect(validate(doc, reg)).toEqual([]);
    const { set, issues } = await compile(doc, reg, fixtureLoaders());
    return { set, issues: issues.map((i) => [i.code, i.severity, i.path]) };
  };
  const at = '/root/children/0/children/2';

  it('a size error, a plain RangeError and any other throw are svs.mapping.align', async () => {
    for (const [name, error] of [
      ['org.example:size', new AlignerSizeError('stub', 50000, 50000)],
      ['org.example:range', new RangeError('invalid array length')],
      ['org.example:plain', new Error('out of memory')],
    ] as const) {
      const r = await report(
        stub(name, async () => {
          throw error;
        }),
      );
      // valid sequences an aligner cannot align: a missing piece, a warning (decision 9)
      expect(r.issues, name).toEqual([['svs.mapping.align', 'warning', at]]);
      expect(find(r.set, 'm')).toBeNull();
    }
  });

  it("only a bad option is svs.params.type, at the mapping's options", async () => {
    const r = await report(
      stub('org.example:options', async () => {
        throw new AlignerOptionsError('option x must be an int');
      }),
    );
    expect(r.issues).toEqual([['svs.params.type', 'error', `${at}/params/options`]]);
  });

  it('an aligner that aligns nothing is svs.mapping.align from the entry; compile adds no svs.mapping.empty', async () => {
    const r = await report(stub('org.example:none', async () => null));
    expect(r.issues).toEqual([['svs.mapping.align', 'warning', at]]);
    // An entry of the host's that returns null has reported its reason itself: compile says nothing more.
    const reg = registries(false);
    const own: AlignerEntry = { name: 'org.example:own', options: SimpleParamsSchema({}), align: async () => null };
    register(reg, { aligners: [own] });
    const quiet = await compile(mapDoc('org.example:own'), reg, fixtureLoaders());
    expect(quiet.issues).toEqual([]);
    expect(find(quiet.set, 'm')).toBeNull();
    // Not aligned is unavailable, so what reads the mapping is left out quietly.
    const index = indexDoc(mapDoc('org.example:own'), reg);
    const molecules = ['a', 'c'].map((k) => ({ key: k, residues: String(index.nodes.get(k)?.params.value) }));
    const ctx: StageContext = {
      signal: new AbortController().signal,
      registries: reg,
      set: assembleSet(null, { molecules }, []).set,
      sources: new Map(),
      unavailable: new Set(),
      source: () => null,
      report: () => {},
    };
    const m = index.nodes.get('m');
    expect(m && (await loadNode(m.node, index, ctx, fixtureLoaders())).unavailable).toBe(true);
  });

  it('an entry that throws is svs.mapping.align and never rejects compile; an abort does', async () => {
    const reg = registries(false);
    const control = new AbortController();
    let aborting = false;
    const own: AlignerEntry = {
      name: 'org.example:throws',
      options: SimpleParamsSchema({}),
      align: async () => {
        if (aborting) control.abort();
        throw new Error('boom');
      },
    };
    register(reg, { aligners: [own] });
    const r = await compile(mapDoc('org.example:throws'), reg, fixtureLoaders());
    expect(r.issues.map((i) => [i.code, i.path])).toEqual([['svs.mapping.align', at]]);
    aborting = true;
    await expect(compile(mapDoc('org.example:throws'), reg, fixtureLoaders(), control.signal)).rejects.toThrow();
  });

  it('a given provenance keeps its fields and gains the aligner and its options under `parameters`', async () => {
    const reg = registries(false);
    const given = { source: 'curator', note: 'checked by hand', parameters: { run: 7 } };
    const { set } = await compile(mapDoc('biomsa', given), reg, fixtureLoaders());
    const m = set.mappings[find(set, 'm')?.id ?? -1];
    expect([set.strings.get(m?.provenance.source ?? -1), set.strings.get(m?.provenance.note ?? -1)]).toEqual([
      'curator',
      'checked by hand',
    ]);
    expect(m?.provenance.parameters).toMatchObject({ run: 7, aligner: 'biomsa', options: expect.anything() });
    const plain = (await compile(mapDoc('biomsa'), reg, fixtureLoaders())).set;
    const p = plain.mappings[find(plain, 'm')?.id ?? -1];
    expect([plain.strings.get(p?.provenance.source ?? -1), p?.provenance.parameters]).toMatchObject([
      'biomsa',
      { aligner: 'biomsa' },
    ]);
  });
});

// --- Routed from MD1 and VW1 -----------------------------------------------------------------------------------------

describe('what the model reports reaches the document', () => {
  it("an A3M's HH-suite records are left out with one svs.source.skipped note (Fasta.alignment's issues)", async () => {
    const a3m = '>query\nMKVLAAG\n>ss_pred\nCCHHHCC\n>Consensus\nMKVLAAG\n>hit\nMKVL-AG\n';
    const b = createSeqViewSpec();
    const aln = b.data.alignment({ ref: 'aln', source: { value: a3m, format: 'a3m' }, flavor: 'a3m' });
    b.view({ ref: 'v', label: 'V' }).block({ alignment: aln }).track({ label: 'L' }).rows();
    const { set, issues } = await compile(b.getDocument(), registries(false), fixtureLoaders());
    expect(issues.map((i) => [i.code, i.severity])).toEqual([['svs.source.skipped', 'info']]);
    expect(set.alignments[find(set, 'aln')?.id ?? -1]?.rowCount).toBe(2);
  });

  it('a provider the model cannot run is a warning and leaves its node unavailable (svs.provider.alphabet)', async () => {
    const b = createSeqViewSpec();
    const aln = b.data.alignment({ ref: 'aln', value: ['XJXJ', 'XJXJ'], alphabet: 'generic' });
    b.data.computed({ ref: 'cons', provider: 'column_stats', alignment: aln, measure: 'conservation' });
    b.view({ ref: 'v', label: 'V' }).block({ alignment: aln }).track({ label: 'L' }).rows();
    const { set, issues } = await compile(b.getDocument(), registries(false), fixtureLoaders());
    expect(issues.map((i) => [i.code, i.severity])).toEqual([['svs.provider.alphabet', 'warning']]);
    expect(find(set, 'cons')).toBeNull();
  });

  it('an empty table needs no columns, in either form', async () => {
    const b = createSeqViewSpec();
    const s = b.data.sequence({ ref: 's', value: 'MKV' });
    for (const geometry of ['ranges', 'values', 'pairs'] as const) {
      b.data.annotation({ ref: `rows-${geometry}`, geometry, target: s, value: [] });
      b.data.annotation({ ref: `cols-${geometry}`, geometry, target: s, value: {} });
    }
    const doc = done(b);
    expect(validate(doc, registries(false))).toEqual([]);
    const { set, issues } = await compile(doc, registries(false), fixtureLoaders());
    expect(issues).toEqual([]);
    expect(annotation(set, 'rows-pairs')?.count).toBe(0);
  });
});

// --- V4: readers -----------------------------------------------------------------------------------------------------

describe('UniProt features', () => {
  // The REST entry's shape: an object with the original and the alternatives; the Proteins API's: letters.
  const rest = {
    primaryAccession: 'P04637',
    features: [
      {
        type: 'Natural variant',
        location: { start: { value: 175, modifier: 'EXACT' }, end: { value: 175, modifier: 'EXACT' } },
        description: 'in LFS',
        featureId: 'VAR_006932',
        alternativeSequence: { originalSequence: 'R', alternativeSequences: ['H', 'C'] },
      },
      {
        type: 'Natural variant',
        location: { start: { value: 248, modifier: 'EXACT' }, end: { value: 249, modifier: 'EXACT' } },
        description: 'a deletion',
        featureId: 'VAR_X',
        alternativeSequence: { originalSequence: 'RR' },
      },
      { type: 'Domain', location: { start: { value: 94 }, end: { value: 292 } }, description: 'DNA-binding' },
    ],
  };
  const proteins = {
    features: [
      {
        type: 'VARIANT',
        begin: '175',
        end: '175',
        wildType: 'R',
        alternativeSequence: 'H',
        consequenceType: 'missense',
      },
    ],
  };
  const read = (value: unknown) => {
    const r = UniprotService.after?.({ kind: 'object', value }, 'P04637', 'features');
    if (r?.kind !== 'table') throw new Error('not a table');
    return Object.fromEntries(Object.entries(r.columns).map(([k, v]) => [k, Array.from(v)]));
  };

  it("one row per alternative, the rows of one feature sharing its id, the original kept, a deletion's alt empty", () => {
    expect(read(rest)).toEqual({
      start: [174, 174, 247, 93],
      stop: [175, 175, 249, 292],
      category: ['Natural variant', 'Natural variant', 'Natural variant', 'Domain'],
      label: ['in LFS', 'in LFS', 'a deletion', 'DNA-binding'],
      id: ['VAR_006932', 'VAR_006932', 'VAR_X', null],
      alt: ['H', 'C', '', null],
      original: ['R', 'R', 'RR', null],
    });
  });

  it("the Proteins API's letters stay one row, the wild type its original; no alt, no column", () => {
    expect(read(proteins)).toMatchObject({ alt: ['H'], original: ['R'], label: ['missense'] });
    const plain = read({ features: [{ type: 'Domain', location: { start: { value: 1 }, end: { value: 2 } } }] });
    expect(Object.keys(plain)).toEqual(['start', 'stop', 'category', 'label']);
  });

  it('loads through the service: alt holds letters, never [object Object]', async () => {
    const b = createSeqViewSpec();
    const p = b.data.sequence({ ref: 'p', value: 'M'.repeat(300) });
    b.data.annotation({
      ref: 'f',
      geometry: 'ranges',
      target: p,
      source: { service: 'uniprot', id: 'P04637', resource: 'features' },
    });
    const files = { 'https://rest.uniprot.org/uniprotkb/P04637.json': JSON.stringify(rest) };
    const { set, issues } = await compile(done(b), registries(false), loaders(files));
    expect(issues).toEqual([]);
    const f = annotation(set, 'f');
    const field = (name: string) => f?.fields.find((x) => set.strings.get(x.name) === name);
    expect(text(set, f?.alt)).toEqual(['H', 'C', '', null]);
    expect(text(set, field('original')?.data as Int32Array)).toEqual(['R', 'R', 'RR', null]);
    // The two alternatives of one feature are one item.
    expect([f?.count, f?.itemCount]).toEqual([4, 3]);
  });

  it("variants past the sequence's end: those starting past it left out, those running past it clipped, one note", async () => {
    // Recorded from the Proteins API (variation-P04637-ends.json): stop-lost extensions and frameshifts from 394, and
    // deletions running from 377 and 385 to 400, on P04637's 393 residues.
    const json = readFileSync(new URL('../fixtures/variation-P04637-ends.json', import.meta.url), 'utf8');
    const b = createSeqViewSpec();
    const p = b.data.sequence({ ref: 'p', value: JSON.parse(json).sequence });
    b.data.annotation({
      ref: 'v',
      geometry: 'ranges',
      target: p,
      source: { service: 'uniprot', id: 'P04637', resource: 'variants' },
    });
    const files = { 'https://www.ebi.ac.uk/proteins/api/variation/P04637': json };
    const { set, issues } = await compile(done(b), registries(false), loaders(files));
    expect(issues.map((i) => [i.code, i.severity, i.message])).toEqual([
      [
        'svs.source.skipped',
        'info',
        "uniprot P04637 variants: 3 start past the sequence's end (393) and are left out; 2 run past it and are clipped",
      ],
    ]);
    const v = annotation(set, 'v');
    expect([Array.from(v?.start ?? []), Array.from(v?.stop ?? [])]).toEqual([
      [174, 174, 392, 376, 384],
      [175, 175, 393, 393, 393],
    ]);
  });

  it('a response without the sequence keeps every variant as given', () => {
    const r = UniprotService.after?.(
      { kind: 'object', value: { features: [{ type: 'VARIANT', begin: '394', end: '394', consequenceType: 'x' }] } },
      'P04637',
      'variants',
    );
    expect(r?.kind === 'table' && Array.from(r.columns.start as number[])).toEqual([393]);
  });
});

describe('delimited files: text cells, trimmed, typed after `fields`', () => {
  const csvDoc = (value: string, extra: Record<string, unknown> = {}, format = 'csv') => {
    const b = createSeqViewSpec();
    const p = b.data.sequence({ ref: 'p', value: 'MKVLAAGIVGLLLAQ' });
    b.data.annotation({
      ref: 'f',
      geometry: 'ranges',
      target: p,
      source: { host_file: 'table', format },
      ...extra,
    } as never);
    b.view({ ref: 'v', label: 'V' }).track({ label: 'L' }).letters();
    return { doc: b.getDocument(), files: { table: value } };
  };
  const load = async (value: string, extra: Record<string, unknown> = {}, format = 'csv') => {
    const { doc, files } = csvDoc(value, extra, format);
    return compile(doc, registries(false), loaders(files));
  };

  it('cells are trimmed, so a space after each comma leaves the columns numeric; quoted cells keep theirs', async () => {
    const { set, issues } = await load('start, stop, category, label\n0, 5, 01, " a b "\n7, 9, 1 ,c\n');
    expect(issues).toEqual([]);
    const f = annotation(set, 'f');
    expect([Array.from(f?.start ?? []), Array.from(f?.stop ?? [])]).toEqual([
      [0, 7],
      [5, 9],
    ]);
    expect(text(set, f?.category)).toEqual(['01', '1']);
    expect(text(set, f?.rowLabel)).toEqual([' a b ', 'c']);
  });

  it('ids, categories and labels stay text however they look; positions and values are numbers', async () => {
    const tsv = 'start\tstop\tid\tcategory\tscore\n0\t2\t007\t1E5\t0.50\n3\t4\t008\t1E5\t\n';
    const { set, issues } = await load(tsv, {}, 'tsv');
    expect(issues).toEqual([]);
    const f = annotation(set, 'f');
    expect(text(set, f?.category)).toEqual(['1E5', '1E5']);
    const score = f?.fields.find((x) => set.strings.get(x.name) === 'score');
    expect(Array.from(score?.data as Float64Array)).toEqual([0.5, Number.NaN]);
    expect(text(set, f?.itemKey)).toEqual(['007', '008']);
  });

  it('a repeated header name is svs.source.format, and the file is not read', async () => {
    const { set, issues } = await load('start,stop,start\n0,1,2\n');
    expect(issues.map((i) => [i.code, i.severity])).toEqual([['svs.source.format', 'error']]);
    expect(issues[0].message).toContain('"start"');
    expect(find(set, 'f')).toBeNull();
  });

  it("`fields` read the file's columns as read: a swap, a fan-out, an ignored column", async () => {
    const swapped = await load('start,stop,category,label\n0,2,cat,lab\n', {
      fields: { category: 'label', label: 'category' },
    });
    const f = annotation(swapped.set, 'f');
    expect([text(swapped.set, f?.category), text(swapped.set, f?.rowLabel)]).toEqual([['lab'], ['cat']]);
    const fan = await load('start,stop,name\n0,2,x\n', { fields: { category: 'name', label: 'name', name: null } });
    const g = annotation(fan.set, 'f');
    expect([text(fan.set, g?.category), text(fan.set, g?.rowLabel)]).toEqual([['x'], ['x']]);
  });

  /** The reader's rules a character at a time (the reader slices): the oracle for the random files below. */
  function oracle(s: string, sep: string): string[][] {
    const rows: string[][] = [];
    let [row, cell, quoted, wasQuoted] = [[] as string[], '', false, false];
    const cellDone = () => {
      row.push(wasQuoted ? cell : cell.trim());
      [cell, wasQuoted] = ['', false];
    };
    for (let i = 0; i <= s.length; i++) {
      const c = s[i];
      if (quoted && c !== undefined) {
        if (c !== '"') cell += c;
        else if (s[i + 1] === '"') cell += s[i++];
        else quoted = false;
      } else if (c === '"' && !wasQuoted && /^[ \t]*$/.test(cell)) [cell, quoted, wasQuoted] = ['', true, true];
      else if (c === sep) cellDone();
      else if (c === '\n' || c === '\r' || c === undefined) {
        if (c === '\r' && s[i + 1] === '\n') i++;
        cellDone();
        if (row.length > 1 || row[0] !== '') rows.push(row);
        row = [];
      } else cell += c;
    }
    return rows;
  }

  it('reads 3,000 random files as the character-at-a-time rules do: quotes, line ends, blanks, ragged rows', () => {
    let seed = 7;
    const random = () => {
      seed = (seed * 1103515245 + 12345) & 0x7fffffff;
      return seed / 0x7fffffff;
    };
    const pieces = ['a', 'bc', '1', ' ', '\t', ',', '"', '""', '\n', '\r\n', '\r', 'x y'];
    for (let k = 0; k < 3000; k++) {
      const sep = k % 2 ? ',' : '\t';
      const s = Array.from({ length: Math.floor(random() * 30) }, () => pieces[Math.floor(random() * 12)]).join('');
      const rows = oracle(s, sep);
      const reported: string[] = [];
      const ctx = { report: (i: Issue) => reported.push(i.message) } as never;
      const r = KernelFormats.find((f) => f.name === (sep === ',' ? 'csv' : 'tsv'))?.read(
        new TextEncoder().encode(s),
        {} as never,
        ctx,
      );
      const names = rows[0]?.map((x) => x.trim());
      const ok = names && new Set(names).size === names.length && rows.every((x) => x.length === names.length);
      if (!ok) {
        expect([JSON.stringify(s), reported.length]).toEqual([JSON.stringify(s), 1]);
        continue;
      }
      const columns = Object.fromEntries(names.map((name, c) => [name, rows.slice(1).map((x) => x[c] || null)]));
      expect([JSON.stringify(s), r]).toEqual([
        JSON.stringify(s),
        { kind: 'table', columns, length: rows.length - 1, text: true },
      ]);
    }
  });
});

describe('FASTA and Newick that cannot be read are reported, not guessed', () => {
  const one = async (build: (b: ReturnType<typeof createSeqViewSpec>) => void) => {
    const b = createSeqViewSpec();
    build(b);
    b.view({ ref: 'v', label: 'V' }).track({ label: 'L' }).letters();
    return compile(b.getDocument(), registries(false), fixtureLoaders());
  };

  it('a header with no sequence is a note and a polymer with no residues; no records is an error', async () => {
    const r = await one((b) => b.data.sequence({ ref: 'a', source: { value: '>a\n', format: 'fasta' } }));
    expect(r.issues.map((i) => [i.code, i.severity])).toEqual([
      ['svs.source.empty_record', 'info'],
      ['svs.params.missing', 'error'],
    ]);
    const none = await one((b) => b.data.sequence({ ref: 'a', source: { value: '   ', format: 'fasta' } }));
    expect(none.issues.map((i) => i.code)).toContain('svs.source.format');
    expect(find(none.set, 'a')).toBeNull();
  });

  it('non-ASCII residues are svs.sequence.alphabet', async () => {
    const r = await one((b) => b.data.sequence({ ref: 'a', source: { value: '>a\nMKŁV\n', format: 'fasta' } }));
    expect(codes(errors(r.issues))).toEqual(['svs.sequence.alphabet']);
  });

  it('a Newick tree that does not parse is svs.source.format; the alignment still loads without its tree', async () => {
    const r = await one((b) => {
      const aln = b.data.alignment({ ref: 'aln', value: ['MKV', 'MKL'], labels: ['A', 'B'] });
      aln.tree({ value: '((A,B),(C,D' });
    });
    // One issue, the parse's: no second "no tree", and none for the source giving no tree.
    expect(r.issues.map((i) => [i.code, i.severity])).toEqual([['svs.source.format', 'error']]);
    const fromFile = await one((b) => {
      const aln = b.data.alignment({ ref: 'aln', value: ['MKV', 'MKL'], labels: ['A', 'B'] });
      aln.tree({ source: { value: '((A,B),(C,D', format: 'newick' } } as never);
    });
    expect(fromFile.issues.map((i) => [i.code, i.severity])).toEqual([['svs.source.format', 'error']]);
    expect(r.set.alignments[find(r.set, 'aln')?.id ?? -1]?.tree).toBeNull();
  });
});

describe('Clustal', () => {
  const aln = readFileSync(new URL('../fixtures/pf00042-seed-8.aln', import.meta.url), 'utf8');
  const load = (text: string) => {
    const b = createSeqViewSpec();
    const a = b.data.alignment({ ref: 'aln', source: { host_file: 'seed.aln', format: 'clustal' } });
    b.view({ ref: 'v', label: 'V' }).block({ alignment: a }).track({ label: 'L' }).rows();
    return compile(done(b), registries(false), loaders({ 'seed.aln': text }));
  };

  it('an .aln is an aligned alignment, its rows by name, its blocks joined', async () => {
    const { set, issues } = await load(aln);
    expect(issues).toEqual([]);
    const a = set.alignments[find(set, 'aln')?.id ?? -1];
    expect([a?.rowCount, a?.columnCount]).toEqual([8, 118]);
    expect([a?.labelText(0 as RowIndex), a?.labelText(7 as RowIndex)]).toEqual([
      'HBA_HUMAN/27-137',
      'A0A2I0WAL8_9ASPA/106-197',
    ]);
  });

  it('a file that is not Clustal is svs.source.format, saying what to write, and the alignment is unavailable', async () => {
    const { set, issues } = await load('>a\nMK\n');
    expect(issues.map((i) => [i.code, i.severity])).toEqual([['svs.source.format', 'error']]);
    expect(issues[0].message).toContain('starts with "CLUSTAL"');
    expect(find(set, 'aln')).toBeNull();
  });
});

describe('Stockholm tags', () => {
  const sto = (rows: string[], tags: string[]) => `# STOCKHOLM 1.0\n${rows.join('\n')}\n${tags.join('\n')}\n//\n`;
  const format = KernelFormats.find((f) => f.name === 'stockholm');
  const read = (data: string, tag: string, geometry: 'ranges' | 'values' | 'pairs' | 'rows' | 'matrix') => {
    const issues: Issue[] = [];
    const ctx = { report: (i: Issue) => issues.push(i) } as never;
    const request = { record: null, tag, categoryName: null, blockHeader: null, blockIndex: 0 } as never;
    const r = format?.read(new TextEncoder().encode(data), { ...(request as object), geometry } as never, ctx);
    return { r: r as Awaited<NonNullable<typeof r>>, issues };
  };

  it('GR values at gap columns are dropped, one string per row, typed columns', () => {
    const data = sto(['r0 AC-GU', 'r1 ACCGU'], ['#=GR r0 PP 12345', '#=GR r1 PP 12345']);
    const { r, issues } = read(data, 'GR:PP', 'values');
    expect(issues).toEqual([]);
    if (r.kind !== 'table') throw new Error('table');
    expect(r.columns.value).toBeInstanceOf(Float64Array);
    expect(Array.from(r.columns.value)).toEqual([1, 2, 4, 5, 1, 2, 3, 4, 5]);
    expect(Array.from(r.columns.target)).toEqual(['0', '0', '0', '0', '1', '1', '1', '1', '1']);
  });

  it('a 2k × 300 GR table builds typed columns with no object per position', () => {
    const [rows, width] = [2000, 300];
    const letters = (k: number) => Array.from({ length: width }, (_, c) => 'ACDE'[(k + c) % 4] as string).join('');
    const ss = (k: number) => Array.from({ length: width }, (_, c) => ((c + k) % 20 < 5 ? 'H' : '-')).join('');
    const names = Array.from({ length: rows }, (_, k) => `r${k}`);
    const data = sto(
      names.map((n, k) => `${n} ${letters(k)}`),
      names.map((n, k) => `#=GR ${n} SS ${ss(k)}`),
    );
    const keys = vi.spyOn(Object, 'keys');
    const { r, issues } = read(data, 'GR:SS', 'ranges');
    const calls = keys.mock.calls.length;
    keys.mockRestore();
    expect(issues).toEqual([]);
    if (r.kind !== 'table') throw new Error('table');
    // The old reader asked `Object.keys` for every item (about 30,000 here); the new one for none of them.
    expect(calls).toBeLessThan(rows);
    expect(r.length).toBe(names.reduce((n, _, k) => n + (ss(k).match(/H+/g) ?? []).length, 0));
    expect([r.columns.start, r.columns.stop].map((c) => c.constructor)).toEqual([Int32Array, Int32Array]);
    // One string per row, not per item.
    const targets = Array.from(r.columns.target) as string[];
    expect(new Set(targets).size).toBe(rows);
    expect(targets[0]).toBe(targets[14]);
    expect(Object.is(targets[0], targets[1])).toBe(true);
  });

  it('an unbalanced WUSS column keeps its issue; pairs of a balanced one are base pairs', () => {
    const rows = ['r0 GCGGAUUUAGCU'];
    const bad = read(sto(rows, ['#=GC SS_cons <<<<....>>>.']), 'GC:SS_cons', 'pairs');
    expect(bad.issues.map((i) => [i.code, i.severity])).toEqual([['model.source.format', 'error']]);
    expect(bad.issues[0].message).toContain('SS_cons');
    expect(bad.issues[0].message).toMatch(/WUSS/);
    const good = read(sto(rows, ['#=GC SS_cons <<<<....>>>>']), 'GC:SS_cons', 'pairs');
    if (good.r.kind !== 'table') throw new Error('table');
    expect([Array.from(good.r.columns.start), Array.from(good.r.columns.b_start)]).toEqual([
      [0, 1, 2, 3],
      [11, 10, 9, 8],
    ]);
    expect(Array.from(good.r.columns.category)).toEqual(Array(4).fill('base_pair'));
  });

  it('a values read with no digit, a geometry a tag cannot be read as, and GS as ranges are svs.source.format', () => {
    const data = sto(['r0 ACGU', 'r1 ACGU'], ['#=GC RF xxxx', '#=GS r0 OS mouse']);
    for (const [tag, geometry] of [
      ['GC:RF', 'values'],
      ['GC:RF', 'matrix'],
      ['GC:RF', 'rows'],
      ['GS:OS', 'ranges'],
      ['GC:MISSING', 'ranges'],
    ] as const) {
      const { r, issues } = read(data, tag, geometry);
      expect(
        issues.map((i) => [i.code, i.severity]),
        `${tag} as ${geometry}`,
      ).toEqual([['svs.source.format', 'error']]);
      expect(r.kind === 'table' && r.length).toBe(0);
    }
    const ok = read(data, 'GS:OS', 'rows');
    expect(ok.r.kind === 'table' && Array.from(ok.r.columns.OS)).toEqual(['mouse']);
  });
});

// --- The tile transport (the client's own tests are `@molstar/seq-core/stream/tiles`'s) ---------------------------------------

describe("a statistic's failure", () => {
  const outage = new Error('503 Service Unavailable');
  const serving = (fail: (control: AbortController) => void = () => {}) => {
    const server = pf00042Server();
    const control = new AbortController();
    const loaders: Loaders = {
      load: fixtureLoaders().load,
      stream: tilesStream(() => (path, signal) => {
        if (!path.startsWith('stats/')) return server.transport(path, signal);
        fail(control);
        return Promise.reject(outage);
      }),
    };
    return { loaders, control };
  };

  it("a failed stats request is svs.source.load, not 'serves no stats'; an abort goes on", async () => {
    const { case9 } = await import('./cases.js');
    const { loaders } = serving();
    const { issues } = await compile(case9(), registries(false), loaders);
    const loads = issues.filter((i) => i.code === 'svs.source.load');
    expect(loads.length).toBe(1);
    expect(loads[0].message).toContain('503');
    expect(issues.some((i) => /serves no stats/.test(i.message))).toBe(false);
    const aborting = serving((c) => c.abort());
    await expect(compile(case9(), registries(false), aborting.loaders, aborting.control.signal)).rejects.toThrow();
  });
});

// --- V7: load caching ------------------------------------------------------------------------------------------------

describe('a source is fetched once and parsed once however many nodes read it', () => {
  it('50 records of one FASTA: 1 fetch, 1 parse', async () => {
    const fasta = Array.from({ length: 50 }, (_, k) => `>seq${k}\nMKVLDEFGHIKLMNPQRSTVWY${'K'.repeat(k)}\n`).join('');
    const reads = vi.fn();
    const counting = {
      ...KernelFormats[0],
      read: (...a: Parameters<(typeof KernelFormats)[0]['read']>) => {
        reads();
        return KernelFormats[0].read(...a);
      },
    };
    const reg = createRegistries();
    register(reg, [{ ...KernelRegistry, formats: [counting, ...KernelFormats.slice(1)] }, MvsRegistry]);
    const b = createSeqViewSpec();
    const src = b.data.source({ ref: 'src', host_file: 'proteome.fasta', format: 'fasta' });
    for (let k = 0; k < 50; k++) b.data.sequence({ ref: `s${k}`, source_ref: src, record: k });
    b.view({ ref: 'v', label: 'V' }).track({ label: 'L' }).letters();
    const l = loaders({ 'proteome.fasta': fasta });
    const { set, issues } = await compile(b.getDocument(), reg, l);
    expect(issues).toEqual([]);
    expect([l.requests.length, reads.mock.calls.length]).toEqual([1, 1]);
    expect(set.molecules.filter(Boolean).length).toBe(50);
    expect(set.molecules[find(set, 's49')?.id ?? -1]?.length).toBe(22 + 49);
  });

  it("an annotation's uri is fetched once, whoever reads it (a lazy alignment's named rows, then the annotation)", async () => {
    const server = tileServer({ rowCount: 100, columnCount: 40 });
    const base = loaders({ 'https://example.org/sites.tsv': 'target\tstart\n3\t5\n9\t7\n' });
    const l: Loaders & { requests: string[] } = {
      ...base,
      requests: base.requests,
      stream: tilesStream((uri) => (uri === 'https://example.org/aln' ? server.transport : null)),
    };
    const b = createSeqViewSpec();
    const aln = b.data.alignment({
      ref: 'aln',
      row_count: 100,
      column_count: 40,
      source: { uri: 'https://example.org/aln', protocol: 'seqstar-tiles/1' },
    });
    b.data.annotation({
      ref: 'sites',
      geometry: 'ranges',
      target: aln,
      uri: 'https://example.org/sites.tsv',
      format: 'tsv',
    } as never);
    b.view({ ref: 'v', label: 'V' }).block({ alignment: aln }).track({ label: 'L' }).rows();
    const { set, issues } = await compile(b.getDocument(), registries(false), l);
    expect(issues).toEqual([]);
    expect(l.requests).toEqual(['https://example.org/sites.tsv']);
    expect(annotation(set, 'sites')?.count).toBe(2);
    // The rows its targets name were fetched whole, and nothing else of the table was.
    expect(server.calls.filter((c) => c.startsWith('rows/')).sort()).toEqual(['rows/0', 'rows/3', 'rows/9']);
  });

  it("an annotation's uri is its source for ctx.source, as a source node's is", async () => {
    const b = createSeqViewSpec();
    const p = b.data.sequence({ ref: 'p', value: 'MKVLAAG' });
    b.data.annotation({
      ref: 'sites',
      geometry: 'ranges',
      target: p,
      uri: 'https://example.org/s.tsv',
      format: 'tsv',
    } as never);
    const index = indexDoc(done(b), registries(false));
    const parsed = { kind: 'table', columns: {}, length: 0 } as const;
    const sources = new Map([
      ['sites#uri', { kind: 'loaded', format: 'tsv', parsed, read: async () => parsed } as const],
    ]);
    expect(sourceParse(index, sources, 'sites')).toBe(parsed);
  });

  it('ST2: the rows mappings and tables name come from one pass over the document, however many alignments are lazy', async () => {
    const b = createSeqViewSpec();
    const alns = [0, 1, 2].map((i) =>
      b.data.alignment({
        ref: `a${i}`,
        row_count: 50,
        column_count: 12,
        source: { uri: `https://example.org/a${i}`, protocol: 'seqstar-tiles/1' },
      }),
    );
    for (const [i, a] of alns.entries())
      b.data.annotation({
        ref: `sites${i}`,
        geometry: 'ranges',
        target: a,
        uri: `https://example.org/s${i}.tsv`,
        format: 'tsv',
      } as never);
    const index = indexDoc(done(b), registries(false));
    const server = tileServer({ rowCount: 50, columnCount: 12 });
    const l: Loaders = { load: async () => new Uint8Array(), stream: tilesStream(() => server.transport) };
    let reads = 0;
    const table = { kind: 'table', columns: { target: [String(3)], start: [1] }, length: 1 } as const;
    const sources = new Map<string, SourceResult>(
      [0, 1, 2].map((i) => [
        `sites${i}#uri`,
        {
          kind: 'loaded',
          format: 'tsv',
          parsed: table,
          read: async () => {
            reads++;
            return table;
          },
        },
      ]),
    );
    const ctx: StageContext = {
      signal: new AbortController().signal,
      registries: registries(false),
      set: null,
      sources,
      unavailable: new Set(),
      source: () => null,
      report: () => {},
    };
    for (const i of [0, 1, 2]) {
      const src = index.nodes.get(String(index.nodes.get(`a${i}`)?.params.source_ref));
      const r = await loadNode(src?.node as DocNode, index, ctx, l);
      sources.set(src?.ref as string, r.source as SourceResult);
    }
    const loaded = await Promise.all(
      [0, 1, 2].map((i) => loadNode(index.nodes.get(`a${i}`)?.node as DocNode, index, ctx, l)),
    );
    // Each annotation's table once, not once per alignment; the rows its targets name, fetched whole per alignment.
    expect(reads).toBe(3);
    expect(server.calls.filter((c) => c.startsWith('rows/')).sort()).toEqual([
      'rows/0',
      'rows/0',
      'rows/0',
      'rows/3',
      'rows/3',
      'rows/3',
    ]);
    expect(loaded.map((r) => r.alignments?.[0].lazy?.rows.size)).toEqual([2, 2, 2]);
  });

  it('ST2 review: a shared scan whose caller aborted is redone for another caller, not its failure', async () => {
    const b = createSeqViewSpec();
    for (const i of [0, 1])
      b.data.alignment({
        ref: `a${i}`,
        row_count: 50,
        column_count: 12,
        source: { uri: `https://example.org/a${i}`, protocol: 'seqstar-tiles/1' },
      });
    b.data.annotation({ ref: 'sites', geometry: 'ranges', uri: 'https://example.org/s.tsv', format: 'tsv' } as never);
    const index = indexDoc(done(b), registries(false));
    const server = tileServer({ rowCount: 50, columnCount: 12 });
    const l: Loaders = { load: async () => new Uint8Array(), stream: tilesStream(() => server.transport) };
    const table = { kind: 'table', columns: { target: ['a1/3'], start: [1] }, length: 1 } as const;
    let open = () => {};
    const gate = new Promise<void>((ok) => {
      open = ok;
    });
    // The table's read honors its caller's signal, as a load does.
    const read = async (_: unknown, c: StageContext) => {
      await gate;
      c.signal.throwIfAborted();
      return table;
    };
    const sources = new Map<string, SourceResult>([
      ['sites#uri', { kind: 'loaded', format: 'tsv', parsed: table, read }],
    ]);
    const context = (signal: AbortSignal): StageContext => ({
      signal,
      registries: registries(false),
      set: null,
      sources,
      unavailable: new Set(),
      source: () => null,
      report: () => {},
    });
    const [first, second] = [new AbortController(), new AbortController()];
    for (const i of [0, 1]) {
      const src = index.nodes.get(String(index.nodes.get(`a${i}`)?.params.source_ref));
      sources.set(
        src?.ref as string,
        (await loadNode(src?.node as DocNode, index, context(second.signal), l)).source as SourceResult,
      );
    }
    const node = (i: number) => index.nodes.get(`a${i}`)?.node as DocNode;
    const a0 = loadNode(node(0), index, context(first.signal), l);
    const a1 = loadNode(node(1), index, context(second.signal), l);
    await new Promise((ok) => setTimeout(ok, 0));
    first.abort();
    open();
    await expect(a0).rejects.toThrow();
    expect([...((await a1).alignments?.[0].lazy?.rows.keys() ?? [])]).toEqual([0, 3]);
  });

  it('named rows are fetched eight at a time, in row order', async () => {
    const [n, rowsNamed] = [100, 20];
    const server = tileServer({ rowCount: n, columnCount: 12 });
    let flying = 0;
    let most = 0;
    const done: number[] = [];
    const transport = async (path: string, signal?: AbortSignal) => {
      const bytes = server.transport(path, signal);
      if (!path.startsWith('rows/')) return bytes;
      flying++;
      most = Math.max(most, flying);
      await new Promise((r) => setTimeout(r, 5 + (Number(path.slice(5)) % 3) * 4));
      flying--;
      done.push(Number(path.slice(5)));
      return bytes;
    };
    const l: Loaders = {
      load: async () => new Uint8Array(),
      stream: tilesStream(() => transport),
    };
    const b = createSeqViewSpec();
    const aln = b.data.alignment({
      ref: 'aln',
      row_count: n,
      column_count: 12,
      source: { uri: 'https://example.org/aln', protocol: 'seqstar-tiles/1' },
    });
    const s = b.data.sequence({ ref: 's', value: 'MKVLAAGIVGLL' });
    for (let k = 1; k <= rowsNamed; k++) {
      b.data.mapping({ ref: `m${k}`, from: `${aln.ref}/${k * 3}`, to: s, method: 'identity' } as never);
    }
    b.view({ ref: 'v', label: 'V' }).block({ alignment: aln }).track({ label: 'L' }).rows();
    const { set } = await compile(b.getDocument(), registries(false), l);
    // The reference row, then the twenty named ones: 21 in all, in batches of eight.
    expect(done.length).toBe(rowsNamed + 1);
    expect(most).toBe(8);
    const lazy = set.alignments[find(set, 'aln')?.id ?? -1];
    expect(lazy?.rows.count).toBe(rowsNamed + 1);
    expect(find(set, 'm20')?.kind).toBeDefined();
  });
});

describe("a FASTA alignment's source (GV-02)", () => {
  it('each of its issues is said once, on the source', async () => {
    const fasta = Array.from({ length: 30 }, (_, k) => `>id${k % 3}\nACGT`).join('\n');
    const b = createSeqViewSpec();
    const s = b.data.source({ ref: 'src', value: fasta, format: 'fasta' });
    const aln = b.data.alignment({ ref: 'aln', source_ref: s });
    b.view({ ref: 'v', label: 'V' }).block({ alignment: aln });
    const { issues } = await compile(b.getDocument(), registries(), fixtureLoaders());
    // Was 54: each repeat, per parse, twice (the default parse and the alignment's); one for the cause (GV-07).
    expect(issues.map((i) => `${i.code} ${i.path} ${i.message}`)).toEqual([
      'svs.source.id /root/children/0/children/0 FASTA: 3 id(s) name several records: "id0", "id1", "id2".',
    ]);
  });
});

describe("a file table's bad cells (CF-03)", () => {
  it('say what the same table inline says: one svs.columns.type for the column, and its rows left out', async () => {
    const value = { start: [...Array.from({ length: 200 }, () => 'x'), 1], stop: Array.from({ length: 201 }, () => 3) };
    const tsv = `start\tstop\n${value.start.map((s, i) => `${s}\t${value.stop[i]}`).join('\n')}\n`;
    const issuesOf = async (inline: boolean) => {
      const b = createSeqViewSpec();
      const s = b.data.sequence({ ref: 's', value: 'ACGTACGTAC' });
      b.data.annotation({
        ref: 'a',
        geometry: 'ranges',
        target: s,
        ...(inline ? { value } : { uri: 'https://example.org/a.tsv', format: 'tsv' }),
      } as never);
      b.view({ ref: 'v', label: 'V' }).block({ sequence: s });
      const doc = b.getDocument();
      const reg = registries();
      const all = inline
        ? validate(doc, reg)
        : (await compile(doc, reg, fixtureLoaders({ 'https://example.org/a.tsv': tsv }))).issues;
      const { set } = await compile(doc, reg, fixtureLoaders({ 'https://example.org/a.tsv': tsv }));
      return { said: all.map((i) => `${i.code} ${i.message}`), rows: annotation(set, 'a')?.count ?? null };
    };
    const inline = await issuesOf(true);
    expect(inline.said).toEqual(['svs.columns.type expected an index (int >= 0), got "x" (200 rows)']);
    const file = await issuesOf(false);
    expect(file.said).toEqual(inline.said);
    expect(file.rows).toBe(1);
  });

  it('number cells that are no index (-1, 2.5) say the same as inline too (LA-10)', async () => {
    for (const bad of [-1, 2.5]) {
      const start = [bad, bad, bad, 1];
      const tsv = `start\tstop\n${start.map((s) => `${s}\t3`).join('\n')}\n`;
      const doc = (inline: boolean) => {
        const b = createSeqViewSpec();
        const s = b.data.sequence({ ref: 's', value: 'ACGTACGTAC' });
        const value = { start, stop: start.map(() => 3) };
        b.data.annotation({
          ref: 'a',
          geometry: 'ranges',
          target: s,
          ...(inline ? { value } : { uri: 'https://example.org/a.tsv', format: 'tsv' }),
        } as never);
        b.view({ ref: 'v', label: 'V' }).block({ sequence: s });
        return b.getDocument();
      };
      const reg = registries();
      const said = (all: { code: string; message: string }[]) => all.map((i) => `${i.code} ${i.message}`);
      const inline = said(validate(doc(true), reg));
      expect(inline).toEqual([`svs.columns.type expected an index (int >= 0), got ${bad} (3 rows)`]);
      const file = await compile(doc(false), reg, fixtureLoaders({ 'https://example.org/a.tsv': tsv }));
      expect(said(file.issues), String(bad)).toEqual(inline);
      expect(annotation(file.set, 'a')?.count).toBe(1);
    }
  });
});

describe("the FASTA reader's one parse slot (GV-02, LA-38)", () => {
  const read = KernelFormats.find((f) => f.name === 'fasta')?.read as NonNullable<
    (typeof KernelFormats)[number]['read']
  >;
  const bytesOf = (text: string) => new TextEncoder().encode(text);
  const ctx = { report: () => {} } as never;
  const molecules = { want: 'molecules' } as never;
  const alignments = { want: 'alignments' } as never;

  it('two reads of the same bytes parse once', () => {
    const parse = vi.spyOn(Fasta, 'parse');
    try {
      const a = bytesOf('>x\nACGT\n>y\nACGA');
      read(a, molecules, ctx);
      read(a, alignments, ctx);
      expect(parse).toHaveBeenCalledTimes(1);
    } finally {
      parse.mockRestore();
    }
  });

  it("a new file's parse evicts the last one: only one parse is kept", () => {
    const parse = vi.spyOn(Fasta, 'parse');
    try {
      const a = bytesOf('>x\nACGT');
      const b = bytesOf('>x\nACGT'); // the same text, another identity
      read(a, molecules, ctx);
      read(b, molecules, ctx);
      expect(parse).toHaveBeenCalledTimes(2);
      read(b, alignments, ctx);
      expect(parse).toHaveBeenCalledTimes(2); // b is kept
      read(a, alignments, ctx);
      expect(parse).toHaveBeenCalledTimes(3); // a was evicted by b: parsed again
    } finally {
      parse.mockRestore();
    }
  });
});
