/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// MessagePack (seqspec:format#values): every example case and the full document round trip JSON → MessagePack → JSON to
// the same text and validate alike, plain and with typed columns; the genome document (genome-fixture.ts, Storybook's
// too) built with typed columns validates in both forms and writes JSON as a replacer of every value would, GC Float32
// too; the builder's float32 twin prints as model's; typed arrays in a document come back as the same kinds and
// numbers; a source reads `format: "msgpack"` data from a host file; the CLIs convert and validate either form; a large
// values table's size.

import { spawnSync } from 'node:child_process';
import { mkdtempSync, readFileSync, writeFileSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { fileURLToPath } from 'node:url';
import { createSeqViewSpec } from '@molstar/seq-core/builder';
import { Annotations, encodeMsgPack, encodeTypedArrays, float32Number, type KeyTarget } from '@molstar/seq-core/model';
import { describe, expect, it } from 'vitest';
import {
  compile,
  fromJSON,
  fromMsgPack,
  isMsgPack,
  readDocument,
  type SeqViewSpecDoc,
  toJSON,
  toMsgPack,
  validate,
} from '../../src/format/index.js';
import { CASES, case17 } from './cases.js';
import { fullDoc } from './fixtures.js';
import { genomeDoc, synthGenome } from './genome-fixture.js';
import { fixtureLoaders, registries } from './load-fixtures.js';

const reg = registries();
type J = Record<string, unknown>;
/** Every case and the full document, case 17 small (8 samples over 20 sites); at its full size (600,000 genotypes) slow. */
const DOCS: [string, () => SeqViewSpecDoc][] = [
  ...Object.entries(CASES).map(([name, make]): [string, () => SeqViewSpecDoc] =>
    name === 'case17' ? [name, () => case17(8, 20)] : [name, make],
  ),
  ['fullDoc', () => fullDoc() as never],
];
const EACH: [string, () => SeqViewSpecDoc, string[]][] = [
  ...DOCS.map(([name, make]): [string, () => SeqViewSpecDoc, string[]] => [name, make, []]),
  ['case17 at full size', () => case17(), ['slow']],
];

describe.each(EACH)('%s', (_, make, tags) => {
  it('JSON → MessagePack → JSON: the same text, deep equal, validated alike', { tags }, () => {
    const text = toJSON(make());
    const json = fromJSON(text).doc as SeqViewSpecDoc;
    for (const typedArrays of [true, false]) {
      const bytes = toMsgPack(json, { typedArrays });
      expect(isMsgPack(bytes)).toBe(true);
      const back = fromMsgPack(bytes);
      expect(back).toEqual({ doc: json, issues: [] });
      expect(toJSON(back.doc as SeqViewSpecDoc)).toBe(text);
      expect(validate(back.doc, reg)).toEqual(validate(json, reg));
      expect(readDocument(bytes).doc).toEqual(json);
    }
    expect(readDocument(new TextEncoder().encode(` ${text}`)).doc).toEqual(json);
  });
});

/** Annotation table columns of numbers as typed arrays: Int32 when whole, else Float64, NaN null. */
function typedColumns(doc: SeqViewSpecDoc): SeqViewSpecDoc {
  const typed = (c: unknown) => {
    if (
      !Array.isArray(c) ||
      !c.some((x) => typeof x === 'number') ||
      !c.every((x) => x === null || typeof x === 'number')
    )
      return c;
    if (c.every((x) => x !== null && (x | 0) === x)) return Int32Array.from(c);
    return Float64Array.from(c, (x) => x ?? Number.NaN);
  };
  const walk = (n: SeqViewSpecDoc['root']) => {
    const v = n.params?.value;
    if (n.kind === 'annotation' && v && typeof v === 'object' && !Array.isArray(v)) {
      for (const [k, c] of Object.entries(v as J)) {
        const col = c as J;
        if (Array.isArray(col?.codes)) col.codes = typed(col.codes);
        else (v as J)[k] = typed(c);
      }
    }
    n.children?.forEach(walk);
  };
  walk(doc.root);
  return doc;
}

describe.each(EACH)('%s, its columns typed', (_, make, tags) => {
  it('JSON → MessagePack → JSON and MessagePack → MessagePack: the same text; decoded, validated as the JSON', {
    tags,
  }, () => {
    const text = toJSON(make());
    const json = fromJSON(text).doc as SeqViewSpecDoc;
    const doc = typedColumns(fromJSON(text).doc as SeqViewSpecDoc);
    expect(toJSON(doc)).toBe(text);
    const bytes = toMsgPack(doc);
    const back = fromMsgPack(bytes);
    expect(back).toEqual({ doc, issues: [] });
    expect(toJSON(back.doc as SeqViewSpecDoc)).toBe(text);
    const again = toMsgPack(back.doc as SeqViewSpecDoc);
    expect(again).toEqual(bytes);
    expect(toJSON(fromMsgPack(again).doc as SeqViewSpecDoc)).toBe(text);
    expect(validate(back.doc, reg)).toEqual(validate(json, reg));
  });
});

it('typed columns: Int32 and Float64 columns and codes, in most cases', () => {
  const kinds = (v: unknown): string[] =>
    ArrayBuffer.isView(v)
      ? [v.constructor.name]
      : v !== null && typeof v === 'object'
        ? Object.values(v).flatMap(kinds)
        : [];
  const per = DOCS.map(([, make]) => kinds(fromMsgPack(toMsgPack(typedColumns(make()))).doc));
  expect(per.filter((k) => k.length > 0).length).toBeGreaterThanOrEqual(10);
  expect(new Set(per.flat())).toEqual(new Set(['Int32Array', 'Float64Array']));
});

describe('the genome document, built with typed columns', () => {
  it('validates with 0 errors as JSON and as decoded MessagePack', () => {
    const doc = genomeDoc(synthGenome({ total: 300_000, genes: 300, proteins: 30, variants: 2000, repeats: 1000 }));
    const gc = doc.root.children?.[0].children?.find((n) => n.ref === 'gc.chr1')?.params?.value as J;
    const genes = doc.root.children?.[0].children?.find((n) => n.ref === 'genes')?.params?.value as J;
    expect([gc.value, genes.start, (genes.target as J).codes].map((c) => (c as object).constructor)).toEqual([
      Float64Array,
      Int32Array,
      Uint8Array,
    ]);
    const json = fromJSON(toJSON(doc)).doc;
    const back = fromMsgPack(toMsgPack(doc)).doc;
    expect(
      (back?.root.children?.[0].children?.find((n) => n.ref === 'genes')?.params?.value as J | undefined)?.start,
    ).toEqual(genes.start);
    const errors = (d: unknown) => validate(d, reg).filter((i) => i.severity === 'error');
    expect([errors(json), errors(back)]).toEqual([[], []]);
  });

  it('JSON: the same text as through a replacer of every value, GC as Float64 and Float32; the document untouched', () => {
    // A fifth of the other tests' genome: the replacer's search per value is the run's cost, and size is not the point.
    const doc = genomeDoc(synthGenome({ total: 60_000, genes: 60, proteins: 6, variants: 400, repeats: 200 }));
    const slow = (x: number) => {
      for (let p = 1; ; p++)
        if (p === 9 || Math.fround(Number(x.toPrecision(p))) === x) return Number(x.toPrecision(p));
    };
    const replacer = (_: string, v: unknown) =>
      ArrayBuffer.isView(v)
        ? Array.from(v as Float64Array, (x) => (Number.isNaN(x) ? null : v instanceof Float32Array ? slow(x) : x))
        : v;
    const gcs = (doc.root.children?.[0].children ?? []).flatMap((n) =>
      n.ref?.startsWith('gc.') ? [n.params?.value as J] : [],
    );
    for (const f32 of [false, true]) {
      if (f32) for (const t of gcs) t.value = Float32Array.from(t.value as []);
      const before = gcs[0].value;
      for (const pretty of [false, true])
        expect(toJSON(doc, { pretty }) === JSON.stringify(doc, replacer, pretty ? 2 : undefined)).toBe(true);
      expect(gcs[0].value).toBe(before);
    }
    expect([gcs.length, gcs[0].value instanceof Float32Array]).toEqual([24, true]);
  });

  it("the builder's float32 twin prints as model's `float32Number`", () => {
    const f = new Float32Array(100_000);
    const bits = new Uint32Array(f.buffer);
    for (let i = 0, s = 88172645; i < bits.length; i++) {
      s ^= s << 13;
      s ^= s >>> 17;
      s ^= s << 5;
      bits[i] = s;
    }
    f.set([0, -0, Number.NaN, Number.POSITIVE_INFINITY, 1e-45, 2 ** -96, 2 ** 87, 2 ** 90, 16777217, 9.0000005e12]);
    const b = createSeqViewSpec({ title: 'Twin' });
    b.data.sequence({ ref: 's', value: 'MKVL' });
    b.view({ ref: 'v', label: 'V' })
      .track({})
      .values({ data: 's', domain: f as never });
    const find = (n: SeqViewSpecDoc['root']): unknown => n.params?.domain ?? n.children?.map(find).find((d) => d);
    const domain = find(b.getDocument().root);
    expect(domain).toEqual(Array.from(f, (x) => (Number.isNaN(x) ? null : float32Number(x))));
  });
});

describe('a Uint8Array through the typed-array extension', () => {
  it('encodes as ext 0x12, not bin', () => {
    const bytes = encodeMsgPack(Uint8Array.of(1, 2, 3, 4), encodeTypedArrays);
    expect(Array.from(bytes)).toEqual([0xd6, 0x12, 1, 2, 3, 4]);
  });
});

describe('typed arrays in a document', () => {
  const t: KeyTarget = { kind: -1, id: -1, row: -1 };
  const value = Float64Array.from({ length: 1000 }, (_, i) => (i % 7 === 0 ? Number.NaN : i / 20));
  const start = Int32Array.from({ length: 10 }, (_, i) => 3 * i);
  function doc() {
    const b = createSeqViewSpec({ title: 'Typed' });
    const s = b.data.sequence({ ref: 's', value: 'ACGT'.repeat(250) });
    // The builder keeps table columns given as lossless typed arrays.
    b.data.annotation({ ref: 'gc', geometry: 'values', target: s, value: { value } });
    b.data.annotation({ ref: 'r', geometry: 'ranges', target: s, value: { start, stop: start.map((x) => x + 2) } });
    b.view({ ref: 'v', label: 'V' }).track({}).values({ data: 'gc' });
    return b.getDocument();
  }

  it('come back as the same kinds and numbers, validate and compile alike, and print as JSON arrays', async () => {
    const d = doc();
    const back = fromMsgPack(toMsgPack(d)).doc as SeqViewSpecDoc;
    expect(back).toEqual(d);
    const cols = (x: SeqViewSpecDoc) => x.root.children?.[0].children?.slice(1).map((n) => n.params?.value);
    expect(cols(d)?.map((c) => Object.values(c as object)[0].constructor)).toEqual([Float64Array, Int32Array]);
    expect(cols(back)?.map((c) => Object.values(c as object)[0].constructor)).toEqual([Float64Array, Int32Array]);
    expect(toJSON(back)).toBe(toJSON(d));
    expect(JSON.parse(toJSON(d)).root.children[0].children[1].params.value.value[1]).toBe(0.05);
    expect([validate(d, reg), validate(back, reg)]).toEqual([[], []]);
    const [a, b] = await Promise.all([d, back].map((x) => compile(x, reg, fixtureLoaders())));
    const gc = (r: typeof a) => (r.set.byName('gc', t) ? r.set.annotations[t.id] : null);
    const exact = (r: typeof a) => Array.from({ length: 1000 }, (_, i) => Annotations.valueAt(gc(r) as never, 0, i));
    expect(exact(b)).toEqual(Array.from(value));
    expect(exact(a)).toEqual(exact(b));
  });

  it("a Float32 column (the producer's choice): its shortest decimals in JSON, Float32 through MessagePack", () => {
    const b = createSeqViewSpec({ title: 'Float32' });
    const s = b.data.sequence({ ref: 's', value: 'MKVL' });
    const score = Float32Array.of(0.15, Number.NaN, 0.3, 1 / 3);
    b.data.annotation({ ref: 'sc', geometry: 'values', target: s, value: { value: score } });
    b.view({ ref: 'v', label: 'V' }).track({}).values({ data: 'sc' });
    const d = b.getDocument();
    const text = toJSON(d);
    expect(text).toContain('"value":{"value":[0.15,null,0.3,0.33333334]}');
    const back = fromMsgPack(toMsgPack(d)).doc as SeqViewSpecDoc;
    const col = (back.root.children?.[0].children?.[1].params?.value as J | undefined)?.value as Float32Array;
    expect([col.constructor, Object.is(col[0], score[0]), toJSON(back)]).toEqual([Float32Array, true, text]);
    expect([validate(d, reg), validate(back, reg)]).toEqual([[], []]);
  });

  it('a large values table: 8 bytes a value with the float64 extension, 9 plain; JSON by the digits', () => {
    const n = 100_000;
    const full = Float64Array.from({ length: n }, (_, i) => Math.sin(i / 1000));
    const short = full.map((x) => Math.round(1000 * x) / 1000);
    const size = (v: Float64Array) => [
      toJSON({ value: v } as never).length,
      toMsgPack({ value: v } as never).length,
      toMsgPack({ value: v } as never, { typedArrays: false }).length,
    ];
    const [json, packed, plain] = size(full);
    // Headers: the map, its key and the ext 32 (13 bytes) or array 32 (12, less 8 for sin(0) = 0, an int).
    expect([packed - 8 * n, plain - 9 * n, Math.round(json / n)]).toEqual([13, 4, 20]);
    // Three decimals print shorter than 8 bytes: there JSON is the smaller file.
    expect(size(short)[0]).toBeLessThan(size(short)[1]);
  });
});

describe('a source in `format: "msgpack"`', () => {
  it('reads a table from a host file, either form, typed arrays or plain', async () => {
    const value = [0.15, 0.3, null, 0.1];
    const files: Record<string, Uint8Array> = {
      'columns.msgpack': encodeMsgPack({ value: Float64Array.of(0.15, 0.3, Number.NaN, 0.1) }, encodeTypedArrays),
      'rows.msgpack': encodeMsgPack(value.map((v, index) => ({ index, value: v }))),
    };
    const base = fixtureLoaders();
    const loaders = {
      ...base,
      load: async (r: Parameters<typeof base.load>[0], signal: AbortSignal) =>
        files[r.hostFile ?? ''] ?? base.load(r, signal),
    };
    const b = createSeqViewSpec({ title: 'MessagePack data' });
    const s = b.data.sequence({ ref: 's', value: 'MKVL' });
    for (const f of Object.keys(files)) {
      const ref = f.split('.')[0];
      b.data.annotation({ ref, geometry: 'values', target: s, source: { host_file: f, format: 'msgpack' } });
    }
    b.view({ ref: 'v', label: 'V' }).track({}).values({ data: 'rows' });
    const d = b.getDocument();
    expect(validate(d, reg)).toEqual([]);
    const { set, issues } = await compile(d, reg, loaders);
    expect(issues).toEqual([]);
    const t: KeyTarget = { kind: -1, id: -1, row: -1 };
    for (const ref of ['columns', 'rows']) {
      const a = set.byName(ref, t) ? set.annotations[t.id] : null;
      const at = (i: number) => Annotations.valueAt(a as never, 0, i);
      expect([0, 1, 2, 3].map(at), ref).toEqual([0.15, 0.3, Number.NaN, 0.1]);
    }
    const bad = { ...loaders, load: async () => Uint8Array.of(0xc1) };
    const failed = await compile(d, reg, bad);
    expect(failed.issues.map((i) => i.message)).toContainEqual(expect.stringMatching(/^not MessagePack/));
  });
});

describe('CLIs', { tags: ['slow'] }, () => {
  const dir = mkdtempSync(join(tmpdir(), 'svs-msgpack-'));
  const cli = (name: string, ...args: string[]) =>
    spawnSync(process.execPath, [fileURLToPath(new URL(`../../lib/format/cli/${name}.js`, import.meta.url)), ...args], {
      encoding: 'utf8',
    });
  const path = (name: string) => join(dir, name);

  it('svs-convert: JSON to MessagePack and back to the same text; svs-validate reads either', () => {
    const text = toJSON(CASES.case1());
    writeFileSync(path('a.svsj'), text);
    expect(cli('svs-convert', path('a.svsj'), path('a.svsm')).status).toBe(0);
    expect(isMsgPack(readFileSync(path('a.svsm')))).toBe(true);
    expect(cli('svs-convert', path('a.svsm'), path('b.svsj')).status).toBe(0);
    expect(readFileSync(path('b.svsj'), 'utf8')).toBe(text);
    const ok = cli('svs-validate', path('a.svsm'));
    expect([ok.status, ok.stderr.trim()]).toEqual([0, `${path('a.svsm')}: 0 error(s), 0 warning(s), 0 note(s)`]);
    writeFileSync(path('bad.svsm'), Uint8Array.of(0x81, 0xc1));
    expect([cli('svs-convert', path('bad.svsm'), path('c.svsj')).status, cli('svs-convert').status]).toEqual([1, 2]);
  });
});
