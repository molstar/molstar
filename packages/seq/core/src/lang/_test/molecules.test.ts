/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { parseJson, printJson } from '@molstar/core/repr/codecs/json';
import { parseSexpr, printSexpr } from '@molstar/core/repr/codecs/sexpr';
import { applyPatch, diff } from '@molstar/core/repr/patch';
import { normalize, validateAndNormalize } from '@molstar/core/repr/schema';
import { hashNode, type Node } from '@molstar/core/repr/tree';
import { readFileSync } from 'node:fs';
import { describe, expect, it } from 'vitest';
import { createSeqSpec, modified } from '../builder.js';
import { checkSeqSpec } from '../check.js';
import { SeqLanguage, SeqLanguageVersion, SeqValidateOptions } from '../language.js';
import { moleculesSpec } from './molecules.js';

const codes = (spec: Node) => checkSeqSpec(spec).map((i) => i.code);

/** A small deterministic generator (mulberry32): property tests replay. */
function rng(seed: number) {
  let a = seed >>> 0;
  return () => {
    a = (a + 0x6d2b79f5) >>> 0;
    let t = a;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

describe('molecules, complexes, units and links', () => {
  const spec = moleculesSpec();

  it('are valid with every new param, the host-only note aside (0.5.0)', () => {
    expect(SeqLanguageVersion).toBe('0.5.0');
    expect(checkSeqSpec(spec)).toEqual([]);
  });

  it('print as the description shows them', () => {
    const text = printSexpr(spec, { head: 64, language: SeqLanguage });
    expect(text).toContain(':modified (:position [0 3] :component [MSE CRO] :length [1 3])');
    expect(text).toContain(':residues [NAG NAG BMA]');
    expect(text).toContain(':atoms (:residue [0 0 0] :name [PG O1G N1] :element [P O N])');
    expect(text).toContain('(compound 2MG :name');
    expect(text).toContain('(unit L :of /data/atp :chain D :group PRD_000001)');
    expect(text).toContain(':occupancy [0.6000000238418579 1]');
    expect(text).toContain(':ref model:7');
    expect(text).toContain(':field plddt');
    expect(text).toContain(':form atoms');
    expect(text).toContain(':modified code');
  });

  it('round-trip through JSON (plain and tagged) and S-expressions, every table typed', () => {
    const plain = printJson(spec, { language: SeqLanguage });
    const a = validateAndNormalize(SeqLanguage, parseJson(plain).node as Node, SeqValidateOptions);
    expect(a.issues).toEqual([]);
    expect(hashNode(a.node as Node)).toBe(hashNode(spec));
    const tagged = parseJson(printJson(spec));
    expect(hashNode(normalize(SeqLanguage, tagged.node as Node, SeqValidateOptions))).toBe(hashNode(spec));
    const back = parseSexpr(printSexpr(spec, { head: 64 }));
    expect(back.issues).toEqual([]);
    expect(hashNode(normalize(SeqLanguage, back.node as Node, SeqValidateOptions))).toBe(hashNode(spec));
    // the declared columns read typed: residues a category, positions int32, occupancy float32
    const data = (spec.children ?? []).find((c) => c.kind === 'data');
    const fc = data?.children?.find((c) => c.id === 'fc');
    expect((fc?.params?.modified as Record<string, unknown>).position).toBeInstanceOf(Int32Array);
    const links = data?.children?.find((c) => c.id === 'c')?.children?.find((c) => c.kind === 'pairs');
    expect((links?.params?.rows as Record<string, unknown>).occupancy).toBeInstanceOf(Float32Array);
  });

  it('diff and patch back: each new param set and unset', () => {
    const b = createSeqSpec();
    b.data.sequence('s', { alphabet: 'protein', letters: 'MKTA' });
    const before = b.build();
    const after = moleculesSpec();
    const options = { language: SeqLanguage, validate: SeqValidateOptions, maxOps: Number.POSITIVE_INFINITY };
    const applied = applyPatch(before, diff(before, after), options);
    expect(applied.issues).toEqual([]);
    expect(applied.ok && hashNode(applied.root)).toBe(hashNode(after));
    const back = applyPatch(after, diff(after, before), options);
    expect(back.ok && hashNode(back.root)).toBe(hashNode(before));
  });

  it('take molecules and units as targets, and a complex as a section (EC-48)', () => {
    const b = createSeqSpec();
    const m = b.data.molecule('m', { residues: ['NAG', 'NAG'] });
    let u = '';
    const c = b.data.complex('c', {}, (cx) => {
      u = cx.unit('A', { of: m });
    });
    b.data.ranges('r', { target: m, rows: { start: [0], end: [2] } });
    b.data.values('v', { target: u, values: [1, 2] });
    b.data.pairs('p', { a: m, b: u, rows: { a_start: [0], a_end: [1], b_start: [1], b_end: [2] } });
    b.views.view('v', (v) => v.section('s', { of: c }));
    expect(checkSeqSpec(b.build())).toEqual([]);
    // past a molecule's residues, past a unit's range
    b.data.ranges('far', { target: m, rows: { start: [0], end: [3] } });
    b.data.values('short', { target: u, values: [1] });
    expect(codes(b.build())).toEqual(['seq.lang.ranges.bounds', 'seq.lang.values.length']);
  });

  it('a source names a host’s ref beside a uri, or alone said once (EC-185); never beside a value', () => {
    const b = createSeqSpec();
    b.data.source('both', { format: 'mmcif', uri: 'x.cif', ref: 'm' });
    b.data.source('alone', { format: 'mmcif', ref: 'm' });
    b.data.source('inline', { format: 'json', value: [1], ref: 'm' });
    const issues = checkSeqSpec(b.build());
    expect(issues.map((i) => [i.code, i.severity, i.address])).toEqual([
      ['seq.lang.source.host_only', 'info', '/data/alone'],
      ['seq.lang.source.content', 'error', '/data/inline'],
    ]);
  });

  it('a values field needs a source', () => {
    const b = createSeqSpec();
    const s = b.data.sequence('s', { alphabet: 'dna', letters: 'ACGT' });
    b.data.values('v', { target: s, values: [1, 2, 3, 4], field: 'plddt' });
    expect(codes(b.build())).toEqual(['seq.lang.values.content']);
  });

  it('modified: alternatives at one position pass, a component over several positions stops at the next row', () => {
    const ok = (rows: Parameters<typeof modified>[0]) => {
      const b = createSeqSpec();
      b.data.sequence('s', { alphabet: 'protein', letters: 'MKTAYIAK', modified: modified(rows) });
      return codes(b.build());
    };
    expect(
      ok([
        { position: 2, component: 'ALA' },
        { position: 2, component: 'GLY' },
        { position: 7, component: 'MSE' },
      ]),
    ).toEqual([]);
    expect(
      ok([
        { position: 2, component: 'CRO', length: 3 },
        { position: 5, component: 'MSE' },
      ]),
    ).toEqual([]);
    expect(ok([{ position: 2, component: 'X', length: 0 }])).toEqual(['seq.lang.modified.length']);
    expect(ok([{ position: -1, component: 'MSE' }])).toEqual(['seq.lang.modified.bounds']);
  });

  it('a unit never wraps (EC-13); its positions are its own, from 0: its numbering labels them', () => {
    const b = createSeqSpec();
    const s = b.data.sequence('s', { alphabet: 'protein', letters: 'MKTA' });
    b.data.complex('c', {}, (c) =>
      c.unit('A', { of: s, range: [1, 3] }, (u) => u.numbering('auth', { method: 'labels', labels: ['1', '2'] })),
    );
    expect(codes(b.build())).toEqual([]);
    b.data.complex('d', {}, (c) =>
      c.unit('A', { of: s }, (u) => u.numbering('auth', { method: 'labels', labels: ['1', '2'] })),
    );
    b.data.complex('e', {}, (c) => c.unit('A', { of: s, range: [3, 1] }));
    expect(codes(b.build())).toEqual(['seq.lang.numbering.length', 'seq.lang.unit.range']);
  });

  it('observed runs are on a unit’s own positions, from 0, and the issue says so', () => {
    const b = createSeqSpec();
    const s = b.data.sequence('s', { alphabet: 'protein', letters: 'MKTAYIAK' });
    b.data.complex('c', {}, (c) => c.unit('A', { of: s, range: [2, 6], observed: { start: [0], end: [4] } }));
    expect(codes(b.build())).toEqual([]);
    b.data.complex('d', {}, (c) => c.unit('A', { of: s, range: [2, 6], observed: { start: [2], end: [6] } }));
    const [issue] = checkSeqSpec(b.build());
    expect(issue.code).toBe('seq.lang.unit.observed');
    expect(`${issue.message} ${issue.hint}`).not.toMatch(/range/);
    expect(issue.hint).toMatch(/own positions/);
  });

  it('links: an intra-unit link and several rows at one site (glycoforms) pass; a link under data names no unit', () => {
    const b = createSeqSpec();
    const s = b.data.sequence('s', { alphabet: 'protein', letters: 'MKTACYAK' });
    const g = b.data.molecule('g', { residues: ['NAG'] });
    const h = b.data.molecule('h', { residues: ['MAN'] });
    b.data.complex('c', {}, (c) => {
      c.unit('A', { of: s });
      c.unit('G', { of: g });
      c.unit('H', { of: h });
      c.links('l', {
        rows: {
          a_unit: ['A', 'A', 'A'],
          a_start: [1, 1, 3],
          a_end: [2, 2, 4],
          b_unit: ['G', 'H', 'A'],
          b_start: [0, 0, 5],
          b_end: [1, 1, 6],
          occupancy: [0.7, 0.3, Number.NaN],
        },
      });
    });
    expect(checkSeqSpec(b.build())).toEqual([]);
    b.data.pairs('p', { a: s, rows: { a_unit: ['A'], a_start: [0], a_end: [1], b_start: [2], b_end: [3] } });
    b.data.pairs('q', { rows: { a_start: [0], a_end: [1], b_start: [2], b_end: [3] } } as never);
    expect(codes(b.build())).toEqual(['seq.lang.link.unit', 'seq.lang.link.target']);
  });

  it('a ranges layer takes pairs (each end a mark); a form of atoms needs atoms, a unit read through its molecule', () => {
    const b = createSeqSpec();
    const m = b.data.molecule('m', { residues: ['ATP'] });
    let u = '';
    b.data.complex('c', {}, (c) => {
      u = c.unit('L', { of: m });
    });
    b.views.view('v', (v) => v.section('s', {}, (sec) => sec.block('b', { of: u, form: 'atoms' })));
    expect(codes(b.build())).toEqual(['seq.lang.molecule.form']);
  });

  it('a link’s end is one position, under a complex or a molecule; pairs under data keep their ranges', () => {
    const b = createSeqSpec();
    const s = b.data.sequence('s', { alphabet: 'protein', letters: 'MKTAYIAK' });
    b.data.pairs('p', { a: s, rows: { a_start: [0], a_end: [3], b_start: [4], b_end: [6] } });
    expect(codes(b.build())).toEqual([]);
    b.data.complex('c', {}, (c) => {
      c.unit('A', { of: s });
      c.links('l', { rows: { a_unit: ['A'], a_start: [0], a_end: [1], b_unit: ['A'], b_start: [4], b_end: [6] } });
    });
    b.data.molecule('m', { residues: ['NAG', 'NAG', 'NAG'] }, (m) =>
      m.links('l', { rows: { a_start: [0], a_end: [2], b_start: [2], b_end: [3] } }),
    );
    expect(checkSeqSpec(b.build()).map((i) => [i.code, i.address, i.path?.slice(-2).join('/')])).toEqual([
      ['seq.lang.link.position', '/data/c/l', 'b_end/0'],
      ['seq.lang.link.position', '/data/m/l', 'a_end/0'],
    ]);
  });

  it('a section’s range is bounded on a unit’s own positions and a molecule’s residues, as on a sequence', () => {
    const b = createSeqSpec();
    const s = b.data.sequence('s', { alphabet: 'protein', letters: 'MKTAYIAK' });
    const m = b.data.molecule('m', { residues: ['NAG', 'NAG'] });
    let u = '';
    b.data.complex('c', {}, (c) => {
      u = c.unit('A', { of: s, range: [2, 4] });
    });
    b.views.view('v', (v) => {
      v.section('unit', { of: u, range: [0, 2] });
      v.section('molecule', { of: m, range: [0, 2] });
    });
    expect(codes(b.build())).toEqual([]);
    b.views.view('w', (v) => {
      v.section('unit', { of: u, range: [0, 3] });
      v.section('molecule', { of: m, range: [0, 6] });
    });
    expect(checkSeqSpec(b.build()).map((i) => [i.code, i.address])).toEqual([
      ['seq.lang.section.range', '/views/w/unit'],
      ['seq.lang.section.range', '/views/w/molecule'],
    ]);
  });

  it('block.form’s default is one rule, in the description and the _spec', () => {
    const rule = '`whole` for a molecule or a unit of one, `residues` for a sequence or a unit of one';
    expect(SeqLanguage.kinds.block.params.form.doc).toContain(rule);
    const spec = readFileSync(new URL('../../../_spec/README.md', import.meta.url), 'utf8').replace(/\s+/g, ' ');
    expect(spec.split(rule).length - 1).toBe(2);
  });

  it('modified tables of any rows: in order and in bounds exactly when no issue is said (property)', () => {
    const next = rng(7);
    for (let run = 0; run < 200; run++) {
      const length = 1 + Math.floor(next() * 12);
      const rows = Array.from({ length: Math.floor(next() * 5) }, () => ({
        position: Math.floor(next() * (length + 2)) - 1,
        component: 'MSE',
        length: 1 + Math.floor(next() * 3),
      }));
      const b = createSeqSpec();
      b.data.sequence('s', { alphabet: 'protein', letters: 'A'.repeat(length), modified: modified(rows) });
      const spec = b.build();
      let fits = true;
      rows.forEach((r, i) => {
        const next = rows.slice(i + 1).find((x) => x.position !== r.position);
        if (i > 0 && r.position < rows[i - 1].position) fits = false;
        if (r.position < 0 || r.position + r.length > length) fits = false;
        if (next && r.position + r.length > next.position) fits = false;
      });
      expect(codes(spec).length === 0).toBe(fits);
      // printing and reading back keeps the table
      const back = parseSexpr(printSexpr(spec));
      expect(hashNode(normalize(SeqLanguage, back.node as Node, SeqValidateOptions))).toBe(hashNode(spec));
    }
  });
});
