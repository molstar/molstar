/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { parseJson, printJson, toJson } from '@molstar/core/repr/codecs/json';
import { parseSexpr, printSexpr } from '@molstar/core/repr/codecs/sexpr';
import { languageToMarkdown, normalize, validateAndNormalize } from '@molstar/core/repr/schema';
import { hashNode, type Node } from '@molstar/core/repr/tree';
import { describe, expect, it } from 'vitest';
import { LanguageVersion } from '../../index.js';
import { createSeqSpec } from '../builder.js';
import { checkSeqSpec } from '../check.js';
import { LangIssues } from '../issues.js';
import { SeqLanguage, SeqLanguageVersion } from '../language.js';
import { EXONS, P53_LETTERS, p53Spec } from './p53.js';

describe('the language', () => {
  it('is the package version, and documents every kind', () => {
    expect(SeqLanguageVersion).toBe(LanguageVersion);
    const md = languageToMarkdown(SeqLanguage);
    for (const kind of Object.keys(SeqLanguage.kinds)) expect(md).toContain(kind);
  });
});

describe('the p53 example', () => {
  const spec = p53Spec();

  it('is valid, with handles that are addresses', () => {
    expect(checkSeqSpec(spec)).toEqual([]);
    const b = createSeqSpec();
    const p53 = b.data.sequence('p53', { alphabet: 'protein', letters: 'MEE' });
    expect(p53).toBe('/data/p53');
    const protein = b.views.view('main', (v) => {
      expect(v.section('protein', { of: p53 }).address).toBe('/views/main/protein');
    });
    expect(protein).toBe('/views/main');
  });

  it('prints as the description shows it', () => {
    const text = printSexpr(spec, { head: 64 });
    expect(text).toContain('(spec :version 0.0.0 :title p53');
    expect(text).toContain(`(sequence p53 :alphabet protein :letters ${P53_LETTERS}`);
    expect(text).toContain('(map p53-on-gene');
    expect(text).toContain(':follow /views/main/protein');
    expect(text).toContain('(layer :type ranges :data /data/domains)');
  });

  it('encodes as JSON both ways and back', () => {
    // declared columns plain
    const plain = printJson(spec, { language: SeqLanguage });
    expect(JSON.parse(plain).children[0].children[0].params.letters).toBe(P53_LETTERS);
    expect(plain).not.toContain('$column');
    const a = parseJson(plain);
    expect(a.issues).toEqual([]);
    const read = validateAndNormalize(SeqLanguage, a.node as Node);
    expect(read.issues).toEqual([]);
    expect(hashNode(read.node as Node)).toBe(hashNode(spec));
    // every column tagged
    const tagged = printJson(spec);
    expect(tagged).toContain('"$column":"int32"');
    const b = parseJson(tagged);
    expect(b.issues).toEqual([]);
    expect(hashNode(normalize(SeqLanguage, b.node as Node))).toBe(hashNode(spec));
    expect(toJson(spec, { language: SeqLanguage })).toEqual(JSON.parse(plain));
  });

  it('prints and reads back as S-expressions', () => {
    const back = parseSexpr(printSexpr(spec, { head: 64 }));
    expect(back.issues).toEqual([]);
    expect(hashNode(normalize(SeqLanguage, back.node as Node))).toBe(hashNode(spec));
  });

  it('writes the declared columns as plain lists', () => {
    const json = JSON.parse(printJson(spec, { language: SeqLanguage }));
    expect(json.children[0].children[3].params.exons.start).toEqual(EXONS.start);
  });
});

/** One spec for each code of the catalog: valid for the description, wrong in the language's own way. */
const Misuses: Record<keyof typeof LangIssues.codes, () => Node> = {
  'seq.lang.version.unsupported': () => createSeqSpec({ version: '1.0.0' }).build(),
  'seq.lang.container.id': () => ({
    kind: 'spec',
    params: { version: '0.0.0' },
    children: [{ kind: 'data', id: 'stuff' }],
  }),
  'seq.lang.source.content': () => {
    const b = createSeqSpec();
    b.data.source('s', { format: 'fasta' });
    return b.build();
  },
  'seq.lang.sequence.content': () => {
    const b = createSeqSpec();
    b.data.sequence('s', { alphabet: 'dna', letters: 'ACGT', source: '/data/x' });
    return b.build();
  },
  'seq.lang.sequence.letters': () => {
    const b = createSeqSpec();
    b.data.sequence('s', { alphabet: 'dna', letters: 'ACGU' });
    return b.build();
  },
  'seq.lang.ranges.content': () => {
    const b = createSeqSpec();
    const s = b.data.sequence('s', { alphabet: 'dna', letters: 'ACGT' });
    b.data.ranges('r', { target: s });
    return b.build();
  },
  'seq.lang.ranges.interval': () => {
    const b = createSeqSpec();
    const s = b.data.sequence('s', { alphabet: 'dna', letters: 'ACGT' });
    b.data.ranges('r', { target: s, rows: { start: [3], end: [1] } });
    return b.build();
  },
  'seq.lang.ranges.bounds': () => {
    const b = createSeqSpec();
    const s = b.data.sequence('s', { alphabet: 'dna', letters: 'ACGT' });
    b.data.ranges('r', { target: s, rows: { start: [0], end: [5] } });
    return b.build();
  },
  'seq.lang.values.content': () => {
    const b = createSeqSpec();
    const s = b.data.sequence('s', { alphabet: 'dna', letters: 'ACGT' });
    b.data.values('v', { target: s });
    return b.build();
  },
  'seq.lang.values.length': () => {
    const b = createSeqSpec();
    const s = b.data.sequence('s', { alphabet: 'dna', letters: 'ACGT' });
    b.data.values('v', { target: s, values: [1, 2, 3] });
    return b.build();
  },
  'seq.lang.map.params': () => {
    const b = createSeqSpec();
    const s = b.data.sequence('s', { alphabet: 'dna', letters: 'ACGT' });
    b.data.map('m', { from: s, to: s, method: 'offset' });
    return b.build();
  },
  'seq.lang.map.alphabet': () => {
    const b = createSeqSpec();
    const s = b.data.sequence('s', { alphabet: 'dna', letters: 'ACGT' });
    b.data.map('m', { from: s, to: s, method: 'codon', exons: EXONS, strand: '+' });
    return b.build();
  },
  'seq.lang.ref.target': () => {
    const b = createSeqSpec();
    b.data.ranges('r', { target: '/data/missing', rows: { start: [0], end: [1] } });
    return b.build();
  },
  'seq.lang.layer.data': () => {
    const b = createSeqSpec();
    const s = b.data.sequence('s', { alphabet: 'dna', letters: 'ACGT' });
    b.views.view('v', (v) => {
      v.section('sec', { of: s }).track('t').layer({ type: 'ranges', data: s });
    });
    return b.build();
  },
  'seq.lang.section.follow': () => {
    const b = createSeqSpec();
    const s = b.data.sequence('s', { alphabet: 'dna', letters: 'ACGT' });
    b.views.view('v', (v) => {
      v.section('a', { of: s, follow: '/views/v/b' });
      v.section('b', { of: s, follow: '/views/v/a' });
    });
    return b.build();
  },
};

describe('the misuses', () => {
  it('cover the catalog', () => {
    expect(Object.keys(Misuses).sort()).toEqual(Object.keys(LangIssues.codes).sort());
  });

  for (const [code, make] of Object.entries(Misuses)) {
    it(`${code} is reported, with its severity and a place`, () => {
      const issues = checkSeqSpec(make()).filter((i) => i.code === code);
      expect(issues.length).toBeGreaterThan(0);
      const [issue] = issues;
      expect(issue.severity).toBe(LangIssues.codes[code as keyof typeof LangIssues.codes].severity);
      expect(issue.path).toBeDefined();
      expect(issue.hint).toBeDefined();
    });
  }

  it('name the node by its address', () => {
    const [issue] = checkSeqSpec(Misuses['seq.lang.sequence.letters']());
    expect(issue.address).toBe('/data/s');
    expect(issue.path).toEqual(['children', 0, 'children', 0, 'params', 'letters']);
  });

  it('leave a valid spec without issues, whatever the references', () => {
    const b = createSeqSpec();
    const s = b.data.sequence('s', { alphabet: 'protein', letters: 'MEEP' });
    b.data.values('v', { target: s, values: [1, 2, Number.NaN, 4] });
    expect(b.check()).toEqual([]);
  });

  it('give the description’s codes first, and the language’s on a valid tree only', () => {
    const bad = createSeqSpec();
    bad.data.sequence('s', { alphabet: 'cat' } as never);
    expect(bad.check().map((i) => i.code)).toEqual(['repr.schema.param.enum']);
    const missing = createSeqSpec();
    missing.data.sequence('s', {} as never);
    expect(missing.check().map((i) => i.code)).toEqual(['repr.schema.param.missing']);
  });

  it('are refused where they are written, when the types do not reach', () => {
    const b = createSeqSpec();
    b.data.sequence('s', { alphabet: 'dna', letters: 'A' });
    expect(() => b.data.sequence('s', { alphabet: 'dna', letters: 'A' })).toThrow(/twice/);
    expect(() => b.data.sequence('a b', { alphabet: 'dna', letters: 'A' })).toThrow(/not an id/);
  });
});
