/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Labels as written (seqspec:seqql#labels-and-indices): a label that reads as a number keeps its written form (CM-01).
import { checkText, printQuery } from '@molstar/seq-core/query';
import { describe, expect, it } from 'vitest';
import { SeqSets } from '../../src/model/index.js';
import { compileQuery, createQueryEnv, SeqSelection } from '../../src/model/query/index.js';

const ENV = createQueryEnv();
// IMGT's long CDR3: 111.1 … 111.10, where 111.10 is the tenth insertion, not 111.1.
const labels = [110, '111.1', '111.2', '111.3', '111.4', '111.5', '111.6', '111.7', '111.8', '111.9', '111.10', 112];
const { set } = SeqSets.build({
  molecules: [
    {
      key: 'h',
      residues: 'ACDEFGHIKLMN',
      ownUnit: true,
      numberings: [{ name: 'imgt', method: 'labels', labels, isDefault: true }],
    },
  ],
});
const run = (text: string) => {
  const c = checkText(text, ENV.language);
  const r = compileQuery(c.expanded as never, ENV).run({ set });
  return { printed: printQuery(c.expr as never), summary: SeqSelection.summary(r.selection, { indices: true }) };
};

describe('a label written as a number', () => {
  it('selects the label written, quoted or not, and prints so the summary copies back', () => {
    const tenth = '1 group, 1 position: h/unit 111.10 @10:11';
    expect(run('(range h 111.10)')).toEqual({ printed: '(range h "111.10")', summary: tenth });
    expect(run('(range h "111.10")').summary).toBe(tenth);
    expect(run('(sel "h[111.10]")').summary).toBe(tenth);
    expect(run('(range h 111.1)').summary).toBe('1 group, 1 position: h/unit 111.1 @1:2');
  });

  it('keeps a plain integer a number', () => {
    expect(run('(range h 112)').printed).toBe('(range h 112)');
  });

  it('only a label keeps its written form: a text argument written as a number is refused, whatever its form (LA-23)', () => {
    const codes = (text: string) => checkText(text, ENV.language).issues.map((i) => i.code);
    for (const text of ['(ids f 1.0)', '(ids f 10)', '(positions :in h :where (matches (label) 1.50))'])
      expect(codes(text), text).toEqual(['seqql.type']);
  });

  it("a suffix no label holds is placed by its order, and the set's strings are left as they were (CM-11)", () => {
    const strings = set.strings.size;
    // `111.11` is no label: placed after 111.10, so the range from it starts at 112.
    expect(run('(range h 111.11-112)').summary).toBe('1 group, 1 position: h/unit 112 @11:12');
    expect(run('(range h 110-111.0)').summary).toBe('1 group, 1 position: h/unit 110 @0:1');
    expect([set.strings.find('.11'), set.strings.find('.0'), set.strings.size]).toEqual([-1, -1, strings]);
  });
});

describe('a numbering named by its namespace (CM-02)', () => {
  const molecule = (numberings: { name: string; first: number }[]) =>
    SeqSets.build({
      molecules: [
        {
          key: 'p',
          residues: 'MEEPQSDPSVEPPLSQ',
          ownUnit: true,
          numberings: numberings.map((n) => ({ name: n.name, method: 'offset' as const, first: n.first })),
        },
      ],
    }).set;
  const on = (s: ReturnType<typeof molecule>, text: string) => {
    const c = checkText(text, ENV.language);
    const r = compileQuery(c.expanded as never, ENV).run({ set: s });
    return r.issues.length > 0 ? r.issues.map((i) => `${i.code}: ${i.message}`) : SeqSelection.summary(r.selection);
  };

  it('`uniprot` is the one `uniprot:<acc>` a loader made, in range and sel, and the full name is written in sel', () => {
    const s = molecule([{ name: 'uniprot:P12345', first: 100 }]);
    // Indices 2 to 4, summarized in the default (ordinal) numbering.
    const want = '1 group, 3 positions: p/unit 3-5';
    expect(on(s, '(range p 102-104 :numbering uniprot)')).toBe(want);
    expect(on(s, '(sel "p[uniprot: 102-104]")')).toBe(want);
    expect(on(s, '(sel "p[uniprot:P12345: 102-104]")')).toBe(want);
    expect(on(s, '(range p 102-104 :numbering uniprot:P12345)')).toBe(want);
  });

  it('with several isoforms it is an error that lists them; the full name picks one', () => {
    const s = molecule([
      { name: 'uniprot:P12345', first: 100 },
      { name: 'uniprot:P12345-2', first: 200 },
    ]);
    expect(on(s, '(range p 102-104 :numbering uniprot)')).toEqual([
      expect.stringMatching(
        /^seqql\.numbering\.ambiguous: p\/unit has 2 numberings under "uniprot": uniprot:P12345, uniprot:P12345-2; name one/,
      ),
    ]);
    expect(on(s, '(sel "p[uniprot:P12345-2: 205-207]")')).toBe('1 group, 3 positions: p/unit 6-8');
  });
});

describe('(label) against constants (CM-04)', () => {
  // Labels with suffixes, an offset numbering, names, and a molecule with none (labels are index + 1).
  const { set: s } = SeqSets.build({
    molecules: [
      {
        key: 'h',
        residues: 'ACDEFGHIKLMN',
        ownUnit: true,
        numberings: [
          { name: 'imgt', method: 'labels', labels, isDefault: true },
          { name: 'shift', method: 'offset', first: -3 },
          { name: 'atoms', method: 'names', names: ['CA', 'CB', 'N', 'CA', 'C', 'O', 'CA', 'N', 'CB', 'C', 'O', 'N'] },
        ],
      },
      { key: 'plain', residues: 'ACDEFGHIK', ownUnit: true },
    ],
  });
  const picked = (text: string) => {
    const c = checkText(text, ENV.language);
    const r = compileQuery(c.expanded as never, ENV).run({ set: s });
    return SeqSelection.summary(r.selection, { indices: true });
  };
  /** The same test through `matches`, which reads each label's text: the probe must agree with it. */
  const viaText = (q: string, values: string[]) => {
    const re = values.map((v) => v.replace(/[.]/g, '\\\\.')).join('|');
    return picked(`(positions :in (or h plain) :where (matches (label${q}) "^(${re})$" :case true))`);
  };
  const cases: [string, string[]][] = [
    ['', ['111.10']],
    ['', ['111.1', '112', '110']],
    ['', ['0111.1', '111.1x', '', 'abc']],
    [' :numbering shift', ['-3', '0', '1', '8']],
    [' :numbering atoms', ['CA', 'N']],
    [' :numbering nope', ['1']],
  ];
  it("= and in agree with the label's text, numbering by numbering", () => {
    for (const [q, values] of cases) {
      const eq = picked(
        `(positions :in (or h plain) :where (in (label${q}) ${values.map((v) => JSON.stringify(v)).join(' ')}))`,
      );
      expect(eq, `${q} ${values}`).toBe(viaText(q, values));
      const one = values[0] as string;
      expect(picked(`(positions :in (or h plain) :where (= (label${q}) ${JSON.stringify(one)}))`)).toBe(
        viaText(q, [one]),
      );
      expect(picked(`(positions :in (or h plain) :where (= ${JSON.stringify(one)} (label${q})))`)).toBe(
        viaText(q, [one]),
      );
    }
  });

  it('!= is its complement where the label is there', () => {
    const not = picked('(positions :in h :where (not (matches (label) "^111\\\\.10$" :case true)))');
    expect(picked('(positions :in h :where (!= (label) "111.10"))')).toBe(not);
    expect(not).toMatch(/11 positions/);
  });
});

describe('on a cyclic system (CM-10)', () => {
  const { set: ring } = SeqSets.build({
    molecules: [
      { key: 'ring', residues: 'ACGT'.repeat(10), ownUnit: true, cyclic: true },
      { key: 'line', residues: 'ACGT'.repeat(10), ownUnit: true },
    ],
  });
  const on = (text: string) => {
    const c = checkText(text, ENV.language);
    const r = compileQuery(c.expanded as never, ENV).run({ set: ring });
    return r.issues.length > 0
      ? r.issues.map((i) => `${i.code}: ${i.message}`)
      : SeqSelection.summary(r.selection, { indices: true });
  };

  it('a range, a slice and an expand wrap the origin', () => {
    const across = '1 group, 6 positions: ring/unit 1-3 @0:3, 38-40 @37:40';
    expect(on('(range ring 38-3)')).toBe(across);
    expect(on('(slice ring 37 3)')).toBe(across);
    expect(on('(expand (slice ring 0 2) 2)')).toBe('1 group, 6 positions: ring/unit 1-4 @0:4, 39-40 @38:40');
    expect(on('(expand (slice ring 0 2) 30)')).toBe('1 group, 40 positions: ring/unit 1-40 @0:40');
  });

  it("the sel shorthand's slice wraps the origin too; only equal ends are refused (LA-07)", () => {
    expect(on('(sel "ring@37:3")')).toBe('1 group, 6 positions: ring/unit 1-3 @0:3, 38-40 @37:40');
    expect(checkText('(sel "ring@3:3")', ENV.language).issues.map((i) => `${i.code}: ${i.message}`)).toEqual([
      expect.stringMatching(/^seqql\.sel\.syntax: .*equal/),
    ]);
  });

  it("a start past a ring's end is outside it, not a start after its stop (LA-21)", () => {
    expect(on('(slice ring 45 3)')).toEqual([
      expect.stringMatching(/^seqql\.slice\.range: \(slice ring\/unit 45 3\) is outside it/),
    ]);
  });

  it('on a linear system a start after the stop is an error that says to write the start first', () => {
    expect(on('(slice line 37 3)')).toEqual([
      expect.stringMatching(
        /^seqql\.slice\.range: .*starts after its stop, which wraps only a cyclic sequence's origin.*; write the start first/,
      ),
    ]);
  });
});
