/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A seeded fuzz of 2,000 mutations over the acceptance cases of cases.ts (the heavy lane; the quick lane runs 100). Whatever a mutation does to a
// document (a key deleted, a value replaced by junk, children turned into an object, a kind swapped, a node
// duplicated), `validate`, `normalize`, `diff` and `applyPatch` report issues and never throw; and when both documents
// validate, `applyPatch(structuredClone(a), diff(a, b))` gives b.

import { describe, expect, it } from 'vitest';
import {
  applyPatch,
  canonical,
  diff,
  normalize,
  type Patch,
  type SeqViewSpecDoc,
  snapshotDoc,
  validate,
} from '../../src/format/index.js';
import { CASES, case17 } from './cases.js';
import { registries } from './load-fixtures.js';

type J = any;
/** A document as normalized text: keys sorted, and no difference between absent children and none. */
const same = (d: SeqViewSpecDoc) =>
  canonical(normalize(d, reg).root).replaceAll(',"children":[]', '').replaceAll('"children":[],', '');
const isError = (i: { severity: string }) => i.severity === 'error';
const reg = registries();

function random(seed: number) {
  let s = seed >>> 0;
  return () => {
    s = (Math.imul(s, 1664525) + 1013904223) >>> 0;
    return s / 4294967296;
  };
}

/** Every key of every object and array in a document, as (container, key). */
function slots(v: J, out: [J, string][] = []): [J, string][] {
  if (v !== null && typeof v === 'object') {
    for (const k of Object.keys(v)) {
      out.push([v, k]);
      slots(v[k], out);
    }
  }
  return out;
}

const JUNK = () => [
  null,
  0,
  -1,
  1.5,
  '',
  'x',
  'a b',
  true,
  [],
  {},
  [null],
  [{}],
  { a: 1 },
  { kind: 'track' },
  ['x', 1],
];
const KINDS = [
  'root',
  'data',
  'source',
  'sequence',
  'view',
  'section',
  'block',
  'group',
  'track',
  'layer',
  'focus',
  'widget',
];

/** One random mutation of `doc` in place; what it did, for the failure message. */
function mutate(doc: J, rnd: () => number): string {
  const all = slots(doc);
  const [parent, key] = all[Math.floor(rnd() * all.length)];
  const at = `${Array.isArray(parent) ? '[]' : '{}'}.${key}`;
  const r = rnd();
  if (r < 0.35) {
    const junk = JUNK();
    const v = junk[Math.floor(rnd() * junk.length)];
    parent[key] = v;
    return `${at} = ${JSON.stringify(v)}`;
  }
  if (r < 0.6) {
    if (Array.isArray(parent)) parent.splice(Number(key), 1);
    else delete parent[key];
    return `delete ${at}`;
  }
  if (r < 0.7 && Array.isArray(parent) && parent[Number(key)] !== null && typeof parent[Number(key)] === 'object') {
    parent.push(structuredClone(parent[Number(key)]));
    return `duplicate ${at}`;
  }
  if (r < 0.8) {
    parent[key] = [parent[key]];
    return `wrap ${at}`;
  }
  const v = parent[key];
  if (v !== null && typeof v === 'object' && typeof v.kind === 'string') {
    v.kind = KINDS[Math.floor(rnd() * KINDS.length)];
    return `kind of ${at} = ${v.kind}`;
  }
  if (typeof v === 'string') {
    parent[key] = v.length > 1 ? v.slice(1) : `${v}~`;
    return `edit ${at}`;
  }
  parent[key] = typeof v === 'number' ? -v - 1 : null;
  return `negate ${at}`;
}

/** `mutations` seeded mutations, cycling through the cases; the failures, and how many diffs applied, of valid pairs. */
function fuzz(mutations: number, seed: number) {
  // Case 17 small (8 samples over 20 sites): its 600,000 genotypes would make its mutations' clones, validations and
  // diffs most of the run, and its shapes (an inline matrix, states, `order: "data"`, rows keyed by state) stay in.
  const makes = Object.entries(CASES).map(([name, make]) => (name === 'case17' ? () => case17(8, 20) : make));
  const rnd = random(seed);
  const failures: string[] = [];
  let applied = 0;
  let both = 0;
  for (let n = 0; n < mutations; n++) {
    const a: SeqViewSpecDoc = makes[n % makes.length]();
    const b: SeqViewSpecDoc = structuredClone(a);
    const what = mutate(b, rnd);
    const run = (name: string, f: () => unknown) => {
      try {
        return f();
      } catch (e) {
        failures.push(
          `#${n} ${what}: ${name} threw ${e instanceof Error ? e.stack?.split('\n').slice(0, 3).join(' | ') : e}`,
        );
        return undefined;
      }
    };
    const found = run('validate(b)', () => validate(b, reg)) as ReturnType<typeof validate> | undefined;
    run('normalize(b)', () => normalize(b, reg));
    // Both valid: the diff must apply both ways (validated once each; a is a case, valid).
    const valid = !!found && !found.some(isError) && !validate(a, reg).some(isError);
    run('snapshotDoc(b)', () => {
      for (const s of b.snapshots ?? []) snapshotDoc(b, (s as J)?.ref, reg);
    });
    for (const [x, y, label] of [
      [a, b, 'a to b'],
      [b, a, 'b to a'],
    ] as const) {
      const p = run(`diff(${label})`, () => diff(x, y, reg)) as Patch | null | undefined;
      if (p === undefined) continue;
      const r = run(`applyPatch(${label})`, () => applyPatch(structuredClone(x), p ?? ({ ops: 0 } as never), reg)) as
        | ReturnType<typeof applyPatch>
        | undefined;
      if (!r || p === null) continue;
      applied++;
      if (valid) {
        both++;
        const errs = r.issues.filter((i) => i.severity === 'error');
        if (errs.length > 0)
          failures.push(`#${n} ${what}: ${label}: the diff does not apply: ${errs[0].code} ${errs[0].path}`);
        else if (same(r.doc) !== same(y)) {
          failures.push(`#${n} ${what}: ${label}: applyPatch(a, diff(a, b)) is not b`);
        }
      }
    }
    expect(found === undefined || Array.isArray(found)).toBe(true);
  }
  return { failures, applied, both };
}

describe('fuzz', () => {
  it('100 mutations of the acceptance cases never throw, and a diff of two valid documents applies', () => {
    const { failures, applied, both } = fuzz(100, 20261011);
    expect(failures).toEqual([]);
    expect(applied).toBeGreaterThan(50);
    expect(both).toBeGreaterThan(2);
  });

  it('2000 mutations of the acceptance cases never throw, and a diff of two valid documents applies', {
    tags: ['heavy'],
  }, () => {
    const { failures, applied, both } = fuzz(2000, 20261011);
    expect(failures).toEqual([]);
    expect(applied).toBeGreaterThan(1000);
    expect(both).toBeGreaterThan(50);
  });
});
