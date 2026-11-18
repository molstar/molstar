/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readdirSync, readFileSync } from 'node:fs';
import { PCG } from '@molstar/core/data/util/hash-functions';
import { describe, expect, it } from 'vitest';
import { addressFor, resolveAddressIn } from '../address.js';
import { AddressIssues } from '../address-issues.js';
import { parseLink, printLink } from '../links.js';
import { location, locationsEqual } from '../location.js';
import { Model, Spec, Systems } from './fixture.js';

// Links in Mol*'s grammar, both ways; the addressing catalog against the codes the sources use; and the way back from
// a Location to an address as a seeded property: whatever the numbering asked, the address reads back to it.

describe('links', () => {
  it('print and parse both ways, an address percent-encoded as a whole', () => {
    const address = '/data/p53/(item :label "DNA binding & more" :id [a=1 #2])';
    const text = printLink({
      actions: [
        { action: 'focus', address: '/refs/dbd' },
        { action: 'highlight', address },
      ],
      channel: 'link',
      other: [['structure-focus', 'A:94']],
    });
    expect(text).toBe(
      '!structure-focus=A%3A94&seq-focus=%2Frefs%2Fdbd&seq-highlight=%2Fdata%2Fp53%2F(item%20%3Alabel%20%22DNA%20binding%20%26%20more%22%20%3Aid%20%5Ba%3D1%20%232%5D)&seq-channel=link',
    );
    expect(parseLink(text)).toEqual({
      link: {
        actions: [
          { action: 'focus', address: '/refs/dbd' },
          { action: 'highlight', address },
        ],
        channel: 'link',
        other: [['structure-focus', 'A:94']],
      },
      issues: [],
    });
  });

  it('read the short form a writer types, and another prefix', () => {
    expect(parseLink('!seq-focus=/refs/dbd')?.link?.actions).toEqual([{ action: 'focus', address: '/refs/dbd' }]);
    expect(parseLink('!x-select=/a', { prefix: 'x-' })?.link?.actions).toEqual([{ action: 'select', address: '/a' }]);
    expect(parseLink('https://example.org')).toBeUndefined();
  });

  it('leave out an unknown action with a warning, and refuse a bad address, bad encoding and a long link', () => {
    const unknown = parseLink('!seq-explode=/a&seq-select=/b')!;
    expect(unknown.link?.actions).toEqual([{ action: 'select', address: '/b' }]);
    expect(unknown.issues.map((i) => [i.code, i.offset])).toEqual([['seq.address.link.action', 1]]);
    expect(parseLink('!seq-focus=data/p53')!.link).toBeUndefined();
    expect(parseLink('!seq-focus=%E0%A4%A')!.issues.map((i) => i.code)).toEqual(['seq.address.link.syntax']);
    expect(parseLink(`!seq-focus=/${'a'.repeat(100)}`, { maxLength: 50 })!.issues.map((i) => i.code)).toEqual([
      'seq.address.link.limit',
    ]);
    expect(parseLink('!seq-focus=/a&seq-focus=/b', { maxActions: 1 })!.issues.map((i) => i.code)).toEqual([
      'seq.address.link.limit',
    ]);
  });
});

describe('seq.address: issue codes', () => {
  it('every code a source uses is cataloged, and every cataloged code is used', () => {
    const dir = new URL('..', import.meta.url);
    const source = readdirSync(dir)
      .filter((f) => f.endsWith('.ts') && f !== 'address-issues.ts')
      .map((f) => readFileSync(new URL(f, dir), 'utf8'))
      .join('\n');
    const used = new Set(source.match(/'seq\.address\.[a-z0-9_.]+'/g)?.map((s) => s.slice(1, -1)));
    const cataloged = Object.keys(AddressIssues.codes);
    for (const c of used) expect(cataloged, c).toContain(c);
    for (const c of cataloged) expect(used.has(c), c).toBe(true);
  });
});

describe('the way back, as a property', () => {
  const cases: [system: keyof typeof Systems, numberings: string[]][] = [
    ['p53', ['ordinal', 'index', 'uniprot', 'dup']],
    ['plasmid', ['ordinal', 'gapped']],
    ['ab', ['ordinal', 'kabat', 'imgt', 'auth']],
    ['offset', ['ordinal', 'skip', 'dot']],
  ];

  it('reads back to the same Location, in any numbering asked', () => {
    const rng = new PCG(7);
    let checked = 0;
    for (let n = 0; n < 200; n++) {
      const [key, numberings] = cases[rng.int() % cases.length];
      const s = Systems[key];
      const ranges: number[] = [];
      for (let i = rng.int() % 3; i >= 0; i--) {
        const a = rng.int() % s.length;
        const b = rng.int() % s.length;
        ranges.push(Math.min(a, b), Math.max(a, b) + 1);
      }
      const target = location({ parts: [{ system: s.address, ranges }] });
      const numbering = numberings[rng.int() % numberings.length];
      const a = addressFor(Spec, Model, { type: 'location', location: target }, { numbering, prefer: ['ids'] });
      expect(a, `${s.address} ${ranges} in ${numbering}`).toBeDefined();
      const back = resolveAddressIn(Spec, Model, a!).result;
      expect(back?.type === 'location' && locationsEqual(back.location, target), a).toBe(true);
      checked++;
    }
    expect(checked).toBe(200);
  });
});
