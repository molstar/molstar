/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Extensions of the query language (extension-openings.md; properties.md): a property by a function and by a values
// annotation's field, a form of a role, registered through the registry and composed by the context; one missing, one
// not namespaced, two under one name; a filter giving what its prefix does not hold refused by its role; a values
// channel pending while its tiles have not arrived, staging its system.

import type { Node } from '@molstar/core/repr/tree';
import { composeSeqContext, queryExtensionsOf } from '@molstar/seqx-core/context';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import {
  type PropertyColumn,
  type QueryExtensionEntry,
  type QueryFormEntry,
  type QueryPropertyEntry,
  rangesLocation,
} from '@molstar/seqx-core/query';
import { describe, expect, it } from 'vitest';
import { codes, ex, glycanSpec, loc, prop, specRun } from './structures.js';

const column = (n: number, v: (i: number) => number): PropertyColumn => ({ length: n, get: (i) => v(i) });

/** pLDDT by a function: ten times a residue's position on its system. */
const plddt: QueryPropertyEntry = {
  kind: 'query-property',
  name: 'af:plddt',
  aliases: ['af:confidence'],
  level: 'residue',
  type: 'number',
  read: (model, system) => {
    const s = model.system(system);
    return s ? column(s.length, (i) => i * 10) : undefined;
  },
};

/** A generator: its prefix's nodes, as they are. */
const same: QueryFormEntry = {
  kind: 'query-form',
  name: 'af:same',
  role: 'generator',
  params: [],
  result: 'nodes',
  resolve: (ctx) => ctx.prefix,
};

/** A filter that is not one: it gives the whole of `:of` whatever its prefix. */
const leaky: QueryFormEntry = {
  kind: 'query-form',
  name: 'af:leaky',
  role: 'filter',
  params: ['of'],
  result: 'location',
  resolve: (ctx) => ctx.evaluate(ctx.form.params!.of),
};

const ext = (...entries: QueryExtensionEntry[]) => composeSeqContext({ registry: entries }).query;

describe('query extensions', () => {
  const g = specRun(glycanSpec(), ext(plddt, same, leaky));
  const D = '/data/c/D';

  it('a property by a function: an elements field (quoted) and a filter’s, by its name or alias', () => {
    expect(loc(g.run(`${D}/(elements :"af:plddt" 30)`), D)).toEqual([3, 4]);
    const f = g.form(`${D}/(elements)`, { kind: 'filter', params: { test: ex('>=', prop('af:confidence'), 60) } });
    expect(loc(f, D)).toEqual([6, 8]);
  });

  it('a form by its role; a filter giving what its prefix does not hold refused', () => {
    expect(codes(g.run(`${D}/(af:same)`))).toEqual([]);
    expect(codes(g.run(`${D}/(range :from 1 :to 2)/(af:leaky :of /data/c/D/(all))`))).toEqual(['seq.query.form.type']);
  });

  it('missing: a namespaced form or property nobody registered', () => {
    const bare = specRun(glycanSpec());
    expect(codes(bare.run(`${D}/(af:same)`))).toEqual(['seq.query.extension.missing']);
    expect(
      codes(bare.form(`${D}/(elements)`, { kind: 'filter', params: { test: ex('>', prop('af:plddt'), 1) } })),
    ).toEqual(['seq.query.extension.missing']);
  });

  it('composition: a name not namespaced, a built-in shadowed, two entries under one name, every one said', () => {
    const bad = { ...plddt, name: 'plddt' };
    expect(() => queryExtensionsOf([bad])).toThrow(/seq.context.extension.namespace/);
    expect(() => queryExtensionsOf([plddt, { ...plddt }])).toThrow(/af:plddt is registered twice/);
    expect(() => queryExtensionsOf([{ ...same, name: 'elements' }, bad])).toThrow(/elements.*plddt/);
    // the same entry twice is one
    expect(queryExtensionsOf([plddt, plddt]).properties.defs).toEqual(queryExtensionsOf([plddt]).properties.defs);
  });
});

/** A sequence with a values channel read from a source, and a store whose one tile is loaded when told. */
function valuesSetup(loaded: { value: boolean }, extensions?: QueryExtensionEntry[]) {
  const b = createSeqSpec();
  const host = b.data.source('host', { format: 'mmcif', uri: 'host.cif' });
  const fc = b.data.sequence('fc', { alphabet: 'protein', letters: 'ACDEFG' });
  b.data.values('plddt', { target: fc, source: host, field: 'plddt' });
  const spec: Node = b.build();
  const data = {
    coverage: () => ({ systems: [{ name: 'fc', length: 6 }], shapes: [], levels: [1] }),
    keysOf: () => ['t0'],
    get: () =>
      loaded.value
        ? {
            status: 'loaded',
            tables: { values: { form: 'positions', start: 0, value: Float32Array.of(50, 90, 20, 80, 95, 10) } },
          }
        : { status: 'pending' },
  };
  return specRun(spec, extensions && ext(...extensions), { data: data as never });
}

describe('values channels', () => {
  const test = (t: Node) => ({ kind: 'filter', params: { test: t } });

  it('pending while the tile has not arrived: the system staged, then read, as a query over data read by region (EC-229)', () => {
    const loaded = { value: false };
    const v = valuesSetup(loaded);
    const high = test(ex('>=', ex('value', '/data/plddt'), 70));
    const waiting = v.form('/data/fc/(elements)', high);
    expect([waiting.complete, codes(waiting)]).toEqual([false, ['seq.query.property.pending']]);
    expect(waiting.pending).toEqual(rangesLocation('/data/fc', [0, 6]));
    loaded.value = true;
    expect(loc(valuesSetup(loaded).form('/data/fc/(elements)', high), '/data/fc')).toEqual([1, 2, 3, 5]);
  });

  it('an extension’s property by a values annotation’s field, by declaration alone', () => {
    const byValues: QueryPropertyEntry = {
      kind: 'query-property',
      name: 'af:score',
      level: 'residue',
      type: 'number',
      read: { values: '/data/plddt', field: 'plddt' },
    };
    const v = valuesSetup({ value: true }, [byValues]);
    expect(loc(v.form('/data/fc/(elements)', test(ex('<', prop('af:score'), 30))), '/data/fc')).toEqual([2, 3, 5, 6]);
  });

  it('b-factor without a channel: lacking, a warning, nothing (MC writes the channel)', () => {
    const g = specRun(glycanSpec());
    const r = g.form('/data/c/D/(elements :atom C1)', test(ex('>', prop('b-factor'), 60)));
    expect([loc(r, '/data/c/D/(atoms)'), codes(r)]).toEqual([[], ['seq.query.property.absent']]);
  });
});
