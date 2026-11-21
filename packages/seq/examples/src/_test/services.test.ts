/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { existsSync, readFileSync } from 'node:fs';
import { resolveSpec } from '@molstar/seqx-core/resolve';
import { describe, expect, it } from 'vitest';
import { PackageFolder } from '../node-assets.js';
import { AlignedSegments, EntityLength, EntityStart } from '../p53/structure.js';
import { build as composed } from '../services/composed.js';
import { Recorded, ServiceHosts, Services, serviceOf, serviceUrl } from '../services/endpoints.js';
import { build as predictions, labScores } from '../services/predictions.js';
import { build as addressing } from '../services/addressing.js';
import { chainsOf, parseSifts, siftsSegments } from '../services/sifts.js';

// The services' fixtures and the SIFTS reading (examples.md, Showcases): what is recorded is where the locations say, the
// answers parse into the maps the examples' own hand-written p53 structure view holds, and the specs resolve.

const read = (file: string) => readFileSync(`${PackageFolder}${file}`, 'utf8');
const sifts = (entry: string) => parseSifts(read(Recorded.get(serviceUrl.sifts(entry))!));

describe('the endpoints', () => {
  it('every recorded location is a service’s, on a host its policy lists, and its file is there', () => {
    for (const [url, file] of Recorded) {
      expect(serviceOf(url), url).toBeDefined();
      expect(ServiceHosts, url).toContain(new URL(url).host);
      expect(existsSync(`${PackageFolder}${file}`), file).toBe(true);
    }
    expect(Services.map((s) => s.id)).toEqual(['uniprot', 'pdbe', 'alphafold']);
    expect(serviceOf('https://example.org/x')).toBeUndefined();
  });
});

describe('SIFTS', () => {
  it('2OCJ’s four chains map P04637 94–312 onto the entity’s 1–219: the segments the hand-written p53 view holds', () => {
    const m = sifts('2OCJ');
    expect(chainsOf(m, 'P04637')).toEqual(['A', 'B', 'C', 'D']);
    for (const s of siftsSegments(m)) {
      expect([s.from, s.to, s.length]).toEqual([EntityStart, AlignedSegments.to[0], EntityLength]);
    }
    expect(m[0].authStart).toBeNull();
  });

  it('1FIN: CDK2 whole on chains A and C, the cyclin A construct 173–432 on B and D (the entity’s 1–260)', () => {
    const s = siftsSegments(sifts('1FIN'));
    const by = (acc: string) => s.filter((x) => x.accession === acc);
    expect(by('P24941').map((x) => [x.chain, x.from, x.to, x.length])).toEqual([
      ['A', 0, 0, 298],
      ['C', 0, 0, 298],
    ]);
    expect(by('P20248').map((x) => [x.chain, x.from, x.to, x.length])).toEqual([
      ['B', 172, 0, 260],
      ['D', 172, 0, 260],
    ]);
  });
});

describe('the specs', () => {
  it.each([
    ['composed', composed],
    ['predictions', predictions],
    ['addressing', addressing],
  ] as const)('%s: resolves with and without SIFTS’s answer, the maps one a chain, ranked first', (_, build) => {
    for (const answer of [undefined, sifts('2OCJ')]) {
      const spec = build(answer ? { sifts: answer } : {});
      const model = resolveSpec(spec);
      // a source's data is not held here: only what the spec says itself is checked
      expect(model.issues.filter((i) => i.severity === 'error').map((i) => i.code)).toEqual([]);
      const maps = [...model.describedMaps!.values()];
      expect(maps.map((d) => d.to)).toEqual(
        answer ? ['A', 'B', 'C', 'D'].map((c) => `/data/2ocj/(node :key complex)/${c}`) : [],
      );
      expect(maps.every((d) => (d.node.params as { priority?: number }).priority === 1)).toBe(true);
    }
  });

  it('the lab’s scores are generated, in [0, 1], one a position of P04637, 1 at the hotspots', () => {
    const v = labScores();
    expect(v).toHaveLength(393);
    expect(Math.min(...v)).toBeGreaterThanOrEqual(0);
    expect(Math.max(...v)).toBe(1);
    for (const label of [175, 248, 273]) expect(v[label - 1]).toBe(1);
    expect(labScores()).toEqual(v);
  });
});
