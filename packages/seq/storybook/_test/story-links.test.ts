/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readdirSync, readFileSync, statSync } from 'node:fs';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';
import { describe, expect, it } from 'vitest';
import { FeatureStories, FeatureTries, Figures, idOf } from '../src/_stories/examples/p53/where.js';
import {
  Figures as PlasmidFigures,
  idOf as plasmidId,
  Tries as PlasmidTries,
} from '../src/_stories/examples/plasmid/where.js';
import { sanitize, storyId, storyNameFromExport } from '../src/_stories/story-id.js';

// A moved story must not leave a dead link: every `?path=/story/seq-next-…` in the seq packages (the docs pages, the
// stories) names a story that exists, read from the stories' titles and exports as Storybook reads them.

const seq = join(dirname(fileURLToPath(import.meta.url)), '..', '..');

function files(dir: string, out: string[] = []): string[] {
  for (const name of readdirSync(dir)) {
    if (name === 'node_modules' || name === 'dist') continue;
    const p = join(dir, name);
    if (statSync(p).isDirectory()) files(p, out);
    else out.push(p);
  }
  return out;
}

const all = files(seq);

/** The ids of the stories: each `*.stories.ts` file's title and the names it exports (`export const`). */
function storyIds(): Set<string> {
  const ids = new Set<string>();
  for (const f of all.filter((x) => x.endsWith('.stories.ts'))) {
    const text = readFileSync(f, 'utf8');
    const title = /title:\s*'([^']+)'/.exec(text)?.[1];
    if (!title) continue;
    for (const m of text.matchAll(/^export const (\w+)\s*(?::[^=]+)?=/gm)) ids.add(storyId(title, m[1]));
  }
  return ids;
}

describe('story ids', () => {
  it('are made as Storybook makes them', () => {
    expect(storyId('Features/Layers', 'ValuesAsALineAndAnArea')).toBe(
      'seq-next-features-layers--values-as-a-line-and-an-area',
    );
    expect(storyId('Features/Navigators', 'GeneRegionOnChr17')).toBe(
      'seq-next-features-navigators--gene-region-on-chr-17',
    );
    expect(storyNameFromExport('Tracks')).toBe('Tracks');
    expect(sanitize('Seq* (next)/Examples/p53/Tracks')).toBe('seq-next-examples-p53-tracks');
  });
});

describe('the links to stories', () => {
  const ids = storyIds();

  it('find the stories (the ids are read from the files)', () => {
    expect(ids.has('seq-next-features-selection--tracks')).toBe(true);
    expect(ids.has('seq-next-examples-p53-tracks--tracks')).toBe(true);
    expect(ids.size).toBeGreaterThan(40);
  });

  it('every ?path=/story/seq-next-… in the packages names a story that exists', () => {
    const dead: string[] = [];
    let found = 0;
    for (const f of all.filter((x) => /\.(mdx|ts|md)$/.test(x) && !x.endsWith('story-links.test.ts'))) {
      const text = readFileSync(f, 'utf8');
      for (const m of text.matchAll(/\?path=\/story\/(seq-next-[a-z0-9-]+)/g)) {
        found++;
        if (!ids.has(m[1])) dead.push(`${f.slice(seq.length)}: ${m[1]}`);
      }
    }
    expect(found).toBeGreaterThan(10);
    expect(dead).toEqual([]);
  });

  it('the p53 index’s stories exist', () => {
    const missing = [...Object.values(Figures), ...Object.values(FeatureStories), ...FeatureTries.map((t) => t.where)]
      .map(idOf)
      .filter((id) => !ids.has(id));
    expect(missing).toEqual([]);
  });

  it('the Plasmid index’s stories exist', () => {
    const missing = [...Object.values(PlasmidFigures), ...PlasmidTries.map((t) => t.where)]
      .map(plasmidId)
      .filter((id) => !ids.has(id));
    expect(missing).toEqual([]);
  });
});
