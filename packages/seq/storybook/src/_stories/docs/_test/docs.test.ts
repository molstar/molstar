/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readdirSync, readFileSync, statSync } from 'node:fs';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';
import { checkSeqSpec } from '@molstar/seqx-core/lang';
import { GlobinsViews, P53Views } from '@molstar/seqx-examples';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { openExample } from '../../examples/_test/draw-checks.js';
import { sanitize, TitlePrefix } from '../../story-id.js';
import * as smallSpec from '../small-spec.js';

// The Docs pages: every page they link exists (a docs page by its title, as Storybook makes its id; the stories'
// links are `_test/story-links.test.ts`'s), every live example a page embeds is a story of the live file, and the views
// they show check clean; the Language page's own spec opens and draws headless.

const docs = join(dirname(fileURLToPath(import.meta.url)), '..');
const seq = join(docs, '..', '..', '..', '..');

function files(dir: string, out: string[] = []): string[] {
  for (const name of readdirSync(dir)) {
    if (name === 'node_modules' || name === 'dist' || name === 'lib') continue;
    const p = join(dir, name);
    if (statSync(p).isDirectory()) files(p, out);
    else out.push(p);
  }
  return out;
}

const all = files(seq);
const mdx = all.filter((f) => f.endsWith('.mdx'));
const docsPages = readdirSync(docs).filter((f) => f.endsWith('.mdx'));

/** A docs page's id, as Storybook makes it from its `<Meta title>`. */
const pageId = (title: string) => `${sanitize(`${TitlePrefix}/${title}`)}--docs`;

describe('the Docs pages', () => {
  it('are eight, each titled under Docs/', () => {
    const titles = docsPages.map((f) => /<Meta title="([^"]+)"/.exec(readFileSync(join(docs, f), 'utf8'))?.[1]);
    expect(titles.length).toBe(8);
    for (const t of titles) expect(t).toMatch(/^Docs\//);
  });

  it('every ?path=/docs/seq-next-… in the packages names a page that exists', () => {
    const pages = new Set<string>();
    for (const f of mdx) {
      const title = /<Meta title="([^"]+)"/.exec(readFileSync(f, 'utf8'))?.[1];
      if (title) pages.add(pageId(title));
    }
    expect(pages.has('seq-next-docs-introduction--docs')).toBe(true);
    const dead: string[] = [];
    let found = 0;
    for (const f of all.filter((x) => /\.(mdx|ts|md)$/.test(x))) {
      for (const m of readFileSync(f, 'utf8').matchAll(/\?path=\/docs\/(seq-next-[a-z0-9-]+--docs)/g)) {
        found++;
        if (!pages.has(m[1])) dead.push(`${f.slice(seq.length)}: ${m[1]}`);
      }
    }
    expect(found).toBeGreaterThan(20);
    expect(dead).toEqual([]);
  });

  it('embed only stories the live file exports', () => {
    const live = readFileSync(join(docs, 'live.stories.ts'), 'utf8');
    const exported = new Set([...live.matchAll(/^export const (\w+)/gm)].map((m) => m[1]));
    const missing: string[] = [];
    for (const f of docsPages)
      for (const m of readFileSync(join(docs, f), 'utf8').matchAll(/<Canvas of=\{Live\.(\w+)\}/g))
        if (!exported.has(m[1])) missing.push(`${f}: ${m[1]}`);
    expect(missing).toEqual([]);
  });

  it('show views that check clean', () => {
    const ids = [
      ...['tracks', 'protein-on-gene', 'text-under-tracks', 'gene-by-region'].map(
        (id) => P53Views.find((v) => v.id === id)!,
      ),
      GlobinsViews.find((v) => v.id === 'lazy')!,
    ];
    for (const v of ids) {
      expect(v).toBeDefined();
      expect(checkSeqSpec(v.build()).filter((i) => i.severity === 'error')).toEqual([]);
    }
  });
});

describe('the Language page’s spec', () => {
  let viewer: SeqViewer | undefined;
  afterEach(() => {
    viewer?.retire();
    viewer = undefined;
    vi.useRealTimers();
  });

  it('checks clean, opens and draws its three tracks headless', async () => {
    expect(checkSeqSpec(smallSpec.build())).toEqual([]);
    viewer = await openExample(
      { dataset: 'docs', id: 'small-spec', title: 'A small spec', summary: '', build: smallSpec.build },
      { wait: 100 },
    );
    expect(viewer.modelIssues.value).toEqual([]);
    const r = viewer.renderSvg({ idPrefix: 'd-' });
    await vi.advanceTimersByTimeAsync(400);
    const svg = await r;
    expect([svg.settled, svg.issues]).toEqual([true, []]);
    expect(svg.svg).toContain('TAD1');
  });
});
