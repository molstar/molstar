/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readFileSync } from 'node:fs';
import { checkSeqSpec } from '@molstar/seqx-core/lang';
import { recordDrawLists } from '@molstar/seqx-core/render';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { openExample } from '../../../examples/_test/draw-checks.js';
import { BuilderCases, shown } from '../builder-compose.js';
import * as check from '../check.js';

// The Builder stories, headless on a fake clock, opened as the stories open them (a viewer of the case's view): the
// code a story shows is its module's source, whole; every view opens and draws; `check()`'s issues are the ones the
// view's model says; the data panels print the spec the view loads.

let viewer: SeqViewer | undefined;

afterEach(() => {
  viewer?.retire();
  viewer = undefined;
  vi.useRealTimers();
});

const sourceOf = (id: string) => readFileSync(new URL(`../${id}.ts`, import.meta.url), 'utf8');

describe('Builder stories, headless', () => {
  for (const c of BuilderCases)
    it(`${c.id}: the code shown is the module's source; the view opens and draws`, async () => {
      expect(c.source).toBe(shown(sourceOf(c.id)));
      expect(c.source).toContain('createSeqSpec');
      viewer = await openExample(c, { gene: '' });
      expect(viewer.frame!.layout).toBeDefined();
      const svg = viewer.renderSvg({ idPrefix: 'b-' });
      await vi.advanceTimersByTimeAsync(400);
      const r = await svg;
      // the view of `check` is the placeholder of a spec with an issue: it has nothing to draw
      if (c.id === 'check') expect(r.issues.length + viewer.modelIssues.value.length).toBeGreaterThan(0);
      else {
        expect([r.settled, r.issues]).toEqual([true, []]);
        expect(viewer.modelIssues.value).toEqual([]);
        expect(recordDrawLists(r.lists, ['data']).split('\n').length).toBeGreaterThan(3);
      }
    });

  it('the handles story keeps its type error under @ts-expect-error, and builds the spec the types accept', () => {
    const c = BuilderCases.find((x) => x.id === 'handles')!;
    expect(c.source).toContain('@ts-expect-error');
    expect(checkSeqSpec(c.build()).filter((i) => i.severity === 'error')).toEqual([]);
  });

  it('check() lists the address that names nothing, and the panel prints it', () => {
    const issues = check.issues();
    expect(issues.length).toBeGreaterThan(0);
    expect(issues.some((i) => i.severity === 'error' && i.path?.includes('data'))).toBe(true);
    const panel = BuilderCases.find((x) => x.id === 'check')!.panels()[0];
    expect(panel.text).toContain(issues[0].code);
  });

  it('build() as data: the JSON panel parses to the spec, the S-expression prints it', () => {
    const c = BuilderCases.find((x) => x.id === 'as-data')!;
    const [json, sexpr] = c.panels();
    expect(JSON.parse(json.text)).toEqual(JSON.parse(JSON.stringify(c.build())));
    expect(sexpr.text).toContain('hello');
    expect(sexpr.text).toContain('letters');
  });
});
