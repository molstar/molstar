/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { cbs, expectOk, headlessViewer } from '../../_test/harness.js';
import { ApiResult } from '../area.js';

// `annotations` (operations.md §2) and Home's loading actions (`load.action`).

describe('annotations', () => {
  it('lists what a structure carries and shows it as one step', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const listed = await v.api.annotations.list();
    expectOk(listed);
    expect(listed.did).toMatch(/B-factors \(b-factor\)/);
    expectOk(await v.api.annotations.apply({ annotation: 'b-factor' }));
    const themes = v.plugin.managers.structure.hierarchy.current.structures[0].components.flatMap((c) =>
      c.representations.map((r) => (r.cell.transform.params as any).colorTheme.name),
    );
    expect(themes).toContain('uncertainty');
    expect(v.managers.history.steps.value).toHaveLength(2);
    const missing = await v.api.annotations.apply({ annotation: 'plddt' });
    expect(missing).toMatchObject({ ok: false, code: 'not-found' });
  });

  it("takes extensions' providers", async () => {
    const v = await headlessViewer({
      files: [cbs()],
      extensions: [
        {
          id: 'demo',
          title: 'Demo',
          annotations: [
            {
              id: 'demo-score',
              title: 'Demo score',
              isApplicable: () => true,
              apply: async () => ApiResult.ok('Shown.', { nodes: [] }),
            },
          ],
        },
      ],
    });
    expectOk(await v.api.annotations.apply({ annotation: 'demo-score' }));
  });
});

describe("Home's actions", () => {
  it('lists the registry actions that load, and refuses an unknown one with the list', async () => {
    const v = await headlessViewer();
    const titles = v.views.homeActions.value.map((a) => a.id);
    expect(titles.some((id) => /download-structure/.test(id))).toBe(true);
    expect(titles.some((id) => /open-files/.test(id))).toBe(false);
    const r = await v.api.load.action({ action: 'nonsense' });
    expect(r).toMatchObject({ ok: false, code: 'not-found' });
    if (!r.ok) expect(r.hint).toMatch(/download-structure/);
  });
});
