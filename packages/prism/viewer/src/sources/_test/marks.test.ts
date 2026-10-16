/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Representation } from '@molstar/graphics/repr/representation';
import { describe, expect, it, vi } from 'vitest';
import { cbs, expectOk, headlessViewer } from '../../_test/harness.js';
import type { PrismViewer } from '../../viewer.js';

// The marks' paint (marks.md) follows what changed: the representations whose cells changed are painted again, an
// unrelated cell or a canvas setting paints nothing.

const frame = () => new Promise<void>((r) => setTimeout(r, 20));

function representations(v: PrismViewer) {
  return v.plugin.managers.structure.hierarchy.current.structures[0].components.flatMap((c) =>
    c.representations.map((r) => ({ ref: r.cell.transform.ref, repr: r.cell.obj!.data.repr as Representation<any> })),
  );
}

describe('the marks on the representations', () => {
  it('paints again only the representation that changed, and nothing for an unrelated cell or a canvas setting', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    expectOk(
      await v.api.highlight.mark({
        owner: 'plot',
        dimRest: true,
        marks: [{ id: 'a', tone: 'accent', target: 'A:1-10' }],
      }),
    );
    await frame();
    const reprs = representations(v);
    expect(reprs.length).toBeGreaterThan(1);
    const spies = reprs.map((r) => vi.spyOn(r.repr, 'setState'));
    const calls = () => spies.map((s) => s.mock.calls.length);

    // A cell that is no representation, and a canvas setting (a camera focus emits it too).
    await v.plugin.builders.data.rawData({ data: 'nothing', label: 'Unrelated' });
    v.plugin.events.canvas3d.settingsUpdated.next(undefined);
    await frame();
    expect(calls().every((n) => n === 0)).toBe(true);

    // One representation updated: it alone gets its marks back.
    const [first] = reprs;
    await v.plugin
      .build()
      .to(first.ref)
      .update((old: any) => ({ ...old, type: { ...old.type, params: { ...old.type.params, alpha: 0.5 } } }))
      .commit();
    await frame();
    expect(calls()[0]).toBeGreaterThan(0);
    expect(
      calls()
        .slice(1)
        .every((n) => n === 0),
    ).toBe(true);
    expect(first.repr.state.overpaint.layers.length).toBeGreaterThan(0);
  });
});
