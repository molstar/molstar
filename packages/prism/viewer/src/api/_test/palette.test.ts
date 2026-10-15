/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ColorLists } from '@molstar/core/util/color/lists';
import { describe, expect, it } from 'vitest';
import { cbs, crambin, expectOk, headlessViewer } from '../../_test/harness.js';
import type { PrismViewer } from '../../viewer.js';

// The data palette (README §3.7): `canvas.palette` swaps the categorical themes' default for Okabe–Ito and recolors
// what used the old default, as one step that Undo takes back, registry included.

const Okabe = [...ColorLists['okabe-ito'].list];

/** The palettes of the representations colored by a categorical theme, by theme. */
function palettes(v: PrismViewer): { theme: string; colors?: unknown }[] {
  return v.plugin.managers.structure.hierarchy.current.structures.flatMap((s) =>
    s.components.flatMap((c) =>
      c.representations.map((r) => {
        const t = (r.cell.transform.params as { colorTheme: { name: string; params?: { palette?: any } } }).colorTheme;
        return { theme: t.name, colors: t.params?.palette?.params?.list?.colors };
      }),
    ),
  );
}

const chainColored = (v: PrismViewer) => palettes(v).filter((p) => p.theme === 'chain-id');

describe('the data palette', () => {
  it('recolors what used the default, says when it is so already, and Undo takes it back with the registry', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    expectOk(await v.api.representation.color({ target: 'polymer', theme: 'chain-id' }));
    expect(chainColored(v).length).toBeGreaterThan(0);
    expect(chainColored(v).every((p) => JSON.stringify(p.colors) !== JSON.stringify(Okabe))).toBe(true);

    const r = await v.api.canvas.palette({ palette: 'colorblind' });
    expectOk(r);
    expect(chainColored(v).every((p) => JSON.stringify(p.colors) === JSON.stringify(Okabe))).toBe(true);
    expect((await v.api.canvas.palette({ palette: 'colorblind' })) as { unchanged?: boolean }).toMatchObject({
      unchanged: true,
    });

    // A structure loaded now takes it by default.
    expectOk(await v.api.load.files({ files: [crambin()] }));
    expectOk(await v.api.representation.color({ target: 'polymer', theme: 'chain-id', structure: '1CRN' }));
    expect(chainColored(v).every((p) => JSON.stringify(p.colors) === JSON.stringify(Okabe))).toBe(true);

    // The steps: the palette and the second load (coloring by chain is the default preset's: already so). Undo both:
    // the registry's default is the theme's own again.
    expect(v.managers.history.steps.value.map((s) => s.label)).toContain(
      'Chains and entities now use colorblind-safe colors',
    );
    for (let i = 0; i < 2; i++) expectOk(await v.api.history.undo({}));
    expect(chainColored(v).every((p) => JSON.stringify(p.colors) !== JSON.stringify(Okabe))).toBe(true);
    expectOk(await v.api.representation.color({ target: 'polymer', theme: 'chain-id' }));
    expect(chainColored(v).every((p) => JSON.stringify(p.colors) !== JSON.stringify(Okabe))).toBe(true);
  });
});
