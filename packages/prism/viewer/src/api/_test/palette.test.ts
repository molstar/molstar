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

// The data palette (README §3.7): Prism's own list by default; `canvas.palette` swaps the categorical themes' default
// (Prism's, Mol*'s classic, colorblind-safe) and recolors what used the old default, as one step that Undo takes back,
// registry included; the spec gives the lists (viewer-spec.md §2.8).

const Okabe = [...ColorLists['okabe-ito'].list];
const Prism = [...ColorLists['many-distinct-prism'].list];

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
    // Prism's own at start.
    expect(chainColored(v).every((p) => JSON.stringify(p.colors) === JSON.stringify(Prism))).toBe(true);

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
    expect(chainColored(v).every((p) => JSON.stringify(p.colors) === JSON.stringify(Prism))).toBe(true);
    expectOk(await v.api.representation.color({ target: 'polymer', theme: 'chain-id' }));
    expect(chainColored(v).every((p) => JSON.stringify(p.colors) === JSON.stringify(Prism))).toBe(true);
  });

  it('classic is the themes’ own (Mol*’s many-distinct): no palette of ours on what it colors', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    expectOk(await v.api.representation.color({ target: 'polymer', theme: 'chain-id' }));
    expectOk(await v.api.canvas.palette({ palette: 'classic' }));
    const classic = chainColored(v);
    expect(classic.length).toBeGreaterThan(0);
    expect(classic.every((p) => JSON.stringify(p.colors) !== JSON.stringify(Prism))).toBe(true);
    expect(classic.every((p) => JSON.stringify(p.colors) !== JSON.stringify(Okabe))).toBe(true);
  });

  it('takes its lists from the spec (hex strings or numbers), and its choice from the options', async () => {
    const own = ['#112233', 0x445566, '#778899'];
    const v = await headlessViewer({ files: [cbs()], spec: { palettes: { prism: own } } });
    expectOk(await v.api.representation.color({ target: 'polymer', theme: 'chain-id' }));
    expect(
      chainColored(v).every((p) => JSON.stringify(p.colors) === JSON.stringify([0x112233, 0x445566, 0x778899])),
    ).toBe(true);

    const blind = await headlessViewer({
      files: [cbs()],
      spec: { options: { canvas: { palette: 'colorblind' } }, palettes: { colorblind: [0x000000, 0xffffff] } },
    });
    expectOk(await blind.api.representation.color({ target: 'polymer', theme: 'chain-id' }));
    expect(chainColored(blind).every((p) => JSON.stringify(p.colors) === JSON.stringify([0, 0xffffff]))).toBe(true);
  });
});
