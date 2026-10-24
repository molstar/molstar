/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readFileSync } from 'node:fs';
import { resolve } from 'node:path';
import { ShellModel } from '@molstar/prism-ui/shell/shell';
import { DefaultPrismViewerSpec } from '@molstar/prism-viewer/default-spec';
import { createPrismViewerHeadless } from '@molstar/prism-viewer/headless';
import type { SceneNode } from '@molstar/prism-viewer';
import type { PrismViewer } from '@molstar/prism-viewer/viewer';
import { afterEach, describe, expect, it } from 'vitest';
import { AnvilExtension } from '../prism-viewer.js';
import { MembraneRow, membraneUI } from '../prism-viewer-ui.js';

// The membrane's row in the scene tree (parity.md §6; viewer-spec.md §3, `rows`) over a headless viewer, in jsdom:
// its state from the `membrane` area, its eye and its look acting through the person's API. Computed locally.

const fixture = (name: string) =>
  new File([readFileSync(resolve(import.meta.dirname, `../../../../../data/examples/${name}`))], name);

const viewers: PrismViewer[] = [];
afterEach(async () => {
  for (const v of viewers.splice(0)) await v.dispose();
  document.body.innerHTML = '';
});

async function viewer(file: string) {
  const base = DefaultPrismViewerSpec();
  const v = await createPrismViewerHeadless({
    spec: { ...base, extensions: [...(base.extensions ?? []), AnvilExtension] },
  });
  viewers.push(v);
  const r = await v.api.load.files({ files: [fixture(file)] });
  if (!r.ok) throw new Error(r.message);
  return v;
}

async function until(check: () => boolean, ms = 10000) {
  const t0 = Date.now();
  while (!check()) {
    if (Date.now() - t0 > ms) throw new Error('timed out');
    await new Promise((r) => setTimeout(r, 5));
  }
}

function mountRow(v: PrismViewer) {
  const ui = membraneUI({ api: v.api.as({ kind: 'person' }), sources: v.sources, managers: v.managers });
  const row = ui.rows![0];
  // Mounted on the structure's Other group.
  const node = v.views.scene.roots.value[0].children.value.find((c) => c.kind === 'other')!;
  const host = document.createElement('div');
  document.body.append(host);
  row.render.mount(host, { shell: new ShellModel(), subject: node }, undefined as never);
  return { ui, row, node, host };
}

describe("the membrane's row", () => {
  it("is a row in each structure's Other group, from the extension’s UI chunk", async () => {
    const v = await viewer('1grm_updated.cif');
    const ui = await AnvilExtension.ui!({ api: v.api, sources: v.sources, managers: v.managers });
    expect(ui.rows?.map((r) => r.id)).toEqual([MembraneRow]);
    const structure = v.views.scene.roots.value[0];
    const node = structure.children.value.find((c) => c.kind === 'other')!;
    expect(ui.rows![0].when(node)).toBe(true);
    expect(ui.rows![0].when(structure)).toBe(false);
    expect(ui.rows![0].when({ ...node, kind: 'component' } as SceneNode)).toBe(false);
  });

  it('shows and hides the membrane from its eye, and says how thick it is', async () => {
    const v = await viewer('1grm_updated.cif');
    const { host } = mountRow(v);
    await until(() => host.textContent?.includes('Membrane') ?? false);
    expect(host.textContent).toContain('Off');
    const eye = () => host.querySelector('button[data-tip]') as HTMLButtonElement;
    expect(eye().dataset.tip).toBe('Show');
    eye().click();
    await until(() => /Å thick/.test(host.textContent ?? ''));
    expect(eye().dataset.tip).toBe('Hide');
    expect(v.api.membrane!.rowState(v.views.scene.roots.value[0].id).value.look).toMatchObject({
      planes: true,
      rims: true,
    });
    eye().click();
    await until(() => host.textContent?.includes('Off') ?? false);
  });

  it('sets how the membrane looks from its controls', async () => {
    const v = await viewer('1grm_updated.cif');
    const ref = v.views.scene.roots.value[0].id;
    expect((await v.api.membrane!.show({})).ok).toBe(true);
    const { host } = mountRow(v);
    (host.querySelector('[data-tree-row]') as HTMLElement).click();
    await until(() => !!host.querySelector('[role=group]'));
    const label = [...host.querySelectorAll('[id]')].find((e) => e.textContent === 'Rims')!;
    const rims = host.querySelector(`[role=switch][aria-labelledby="${label.id}"]`) as HTMLButtonElement;
    expect(rims.getAttribute('aria-checked')).toBe('true');
    expect(host.textContent).toContain('Fit…');
    rims.click();
    await until(() => v.api.membrane!.rowState(ref).value.look?.rims === false);
  });

  it('drops a pending opacity change when its controls go', async () => {
    const v = await viewer('1grm_updated.cif');
    const ref = v.views.scene.roots.value[0].id;
    expect((await v.api.membrane!.show({})).ok).toBe(true);
    const { host } = mountRow(v);
    const row = host.querySelector('[data-tree-row]') as HTMLElement;
    row.click();
    await until(() => !!host.querySelector('[role=group]'));
    const before = v.api.membrane!.rowState(ref).value.look!.opacity;
    const slider = host.querySelector('input[type=range]') as HTMLInputElement;
    slider.value = '0.9';
    slider.dispatchEvent(new Event('input', { bubbles: true }));
    // Folded within the 150 ms: the call does not run afterwards.
    row.click();
    await until(() => !host.querySelector('[role=group]'));
    await new Promise((r) => setTimeout(r, 300));
    expect(v.api.membrane!.rowState(ref).value.look!.opacity).toBe(before);
  });

  it('is not shown for a structure ANVIL does not apply to', async () => {
    const v = await viewer('test.gro');
    const { host, row, node } = mountRow(v);
    expect(v.api.membrane!.rowState(v.views.scene.roots.value[0].id).value.applicable).toBe(false);
    // No slot in the tree for it (no blank row in Other), and nothing drawn if mounted anyway.
    expect(row.when(node)).toBe(false);
    expect(host.textContent).toBe('');
  });
});
