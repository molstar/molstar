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
import type { PrismViewer } from '@molstar/prism-viewer/viewer';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { AssemblySymmetryExtension } from '../prism-viewer.js';
import { assemblySymmetryUI, SymmetryRow } from '../prism-viewer-ui.js';

// The assembly symmetry's row in the scene tree (parity.md §6; viewer-spec.md §3, `rows`) over a headless viewer, in
// jsdom. RCSB PDB is never reached: `fetch` is stubbed; the row itself never fetches.

const fixture = (name: string) =>
  new File([readFileSync(resolve(import.meta.dirname, `../../../../../data/examples/${name}`))], name);

const Gramicidin = [
  {
    kind: 'Global Symmetry',
    oligomeric_state: 'Homo 2-mer',
    stoichiometry: ['A2'],
    symbol: 'C2',
    type: 'Cyclic',
    clusters: [{ avg_rmsd: 0.4, members: [{ asym_id: 'A', pdbx_struct_oper_list_ids: ['1'] }, { asym_id: 'B' }] }],
    rotation_axes: [{ order: 2, start: [0, 0, -16], end: [0, 0, 16] }],
  },
];

const Monomer = [
  {
    kind: 'Global Symmetry',
    oligomeric_state: 'Monomer',
    stoichiometry: ['A'],
    symbol: 'C1',
    type: 'Asymmetric',
    clusters: [{ members: [{ asym_id: 'A' }] }],
  },
];

const fetchMock = vi.fn(async (_url: string, init?: RequestInit) => {
  const entry = JSON.parse(String(init?.body)).variables.entry_id as string;
  const symmetry = entry === '1GRM' ? Gramicidin : Monomer;
  return new Response(JSON.stringify({ data: { assembly: { rcsb_struct_symmetry: symmetry } } }), { status: 200 });
});

beforeEach(() => {
  fetchMock.mockClear();
  vi.stubGlobal('fetch', fetchMock);
});

const viewers: PrismViewer[] = [];
afterEach(async () => {
  for (const v of viewers.splice(0)) await v.dispose();
  vi.unstubAllGlobals();
  document.body.innerHTML = '';
});

async function viewer(file: string) {
  const base = DefaultPrismViewerSpec();
  const v = await createPrismViewerHeadless({
    spec: { ...base, extensions: [...(base.extensions ?? []), AssemblySymmetryExtension] },
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
  const ui = assemblySymmetryUI({ api: v.api.as({ kind: 'person' }), sources: v.sources, managers: v.managers });
  const host = document.createElement('div');
  document.body.append(host);
  // Mounted on the structure's Other group.
  const node = v.views.scene.roots.value[0].children.value.find((c) => c.kind === 'other')!;
  ui.rows![0].render.mount(host, { shell: new ShellModel(), subject: node }, undefined as never);
  return { ui, host, node };
}

describe("the assembly symmetry's row", () => {
  it("is a row in the structure's Other group; it asks RCSB PDB nothing until its eye is clicked, then shows the axes", async () => {
    const v = await viewer('1grm_updated.cif');
    const { ui, host, node } = mountRow(v);
    expect(ui.rows?.map((r) => r.id)).toEqual([SymmetryRow]);
    expect(ui.rows![0].when(node)).toBe(true);
    await until(() => host.textContent?.includes('Assembly symmetry') ?? false);
    expect(host.textContent).toContain('From RCSB PDB');
    // Opened before any fetch: no list yet, and still no request.
    (host.querySelector('[data-tree-row]') as HTMLElement).click();
    await until(() => host.textContent?.includes('get the symmetries') ?? false);
    expect(fetchMock).not.toHaveBeenCalled();

    (host.querySelector('button[data-tip]') as HTMLButtonElement).click();
    await until(() => host.textContent?.includes('C2 axes') ?? false);
    expect(fetchMock).toHaveBeenCalledTimes(1);
    expect(host.querySelector('[aria-label=Symmetry]')?.textContent).toContain('C2 · global · A2');
    expect(v.api['assembly-symmetry']!.rowState(v.views.scene.roots.value[0].id).value).toMatchObject({
      shown: true,
      selected: 1,
    });

    (host.querySelector('button[data-tip=Hide]') as HTMLButtonElement).click();
    await until(() => host.textContent?.includes('1 symmetry') ?? false);
    expect(v.api['assembly-symmetry']!.rowState(v.views.scene.roots.value[0].id).value.shown).toBe(false);
  });

  it('colors by cluster from its controls', async () => {
    const v = await viewer('1grm_updated.cif');
    const { host } = mountRow(v);
    (host.querySelector('[data-tree-row]') as HTMLElement).click();
    await until(() => host.textContent?.includes('Color by cluster') ?? false);
    const steps = v.managers.history.steps.value.length;
    (
      [...host.querySelectorAll('button')].find((b) => b.textContent === 'Color by cluster') as HTMLButtonElement
    ).click();
    await until(() => v.managers.history.steps.value.length === steps + 1);
    expect(v.managers.history.steps.value[0]?.label).toMatch(/cluster/);
  });

  it('says when RCSB PDB lists only C1, though nothing in the scene changed', async () => {
    const v = await viewer('1cbs_updated.cif');
    const { host } = mountRow(v);
    (host.querySelector('[data-tree-row]') as HTMLElement).click();
    await until(() => !!host.querySelector('button[data-tip]'));
    (host.querySelector('button[data-tip]') as HTMLButtonElement).click();
    await until(() => host.textContent?.includes('no symmetry but C1') ?? false);
    expect(host.textContent).toContain('1 symmetry');
  });

  it('is not shown for a structure that is not a PDB assembly', async () => {
    const v = await viewer('test.gro');
    const { ui, host, node } = mountRow(v);
    expect(v.api['assembly-symmetry']!.rowState(v.views.scene.roots.value[0].id).value.applicable).toBe(false);
    // Its `when` says so too: the Structure panel's rows are of one height, so a row that renders nothing would leave
    // an empty slot in the tree.
    expect(ui.rows![0].when!(node)).toBe(false);
    expect(host.textContent).toBe('');
  });
});
