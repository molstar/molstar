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
import { afterEach, describe, expect, it, vi } from 'vitest';
import { PaeRow, paeUI } from '../pae/ui.js';
import { ModelArchiveExtension } from '../prism-viewer.js';
import { archive, modelCif, Path, paeCif, Zip } from './associated-fixture.js';

// The PAE's row in the Structure panel (parity.md §6, "PAE plot"; viewer-spec.md §3, `rows`), in jsdom: under a
// structure whose PAE can be fetched (Model Archive's associated file) as under one that has it, as volume streaming's
// row is under a structure that can stream. Its click fetches the PAE when it is not here, then opens the card. Every
// request is a stub: no network.

const viewers: PrismViewer[] = [];
afterEach(async () => {
  vi.unstubAllGlobals();
  for (const v of viewers.splice(0)) await v.dispose();
  document.body.innerHTML = '';
});

async function viewer(file: File) {
  const base = DefaultPrismViewerSpec();
  const v = await createPrismViewerHeadless({
    spec: { ...base, extensions: [...(base.extensions ?? []), ModelArchiveExtension] },
  });
  viewers.push(v);
  const r = await v.api.load.files({ files: [file] });
  if (!r.ok) throw new Error(r.message);
  await until(() => v.views.scene.roots.value.length > 0);
  return v;
}

async function until(check: () => boolean, ms = 5000) {
  const t0 = Date.now();
  while (!check()) {
    if (Date.now() - t0 > ms) throw new Error('timed out');
    await new Promise((r) => setTimeout(r, 5));
  }
}

function mountRow(v: PrismViewer) {
  const person = v.api.as({ kind: 'person' });
  // The card's opening, seen: the row asks `view.open` for it.
  const open = vi.fn(async (_: { card: string }) => ({ ok: true }));
  const api = new Proxy(person, {
    get: (t, k) => (k === 'view' ? { ...t.view, open } : Reflect.get(t, k)),
  }) as typeof person;
  const ui = paeUI({ api, sources: v.sources, managers: v.managers });
  const row = ui.rows!.find((r) => r.id === PaeRow)!;
  // The structure's Other group, where the row is.
  const node = v.views.scene.roots.value[0].children.value.find((c) => c.kind === 'other')!;
  const host = document.createElement('div');
  document.body.append(host);
  row.render.mount(host, { shell: new ShellModel(), subject: node }, undefined as never);
  return { row, node, host, open };
}

describe("the PAE's row", () => {
  it('is under a structure whose PAE can be fetched; its click fetches it and opens the card', async () => {
    const v = await viewer(new File([modelCif({ kind: 'zip' })], 'ma-test-0001.cif'));
    const { row, node, host, open } = mountRow(v);
    expect(row.when(node)).toBe(true);
    expect(row.when({ ...node, kind: 'component' } as SceneNode)).toBe(false);
    expect(row.when(v.views.scene.roots.value[0])).toBe(false);
    await until(() => host.textContent?.includes('PAE') ?? false);
    expect(host.textContent).toContain('Not loaded · Model Archive');

    const body = await archive({ [Path]: paeCif() }, true);
    vi.stubGlobal(
      'fetch',
      vi.fn(async (u: unknown) =>
        String(u) === Zip ? new Response(body as BodyInit) : new Response('', { status: 404 }),
      ),
    );
    (host.querySelector('[data-tree-row]') as HTMLElement).click();
    await until(() => open.mock.calls.length > 0);
    expect(open).toHaveBeenCalledWith({ card: 'model-archive.overlay.pae' });
    expect(v.api['model-archive']!.hasPae(v.views.scene.roots.value[0].id)).toBe(true);
    await until(() => host.textContent?.includes('3 residues') ?? false);
    // Here now: a click opens the card without a request.
    vi.mocked(fetch).mockClear();
    (host.querySelector('[data-tree-row]') as HTMLElement).click();
    await until(() => open.mock.calls.length > 1);
    expect(fetch).not.toHaveBeenCalled();
  });

  it('is not under a structure with no PAE anywhere', async () => {
    const cbs = readFileSync(resolve(import.meta.dirname, '../../../../../data/examples/1cbs_updated.cif'));
    const v = await viewer(new File([cbs], '1cbs.cif'));
    const { row, node } = mountRow(v);
    expect(row.when(node)).toBe(false);
  });
});
