/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Paul Pillot <paul.pillot@tandemai.com>
 */

import { Asset } from '@molstar/core/util/assets';
import { createPrismViewer } from '@molstar/prism-viewer/page';
import type { PrismViewer } from '@molstar/prism-viewer/viewer';
import { BondOrdersViewerSpec } from './spec.js';
import './index.html';

// The bond orders example page: a viewer with the bond-orders extension, a header that loads an entry and switches
// the structure's bonds through the `structure.bonds` operation (as the inspector's Bonds field does).

function element<T extends HTMLElement>(selector: string): T {
  const e = document.querySelector<T>(selector);
  if (!e) throw new Error(`The page has no ${selector}.`);
  return e;
}

/** What each button asks for: a provider and its params. */
const Choices: Record<string, { provider: string; params?: Record<string, unknown> }> = {
  default: { provider: 'default' },
  auto: { provider: 'bond-order-perception', params: { mode: 'auto' } },
  forceCompute: { provider: 'bond-order-perception', params: { mode: 'forceCompute' } },
};

async function main() {
  const status = element<HTMLElement>('#status');
  const viewer: PrismViewer = await createPrismViewer({ target: element('#app'), spec: BondOrdersViewerSpec() });
  (window as unknown as { viewer: PrismViewer }).viewer = viewer;

  const say = (r: { ok: boolean; did?: string; message?: string; hint?: string }) => {
    status.textContent = r.ok ? (r.did ?? '') : `${r.message}${r.hint ? ` ${r.hint}` : ''}`;
  };

  const load = async (id: string) => {
    const pdbId = id.trim().toUpperCase();
    status.textContent = `Loading ${pdbId}…`;
    await viewer.api.session.clear();
    try {
      await viewer.run(`Load ${pdbId}`, async (plugin) => {
        // the PDB format: no bond orders (CONECT records at most), so there is something to perceive; the mmCIF and
        // BinaryCIF files carry the ligands' `chem_comp_bond`, whose orders `auto` leaves as they are
        const data = await plugin.builders.data.download(
          { url: Asset.Url(`https://files.rcsb.org/download/${pdbId}.pdb`) },
          { state: { isGhost: true } },
        );
        const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
        await plugin.builders.structure.hierarchy.applyPreset(trajectory, 'default', {
          representationPreset: 'ball-and-stick' as never,
        });
      });
      status.textContent = `${pdbId}: the default bonds. Perceive the bond orders with the buttons or the structure's Bonds field.`;
    } catch (e) {
      status.textContent = `Could not load: ${e instanceof Error ? e.message : String(e)}`;
    }
  };

  for (const button of document.querySelectorAll<HTMLButtonElement>('button[data-bonds]')) {
    button.addEventListener('click', async () => {
      const choice = Choices[button.dataset.bonds!];
      say(await viewer.api.structure.bonds({ provider: choice.provider, params: choice.params ?? {} }));
    });
  }

  for (const button of document.querySelectorAll<HTMLButtonElement>('button[data-entry]')) {
    button.addEventListener('click', () => {
      element<HTMLInputElement>('#source').value = button.dataset.entry!;
      void load(button.dataset.entry!);
    });
  }

  element<HTMLFormElement>('#load').addEventListener('submit', (e) => {
    e.preventDefault();
    void load(element<HTMLInputElement>('#source').value);
  });
  await load(element<HTMLInputElement>('#source').value);
}

main().catch((e) => {
  console.error(e);
  const p = document.createElement('p');
  p.style.cssText = 'padding: 24px';
  p.textContent = `The viewer could not start: ${e instanceof Error ? e.message : String(e)}`;
  document.body.appendChild(p);
});
