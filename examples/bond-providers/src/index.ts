/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Asset } from '@molstar/core/util/assets';
import { createPrismViewer } from '@molstar/prism-viewer/page';
import type { PrismViewer } from '@molstar/prism-viewer/viewer';
import { NaiveOrdersBondProvider } from './providers.js';
import { BondProvidersViewerSpec } from './spec.js';
import './index.html';

// The bond providers example page: a viewer with the example's providers, a header that loads an entry and switches
// the structure's bonds through the `structure.bonds` operation (as the inspector's Bonds field does), and a button
// that unregisters and registers one provider at run time.

function element<T extends HTMLElement>(selector: string): T {
  const e = document.querySelector<T>(selector);
  if (!e) throw new Error(`The page has no ${selector}.`);
  return e;
}

/** What each button asks for: a provider and its params. */
const Choices: Record<string, { provider: string; params?: Record<string, unknown> }> = {
  default: { provider: 'default' },
  none: { provider: 'none' },
  'naive-distance': { provider: 'naive-distance', params: { cutoff: 1.9 } },
  'naive-orders': { provider: 'naive-orders' },
  // naive orders for the ligands, the built-in bonds for the rest and between units
  'per-unit': {
    provider: 'per-unit',
    params: {
      rules: [
        { match: { name: 'entity-type', params: ['non-polymer'] }, provider: { name: 'naive-orders', params: {} } },
      ],
      inter: { provider: { name: 'default', params: {} } },
    },
  },
};

async function main() {
  const status = element<HTMLElement>('#status');
  const viewer: PrismViewer = await createPrismViewer({ target: element('#app'), spec: BondProvidersViewerSpec() });
  (window as unknown as { viewer: PrismViewer }).viewer = viewer;

  const say = (r: { ok: boolean; did?: string; message?: string; hint?: string }) => {
    status.textContent = r.ok ? (r.did ?? '') : `${r.message}${r.hint ? ` ${r.hint}` : ''}`;
  };

  const load = async (id: string) => {
    const pdbId = id.trim().toUpperCase();
    status.textContent = `Loading ${pdbId}…`;
    await viewer.api.session.clear();
    try {
      // ball-and-stick everywhere, so the bonds are what you see
      await viewer.run(`Load ${pdbId}`, async (plugin) => {
        const data = await plugin.builders.data.download(
          { url: Asset.Url(`https://models.rcsb.org/${pdbId}.bcif`), isBinary: true },
          { state: { isGhost: true } },
        );
        const trajectory = await plugin.builders.structure.parseTrajectory(data, 'mmcif');
        await plugin.builders.structure.hierarchy.applyPreset(trajectory, 'default', {
          representationPreset: 'ball-and-stick' as never,
        });
      });
      status.textContent = `${pdbId}: the default bonds. Pick the structure in the Structure panel to see its Bonds field.`;
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

  // Run-time registration: the Bonds field and the operation follow the registry; a structure keeps its bonds until
  // its next update, which falls back to the default bonds.
  const toggle = element<HTMLButtonElement>('#toggle');
  toggle.addEventListener('click', () => {
    const registry = viewer.plugin.bondProviders;
    if (registry.has(NaiveOrdersBondProvider.name)) {
      // the spec's registration counts once: removing it once unregisters it
      registry.remove(NaiveOrdersBondProvider);
      toggle.textContent = 'Register naive orders';
      status.textContent = "'naive-orders' unregistered: the structures that use it keep their bonds until rebuilt.";
    } else {
      registry.add(NaiveOrdersBondProvider);
      toggle.textContent = 'Unregister naive orders';
      status.textContent = "'naive-orders' registered again.";
    }
  });

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
