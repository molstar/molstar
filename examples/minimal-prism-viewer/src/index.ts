/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Asset } from '@molstar/core/util/assets';
import { createPrismViewer } from '@molstar/prism-viewer/page';
import type { PrismViewer } from '@molstar/prism-viewer/viewer';
import { MinimalPrismViewerSpec } from './spec.js';

// The minimal Prism viewer page: the viewer from `MinimalPrismViewerSpec`, and a header of its own (outside the viewer)
// that loads an entry and switches the structure between the cartoon and ball-and-stick presets.

type Look = 'cartoon' | 'ball-and-stick';

const Presets: Record<Look, string> = { cartoon: 'polymer-cartoon', 'ball-and-stick': 'ball-and-stick' };

function element<T extends HTMLElement>(selector: string): T {
  const e = document.querySelector<T>(selector);
  if (!e) throw new Error(`The page has no ${selector}.`);
  return e;
}

/** A PDB ID (RCSB's binary CIF) or the URL of an mmCIF file. */
function sourceUrl(source: string): { url: string; isBinary: boolean } {
  const text = source.trim();
  if (/^[0-9][a-z0-9]{3}$/i.test(text))
    return { url: `https://models.rcsb.org/${text.toUpperCase()}.bcif`, isBinary: true };
  const path = text.split(/[?#]/)[0].toLowerCase();
  return { url: text, isBinary: path.endsWith('.bcif') || path.endsWith('.bcif.gz') };
}

async function main() {
  const status = element<HTMLElement>('#status');
  const look = () => element<HTMLInputElement>('input[name=look]:checked').value as Look;

  const viewer: PrismViewer = await createPrismViewer({ target: element('#app'), spec: MinimalPrismViewerSpec() });
  (window as unknown as { viewer: PrismViewer }).viewer = viewer;
  performance.mark('minimal-prism-viewer:created');

  const load = async (source: string) => {
    const { url, isBinary } = sourceUrl(source);
    status.textContent = `Loading ${url}…`;
    await viewer.api.session.clear();
    try {
      // Raw plugin work as one history step (README §3.2): download, parse, and the default preset with the look.
      await viewer.run(`Load ${source}`, async (plugin) => {
        const data = await plugin.builders.data.download(
          { url: Asset.Url(url), isBinary },
          { state: { isGhost: true } },
        );
        const trajectory = await plugin.builders.structure.parseTrajectory(data, 'mmcif');
        await plugin.builders.structure.hierarchy.applyPreset(trajectory, 'default', {
          representationPreset: Presets[look()] as never,
        });
      });
      performance.mark('minimal-prism-viewer:structure-shown');
      status.textContent = url;
    } catch (e) {
      status.textContent = `Could not load: ${e instanceof Error ? e.message : String(e)}`;
    }
  };

  const show = async (value: Look) => {
    const result = await viewer.api.structure.preset({ preset: Presets[value] });
    if (!result.ok) status.textContent = result.message;
  };

  element<HTMLFormElement>('#load').addEventListener('submit', (e) => {
    e.preventDefault();
    void load(element<HTMLInputElement>('#source').value);
  });
  for (const radio of document.querySelectorAll<HTMLInputElement>('input[name=look]')) {
    radio.addEventListener('change', () => void show(look()));
  }

  await load(element<HTMLInputElement>('#source').value);
}

main().catch((e) => {
  console.error(e);
  const p = document.createElement('p');
  p.style.cssText = 'padding: 24px';
  p.textContent = `The viewer could not start: ${e instanceof Error ? e.message : String(e)}`;
  document.body.appendChild(p);
});
