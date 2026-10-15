/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { onCleanup } from '@molstar/prism-ui/dom/scope';
import { DefaultPrismViewerSpec } from '../default-spec.js';
import { createPrismViewer } from '../page.js';
import type { PrismViewerSpec } from '../spec.js';
import type { PrismViewer } from '../viewer.js';

// The real viewer in a story (plan phase 7: the mocks retired for what is built): `createPrismViewer` over the given
// spec, the repository's example files loaded through the API, then the story's own steps (open a tab, select, …).
// Disposed with the story.

/** The repository's example files, served with the stories (`data/examples`). */
export const ExampleFiles = {
  /** 1CBS: protein, retinoic acid, waters. */
  cbs: new URL('../../../../../data/examples/1cbs_updated.cif', import.meta.url).href,
  /** 1CRN: crambin. */
  crambin: new URL('../../../../../data/examples/1crn.cif', import.meta.url).href,
  /** 1GRM: gramicidin. */
  grm: new URL('../../../../../data/examples/1grm_updated.cif', import.meta.url).href,
  /** A small protein's trajectory: its model (PDB) and coordinates (DCD). */
  trajectoryModel: new URL('../../../../../data/examples/trajectory/protein.pdb', import.meta.url).href,
  trajectoryCoordinates: new URL('../../../../../data/examples/trajectory/protein.dcd', import.meta.url).href,
};

/** The element the viewer's shell lives in (its shadow root): for a story's clicks. */
export function shellRoot(host: Element): ShadowRoot | undefined {
  return host.querySelector('*')?.shadowRoot ?? host.shadowRoot ?? undefined;
}

export interface RealViewerStoryProps {
  /** Example files to load (`ExampleFiles`), in order. */
  readonly load?: readonly string[];
  /** The spec: the viewer's default by default; a story changes what it shows. */
  readonly spec?: (base: PrismViewerSpec) => PrismViewerSpec;
  /** After the loads: open a tab, select, set a background, … */
  readonly setup?: (viewer: PrismViewer) => void | Promise<void>;
  /** The host element's style (its size: the viewer fills it). */
  readonly style?: string;
}

/** A host element with the real viewer in it. */
export function RealViewerStory(props: RealViewerStoryProps) {
  const host = document.createElement('div');
  host.style.cssText = props.style ?? 'height: 100vh';
  let viewer: PrismViewer | undefined;
  let gone = false;
  onCleanup(() => {
    gone = true;
    void viewer?.dispose();
  });
  // Mounted first (the viewer measures its host), then made.
  queueMicrotask(async () => {
    const base = DefaultPrismViewerSpec();
    const made = await createPrismViewer({ target: host, spec: props.spec ? props.spec(base) : base });
    if (gone) return void made.dispose();
    viewer = made;
    (window as unknown as { molstar?: unknown }).molstar = { viewer: made };
    for (const url of props.load ?? []) {
      const r = await made.api.load.url!({ url, format: 'mmcif' });
      if (!r.ok) console.error(`[story] ${url}: ${r.message}`);
    }
    await props.setup?.(made);
  });
  return host;
}
