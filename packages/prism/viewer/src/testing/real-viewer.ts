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
import type { SceneNode } from '../model/scene.js';
import type { PrismViewerSpec } from '../spec.js';
import type { PrismViewer } from '../viewer.js';

// The real viewer in a story: `createPrismViewer` over the given
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
  /** A docking receptor (PDB) and four poses of a ligand docked into it (one SDF, four molecules). */
  dockingReceptor: new URL('../../../../../data/examples/docking/receptor_1.pdb', import.meta.url).href,
  dockingPoses: new URL('../../../../../data/examples/docking/ligands_1.sdf', import.meta.url).href,
  /** A small protein's trajectory: its model (PDB) and coordinates (DCD). */
  trajectoryModel: new URL('../../../../../data/examples/trajectory/protein.pdb', import.meta.url).href,
  trajectoryCoordinates: new URL('../../../../../data/examples/trajectory/protein.dcd', import.meta.url).href,
  /** The same protein as Amber writes it: its topology (PRMTOP) and coordinates (NetCDF). */
  trajectoryTopology: new URL('../../../../../data/examples/trajectory/protein.parm7', import.meta.url).href,
  trajectoryNetcdf: new URL('../../../../../data/examples/trajectory/protein.nc', import.meta.url).href,
  /** A peptide and waters split across a dodecahedral box, 8 frames (multi-frame GRO; its generator next to it). */
  splitPeptide: new URL('../../../../../data/examples/trajectory/periodic/split-peptide.gro', import.meta.url).href,
  /** Particles (`data/examples/particles`, sources and licenses in its README): 200 ribosomes (CryoET Data Portal). */
  cryoetRibosomes: new URL('../../../../../data/examples/particles/cryoet-ribosomes.ndjson', import.meta.url).href,
  /** 300 HIV capsid picks with their poses, one tomogram (RELION STAR). */
  relionSubtomo: new URL('../../../../../data/examples/particles/relion-subtomo.star', import.meta.url).href,
  /** 300 picks in two tomograms (Dynamo table). */
  dynamo: new URL('../../../../../data/examples/particles/dynamo.tbl', import.meta.url).href,
  /** 225 ribosomes (Artiatomi motivelist, EM). */
  artiatomi: new URL('../../../../../data/examples/particles/artiatomi.em', import.meta.url).href,
  /** Five spheres over five frames (Simularium). */
  simularium: new URL('../../../../../data/examples/particles/small.simularium', import.meta.url).href,
};

/** The first scene node the predicate takes, depth first (what a story picks). */
export function findNode(nodes: readonly SceneNode[], is: (n: SceneNode) => boolean): SceneNode | undefined {
  for (const n of nodes) {
    if (is(n)) return n;
    const below = findNode(n.children.value, is);
    if (below) return below;
  }
  return undefined;
}

/**
 * The Structure panel open and a node picked (when given), its ancestors expanded so the tree shows where it is.
 */
export function openAt(v: PrismViewer, is?: (n: SceneNode) => boolean): void {
  v.views.shell?.show('structure');
  if (!is) return;
  const path = (nodes: readonly SceneNode[]): SceneNode[] | undefined => {
    for (const n of nodes) {
      if (is(n)) return [n];
      const below = path(n.children.value);
      if (below) return [n, ...below];
    }
    return undefined;
  };
  const found = path(v.views.scene.roots.value);
  if (!found) return;
  for (const n of found.slice(0, -1)) if (!n.isExpanded.value) n.toggleExpanded();
  v.views.scene.pick(found[found.length - 1]);
}

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
