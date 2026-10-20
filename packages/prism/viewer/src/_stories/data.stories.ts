/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { story } from '@molstar/prism-ui/testing/index';
import mvsFocus1cbs from '../../../../../data/examples/mvs/1cbs-focus.mvsj?raw';
import { RealViewerStory } from '../testing/real-viewer.js';
import type { PrismViewer } from '../viewer.js';

// One story per kind of data (README decision 19): each opens the data as the person would (RCSB's files, as a
// dropped file is) and calls the API for what shows it, then leaves the viewer usable. Most reach the network (RCSB,
// PDBe); tests never do.

const meta: Meta = { title: 'Viewer/Data', parameters: { layout: 'fullscreen' } };
export default meta;

/** The real viewer in a story, the setup run once it is made. */
const viewerStory = (setup: (viewer: PrismViewer) => Promise<void>) => story(() => RealViewerStory({ setup }));

/** The API as a person calls it (the person's own actions; the agent's would ask first for the network). */
const person = (viewer: PrismViewer) => viewer.api.as({ kind: 'person' });

/** Opens bytes as a dropped file is: the format is read from the name (`.bcif` and `.ply` are binary or text alike). */
async function loadFile(viewer: PrismViewer, data: BlobPart, name: string) {
  await person(viewer).load.files({ files: [new File([data], name)] as never });
}

/** Fetches a file and opens it as a dropped file is. */
async function loadUrl(viewer: PrismViewer, url: string, name: string) {
  const response = await fetch(url);
  if (!response.ok) throw new Error(`${url}: ${response.status}`);
  await loadFile(viewer, await response.blob(), name);
}

/** An entry's coordinates from the RCSB model server (binary CIF). */
const loadRcsb = (viewer: PrismViewer, id: string) =>
  loadUrl(viewer, `https://models.rcsb.org/${id.toUpperCase()}.bcif`, `${id.toLowerCase()}.bcif`);

/** A story's operation that must work: a failed result is thrown, so the story shows it rather than a bare canvas. */
async function must(result: Promise<{ ok: boolean; message?: string }>) {
  const r = await result;
  if (!r.ok) throw new Error(`Story step failed: ${r.message}`);
}

// --- Structures -------------------------------------------------------------------------------------------------

/** A small protein, crambin (1CRN, 46 residues), in the default cartoon. */
export const SmallProtein1CRN: StoryObj = {
  render: () => viewerStory((v) => loadRcsb(v, '1CRN')),
};

/** A protein with a ligand: 1CBS (retinoic acid in CRABP II), the ligand focused in the canvas. */
export const ProteinWithLigand1CBS: StoryObj = {
  render: () =>
    viewerStory(async (v) => {
      await loadRcsb(v, '1CBS');
      await must(person(v).focus.set({ target: 'ligand' }));
    }),
};

/** A nucleic acid: the Dickerson dodecamer (1BNA), the DNA in the polymer-and-ligand preset. */
export const NucleicAcid1BNA: StoryObj = {
  render: () =>
    viewerStory(async (v) => {
      await loadRcsb(v, '1BNA');
      await must(person(v).structure.preset({ preset: 'polymer-and-ligand' }));
    }),
};

/** A large complex: GroEL-GroES (1AON, 58,000 atoms, 21 chains), its assembly 1 made explicit. */
export const LargeComplex1AON: StoryObj = {
  render: () =>
    viewerStory(async (v) => {
      await loadRcsb(v, '1AON');
      await person(v).structure.assembly({ id: '1' });
    }),
};

// --- Many structures -------------------------------------------------------------------------------------------

/** Local entries of 1.2k–11.5k atoms (no network): protein–ligand, protein–DNA, multi-chain, an AlphaFold model. */
const ManySources = [
  new URL('../../../../../data/examples/1cbs_full.bcif', import.meta.url).href,
  new URL('../sequence/molstar/_test/fixtures/structures/1tup.bcif', import.meta.url).href,
  new URL('../sequence/molstar/_test/fixtures/structures/2ocj.bcif', import.meta.url).href,
  new URL('../sequence/molstar/_test/fixtures/structures/2xwr.bcif', import.meta.url).href,
  new URL('../sequence/molstar/_test/fixtures/structures/AF-P04637-F1-model_v6.bcif', import.meta.url).href,
  new URL('../../../../../data/examples/7qpd.fw2.cif', import.meta.url).href,
];

/**
 * `count` separate structures in one scene (perf-audit 14): the local entries above cycled, each its own file, opened
 * in one drop (`load.files`) and moved onto a grid (`structure.transform`). Nothing is shared between them, so the
 * scene has `count` structures and several representations each. The times show in the corner.
 */
function manyStructures(count: number) {
  return story(() => {
    const box = document.createElement('div');
    box.style.cssText = 'position: relative; height: 100vh';
    const status = document.createElement('div');
    status.style.cssText =
      'position: absolute; left: 8px; bottom: 8px; z-index: 10; padding: 4px 8px; border-radius: 4px; ' +
      'font: 12px/1.4 ui-monospace, monospace; background: rgba(0, 0, 0, 0.6); color: white; pointer-events: none';
    const show = (text: string) => {
      status.textContent = text;
      console.info(`[story] ${text}`);
    };
    show(`Fetching ${ManySources.length} files…`);
    const viewer = RealViewerStory({
      style: 'height: 100%',
      setup: async (v) => {
        const bytes = await Promise.all(
          ManySources.map(async (url) => {
            const response = await fetch(url);
            if (!response.ok) throw new Error(`${url}: ${response.status}`);
            return response.arrayBuffer();
          }),
        );
        const files = Array.from({ length: count }, (_, i) => {
          const source = ManySources[i % ManySources.length];
          const name = source.slice(source.lastIndexOf('/') + 1);
          return new File([bytes[i % bytes.length]], `${String(i + 1).padStart(3, '0')}-${name}`);
        });
        show(`Loading ${count} structures…`);
        const t0 = performance.now();
        await must(person(v).load.files({ files: files as never }));
        const loaded = performance.now() - t0;
        show(`${count} structures loaded in ${(loaded / 1000).toFixed(1)} s; arranging…`);

        // A grid, 80 Å apart: each structure moved from where it was loaded, its center to its cell
        const t1 = performance.now();
        const structures = v.api.structure.loaded.value;
        const side = Math.ceil(Math.sqrt(structures.length));
        await v.plugin.dataTransaction(async () => {
          for (let i = 0; i < structures.length; i++) {
            const data = v.plugin.state.data.cells.get(structures[i].ref)?.obj?.data as
              | { boundary: { sphere: { center: ArrayLike<number> } } }
              | undefined;
            if (!data) continue;
            const c = data.boundary.sphere.center;
            const translation = [80 * (i % side) - c[0], 80 * Math.floor(i / side) - c[1], -c[2]];
            await must(
              person(v).structure.transform({ structure: structures[i].ref, translation: translation as never }),
            );
          }
        });
        await must(person(v).camera.reset({ animate: false }));
        const arranged = performance.now() - t1;
        const atoms = structures.reduce((a, s) => a + s.atoms, 0);
        show(
          `${structures.length} structures, ${atoms.toLocaleString('en')} atoms: loaded in ${(loaded / 1000).toFixed(1)} s, ` +
            `arranged in ${(arranged / 1000).toFixed(1)} s`,
        );
      },
    });
    box.append(viewer, status);
    return box;
  });
}

/** 50 separate structures (1.2k–11.5k atoms each, about 200k in their default assemblies) in one scene, local files. */
export const ManyStructures50: StoryObj = {
  render: () => manyStructures(50),
};

// --- Models and trajectories ------------------------------------------------------------------------------------

/** An NMR ensemble (ubiquitin, 1D3Z, 10 models) opened on its first model, then `trajectory.frame` shows the fourth. */
export const NmrEnsemble1D3Z: StoryObj = {
  render: () =>
    viewerStory(async (v) => {
      await loadRcsb(v, '1D3Z');
      await must(person(v).trajectory.frame({ index: 3 }));
    }),
};

/** A multi-model file in the PDB format (2L22, 20 models from files.rcsb.org), the last model shown. */
export const MultiModelPdb2L22: StoryObj = {
  render: () =>
    viewerStory(async (v) => {
      await loadUrl(v, 'https://files.rcsb.org/download/2L22.pdb', '2l22.pdb');
      await must(person(v).trajectory.frame({ index: 19 }));
    }),
};

// --- Assemblies -------------------------------------------------------------------------------------------------

/** The asymmetric unit of 1A34 (satellite tobacco mosaic virus): one protein chain and its RNA, the deposited model. */
export const AsymmetricUnit1A34: StoryObj = {
  render: () =>
    viewerStory(async (v) => {
      await loadRcsb(v, '1A34');
      await person(v).structure.assembly({ id: '' });
    }),
};

/** The same entry as its biological assembly 1: sixty copies make the icosahedral capsid. */
export const Assembly1A34: StoryObj = {
  render: () =>
    viewerStory(async (v) => {
      await loadRcsb(v, '1A34');
      await must(person(v).structure.assembly({ id: '1' }));
    }),
};

// --- Maps -------------------------------------------------------------------------------------------------------

/** An EMDB map (EMD-8004) from PDBe's density server, by Home's `download-density` action: an isosurface, no model. */
export const EmdbMap: StoryObj = {
  render: () =>
    viewerStory((v) =>
      must(
        person(v).load.action({
          action: 'download-density',
          params: {
            source: { name: 'pdb-emd-ds', params: { provider: { id: 'emd-8004', server: 'pdbe' }, detail: 3 } },
          },
        }),
      ),
    ),
};

/** An X-ray structure (1TQN) with its cell's 2Fo-Fc and Fo-Fc maps from the density server (`download-density`). */
export const XrayDensityDownload1TQN: StoryObj = {
  render: () =>
    viewerStory(async (v) => {
      await loadRcsb(v, '1TQN');
      await must(
        person(v).load.action({
          action: 'download-density',
          params: { source: { name: 'pdb-xray-ds', params: { provider: { id: '1tqn', server: 'pdbe' }, detail: 3 } } },
        }),
      );
      await must(person(v).focus.set({ target: 'ligand' }));
    }),
};

// --- Particles and shapes ---------------------------------------------------------------------------------------

/**
 * Particles: a RELION STAR file made here (no particle data ships in the repository): 400 picks on a sphere, each with
 * its angles. The list opens at 1 Å radius, so the story makes them larger (`particles.representation`).
 */
export const ParticlesStar: StoryObj = {
  render: () =>
    viewerStory(async (v) => {
      await loadFile(v, relionStar(400, 200), 'particles.star');
      await must(person(v).particles.representation({ type: 'spacefill', params: { sizeFactor: 4 } }));
      // Data that is not a structure does not frame the camera on its own.
      await must(person(v).camera.reset({ animate: false }));
    }),
};

/** A shape: a PLY mesh made here (a torus with vertex colors), as a dropped `.ply` opens in the canvas. */
export const ShapePly: StoryObj = {
  render: () =>
    viewerStory(async (v) => {
      await loadFile(v, plyTorus(48, 24), 'torus.ply');
      await must(person(v).camera.reset({ animate: false }));
    }),
};

// --- MolViewSpec ------------------------------------------------------------------------------------------------

/** A MolViewSpec file (`data/examples/mvs/1cbs-focus.mvsj`: 1CBS with a labelled, zoomed ligand) opened as a file. */
export const MolViewSpec1CBS: StoryObj = {
  render: () => viewerStory((v) => loadFile(v, mvsFocus1cbs, '1cbs-focus.mvsj')),
};

// --- Generated fixtures -----------------------------------------------------------------------------------------

/** RELION STAR text: `count` particles on a sphere of `radius` (Å), the angles pointing out, in two classes. */
function relionStar(count: number, radius: number): string {
  const rows: string[] = [];
  const golden = Math.PI * (3 - Math.sqrt(5));
  for (let i = 0; i < count; i++) {
    const y = 1 - (2 * (i + 0.5)) / count;
    const r = Math.sqrt(1 - y * y);
    const phi = i * golden;
    const x = Math.cos(phi) * r;
    const z = Math.sin(phi) * r;
    const rot = (Math.atan2(z, x) * 180) / Math.PI;
    const tilt = (Math.acos(y) * 180) / Math.PI;
    rows.push(
      [x * radius, y * radius, z * radius, rot, tilt, 0, 1 + (i % 2)]
        .map((n) => (Number.isInteger(n) ? n : n.toFixed(3)))
        .join(' '),
    );
  }
  return [
    'data_particles',
    '',
    'loop_',
    '_rlnCenteredCoordinateXAngst #1',
    '_rlnCenteredCoordinateYAngst #2',
    '_rlnCenteredCoordinateZAngst #3',
    '_rlnAngleRot #4',
    '_rlnAngleTilt #5',
    '_rlnAnglePsi #6',
    '_rlnClassNumber #7',
    ...rows,
    '',
  ].join('\n');
}

/** ASCII PLY text of a torus (`u` steps around, `v` around the tube), colored by the angle around. */
function plyTorus(u: number, v: number, major = 20, minor = 7): string {
  const vertices: string[] = [];
  const faces: string[] = [];
  for (let i = 0; i < u; i++) {
    const a = (2 * Math.PI * i) / u;
    for (let j = 0; j < v; j++) {
      const b = (2 * Math.PI * j) / v;
      const x = (major + minor * Math.cos(b)) * Math.cos(a);
      const y = (major + minor * Math.cos(b)) * Math.sin(a);
      const z = minor * Math.sin(b);
      const red = Math.round(127 + 127 * Math.cos(a));
      const green = Math.round(127 + 127 * Math.cos(a + 2.094));
      const blue = Math.round(127 + 127 * Math.cos(a + 4.189));
      vertices.push(`${x.toFixed(3)} ${y.toFixed(3)} ${z.toFixed(3)} ${red} ${green} ${blue}`);
      const p = i * v + j;
      const right = ((i + 1) % u) * v + j;
      const up = i * v + ((j + 1) % v);
      const diagonal = ((i + 1) % u) * v + ((j + 1) % v);
      faces.push(`3 ${p} ${right} ${diagonal}`, `3 ${p} ${diagonal} ${up}`);
    }
  }
  return [
    'ply',
    'format ascii 1.0',
    'comment Torus',
    `element vertex ${vertices.length}`,
    'property float x',
    'property float y',
    'property float z',
    'property uchar red',
    'property uchar green',
    'property uchar blue',
    `element face ${faces.length}`,
    'property list uchar int vertex_index',
    'end_header',
    ...vertices,
    ...faces,
    '',
  ].join('\n');
}
