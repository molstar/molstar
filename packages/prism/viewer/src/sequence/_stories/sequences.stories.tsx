/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { story } from '@molstar/prism-ui/testing/index';
import { type Structure, Unit } from '@molstar/model/model/structure';
import type { SceneNode } from '../../model/scene.js';
import { ExampleFiles, findNode, openAt, RealViewerStory } from '../../testing/real-viewer.js';
import type { PrismViewer } from '../../viewer.js';
import type { ViewerExtension } from '../../spec.js';
import type { SequenceTrackProvider } from '../tracks.js';

// Each structure's sequence in the structure tree (sequence.md §8.1): a Sequence node under every structure with
// polymers, its eye the band (one structure at a time; Auto until the person picks one), its children the tracks
// (each an eye, its colorings, Remove), its menu Add track and Chains, its inspector the chain filter. The extension
// point for more sequence data is a track provider (§3): the last stories add one. The repository's examples (1CBS,
// 1CRN, 1TUP); no network.

const meta: Meta = { title: 'Viewer/Sequences', parameters: { layout: 'fullscreen' } };
export default meta;

const ok = (r: { ok: boolean; message?: string }) => {
  if (!r.ok) throw new Error(r.message);
};

/** Waits for a scene node (the Sequence rows come once the sequence's runtime is in). */
async function node(v: PrismViewer, is: (n: SceneNode) => boolean, ms = 10_000): Promise<SceneNode | undefined> {
  const t0 = performance.now();
  for (;;) {
    const n = findNode(v.views.scene.roots.value, is);
    if (n || performance.now() - t0 > ms) return n;
    await new Promise((r) => setTimeout(r, 50));
  }
}

/** The Sequence node of the structure at `index` (in load order). */
const sequenceAt = (v: PrismViewer, index: number) => (n: SceneNode) =>
  n.kind === 'sequence' &&
  n.id.startsWith(v.plugin.managers.structure.hierarchy.current.structures[index]?.cell.transform.ref ?? '#');

/** The Structure panel open at a structure's Sequence node, expanded (its tracks shown), picked (its inspector). */
async function openSequence(v: PrismViewer, index: number) {
  const found = await node(v, sequenceAt(v, index));
  if (!found) return;
  openAt(v, (n) => n === found);
  if (!found.isExpanded.value) found.toggleExpanded();
}

const load = async (v: PrismViewer, url: string, format: string) => ok(await v.api.load.url!({ url, format }));

function sequenceStory(o: {
  files: readonly [url: string, format: string][];
  steps?: (v: PrismViewer) => Promise<void>;
  open?: number;
  extensions?: readonly ViewerExtension[];
}): StoryObj {
  return {
    render: () =>
      story(() => (
        <RealViewerStory
          spec={(base) =>
            o.extensions ? { ...base, extensions: [...(base.extensions ?? []), ...o.extensions] } : base
          }
          setup={async (v) => {
            for (const [url, format] of o.files) await load(v, url, format);
            await v.api.camera.reset({ animate: false });
            await v.views.sequence?.load();
            await o.steps?.(v);
            await openSequence(v, o.open ?? 0);
          }}
        />
      )),
  };
}

/**
 * Two structures: each has its Sequence node, last under it. The band follows the first by the rule (Auto, said on
 * its row); the second's eye is off. Click the second's eye to show its sequence.
 */
export const TwoStructures: StoryObj = sequenceStory({
  files: [
    [ExampleFiles.cbs, 'mmcif'],
    [ExampleFiles.crambin, 'mmcif'],
  ],
});

/** The second structure picked by its eye (`sequence.follow`): the band shows crambin; Follow automatically in its menu */
export const PickedByItsEye: StoryObj = sequenceStory({
  files: [
    [ExampleFiles.cbs, 'mmcif'],
    [ExampleFiles.crambin, 'mmcif'],
  ],
  open: 1,
  steps: async (v) => {
    const second = v.plugin.managers.structure.hierarchy.current.structures[1].cell.transform.ref;
    ok(await v.api.sequence!.follow({ structure: second }));
  },
});

/** Tracks as children: the ruler removed (Add track brings it back), secondary structure hidden by its eye */
export const TracksEdited: StoryObj = sequenceStory({
  files: [[ExampleFiles.cbs, 'mmcif']],
  steps: async (v) => {
    ok(await v.api.sequence!.removeTrack!({ track: 'Ruler' }));
    ok(await v.api.sequence!.hideTrack!({ track: 'Secondary structure' }));
  },
});

/**
 * The chain filter in the inspector: 1TUP (p53 bound to DNA), the two DNA strands (A, B) and one of the three p53
 * copies (C; its author's A), named as `sequence.get` lists them
 */
export const ChosenChains: StoryObj = sequenceStory({
  files: [[ExampleFiles.tup, 'bcif']],
  steps: async (v) => ok(await v.api.sequence!.chains!({ chains: 'A, B, C' })),
});

/** Chains › Polymers only: 1TUP without its zinc ions */
export const PolymersOnly: StoryObj = sequenceStory({
  files: [[ExampleFiles.tup, 'bcif']],
  steps: async (v) => ok(await v.api.sequence!.chains!({ polymer: true })),
});

/** Kyte–Doolittle hydropathy by residue (J. Mol. Biol. 157, 105–132, 1982). */
const KyteDoolittle: Readonly<Record<string, number>> = {
  ILE: 4.5,
  VAL: 4.2,
  LEU: 3.8,
  PHE: 2.8,
  CYS: 2.5,
  MET: 1.9,
  ALA: 1.8,
  GLY: -0.4,
  THR: -0.7,
  SER: -0.8,
  TRP: -0.9,
  TYR: -1.3,
  PRO: -1.6,
  HIS: -3.2,
  GLU: -3.5,
  GLN: -3.5,
  ASP: -3.5,
  ASN: -3.5,
  LYS: -3.9,
  ARG: -4.5,
};

/**
 * A values table (Seq*'s `values` geometry) of the structure's residues' hydropathy: each polymer residue by its unit
 * (named as the sequence's complex names it: `label_asym_id`, `@<operator>` off the identity) and `label_seq_id`.
 */
function hydropathy(structure: Structure, complex: string) {
  const t = { target: [] as string[], pos: [] as string[], value: [] as number[] };
  const seen = new Set<string>();
  for (const unit of structure.units) {
    if (!Unit.isAtomic(unit)) continue;
    const h = unit.model.atomicHierarchy;
    const op = unit.conformation.operator;
    let last = -1;
    for (let i = 0; i < unit.elements.length; i++) {
      const e = unit.elements[i];
      const r = h.residueAtomSegments.index[e];
      if (r === last) continue;
      last = r;
      const kd = KyteDoolittle[h.atoms.label_comp_id.value(e)];
      if (kd === undefined) continue;
      const asym = h.chains.label_asym_id.value(h.chainAtomSegments.index[e]);
      const unitName = op.isIdentity ? asym : `${asym}@${op.instanceId}`;
      const seq = h.residues.label_seq_id.value(r);
      if (seen.has(`${unitName}:${seq}`)) continue;
      seen.add(`${unitName}:${seq}`);
      t.target.push(`${complex}/${unitName}`);
      t.pos.push(`${seq}`);
      t.value.push(kd);
    }
  }
  return t;
}

/**
 * An extension's track (sequence.md §3, a track provider): the residues' hydropathy as bars, data that comes with the
 * structure, computed in `build` and carried in the document as a values table (so the hover says each value). Out by
 * default: the Sequence node's menu offers it under Add track.
 */
const Hydropathy: SequenceTrackProvider = {
  id: 'story.sequence.hydropathy',
  label: 'Hydropathy',
  order: 15,
  track: 'story.hydropathy',
  shown: false,
  applies: (structure) => structure.polymerResidueCount > 0,
  describe: () => 'Bars · Kyte–Doolittle',
  build: (ctx) => ({
    data: [
      {
        kind: 'annotation',
        ref: 'story.hydropathy.data',
        params: {
          geometry: 'values',
          numbering: 'label',
          label: 'Hydropathy',
          value: hydropathy(ctx.structure, ctx.complex),
        },
      },
    ],
    track: {
      kind: 'track',
      ref: 'story.hydropathy',
      params: { label: '', height: 2 },
      children: [
        {
          kind: 'layer',
          ref: 'story.hydropathy.values',
          params: {
            type: 'values',
            data: 'story.hydropathy.data',
            encoding: 'bars',
            domain: [-4.5, 4.5],
            color: { theme: 'uniform', params: { color: '#5b8bd0' } },
          },
        },
      ],
    },
  }),
};
const StoryExtension: ViewerExtension = {
  id: 'story',
  title: 'Story tracks',
  sequence: { tracks: [Hydropathy] },
};

/** The extension's track offered: open the Sequence node's menu (…) › Add track › Hydropathy */
export const ExtensionTrackOffered: StoryObj = sequenceStory({
  files: [[ExampleFiles.cbs, 'mmcif']],
  extensions: [StoryExtension],
});

/** The extension's track added (`sequence.track.add`): bars under the letters, its value in the hover; a row of its own */
export const ExtensionTrackAdded: StoryObj = sequenceStory({
  files: [[ExampleFiles.cbs, 'mmcif']],
  extensions: [StoryExtension],
  steps: async (v) => ok(await v.api.sequence!.addTrack!({ track: 'Hydropathy' })),
});
