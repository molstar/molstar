/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { PrismViewer } from '../../viewer.js';

// The poses of the multi-viewer stories: one receptor and a ligand docked into it four ways, `data/examples/docking`.
// The poses are one SDF file of four molecules; each viewer gets one of them as an SDF of its own.

/** The molecules of an SDF file, each as a file of its own (ending with `$$$$`). */
export function splitSdf(text: string): string[] {
  return text
    .split(/^\$\$\$\$[ \t]*\r?\n?/m)
    .map((block) => block.replace(/^(\r?\n)+/, ''))
    .filter((block) => block.trim().length > 0)
    .map((block) => `${block.replace(/\s+$/, '')}\n$$$$\n`);
}

/** The name of a molecule in an SDF: its first line (the title), else `Pose n`. */
export function poseName(sdf: string, index: number): string {
  const title = sdf.split(/\r?\n/, 1)[0]?.trim();
  return title ? `${title} (pose ${index + 1})` : `Pose ${index + 1}`;
}

export interface PoseFiles {
  readonly receptor: string;
  readonly poses: readonly string[];
}

/** Fetches the receptor's address and the poses of the SDF at `posesUrl` (local files served with the stories). */
export async function fetchPoses(receptorUrl: string, posesUrl: string): Promise<PoseFiles> {
  const response = await fetch(posesUrl);
  if (!response.ok) throw new Error(`${posesUrl}: ${response.status}`);
  return { receptor: receptorUrl, poses: splitSdf(await response.text()) };
}

/**
 * The receptor and one pose, in a viewer: the receptor as a cartoon, the ligand as balls and sticks, the camera on the
 * ligand (every viewer's pose is in the receptor's frame, so their cameras can be one).
 */
export async function showPose(viewer: PrismViewer, files: PoseFiles, index: number): Promise<void> {
  const receptor = await viewer.api.load.url!({ url: files.receptor, format: 'pdb' });
  if (!receptor.ok) throw new Error(receptor.message);
  const sdf = files.poses[index];
  if (sdf === undefined) throw new Error(`There is no pose ${index + 1}.`);
  const ligand = await viewer.api.load.data!({ data: sdf, format: 'sdf', label: poseName(sdf, index) });
  if (!ligand.ok) throw new Error(ligand.message);
}
