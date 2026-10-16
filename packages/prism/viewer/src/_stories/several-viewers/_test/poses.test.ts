/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as fs from 'node:fs';
import * as path from 'node:path';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { headlessViewer, targetsOf } from '../../../_test/harness.js';
import { fetchPoses, poseName, showPose, splitSdf } from '../poses.js';

// The docking poses of the multi-viewer stories, from the repository's example files (the stories serve the same ones).

const dir = path.resolve(import.meta.dirname, '../../../../../../../data/examples/docking');
const sdfText = fs.readFileSync(path.join(dir, 'ligands_1.sdf'), 'utf8');

afterEach(() => vi.unstubAllGlobals());

describe('splitSdf', () => {
  it('gives each molecule of the file as an SDF of its own', () => {
    const poses = splitSdf(sdfText);
    expect(poses).toHaveLength(4);
    for (const p of poses) {
      expect(p.endsWith('$$$$\n')).toBe(true);
      expect(p.match(/\$\$\$\$/g)).toHaveLength(1);
      expect(p).toMatch(/V2000/);
    }
  });

  it('reads one with no end mark, and nothing from nothing', () => {
    expect(splitSdf('x\n  \n\n  0  0  0  0  0  0  0  0  0  0999 V2000\nM  END')).toHaveLength(1);
    expect(splitSdf('')).toEqual([]);
    expect(splitSdf('\n\n')).toEqual([]);
  });

  it('names a pose by its title', () => {
    expect(poseName('forcefield_3904\n', 0)).toBe('forcefield_3904 (pose 1)');
    expect(poseName('\n', 2)).toBe('Pose 3');
  });
});

describe('showPose', () => {
  it('loads the receptor and one pose, as two structures', async () => {
    const files: Record<string, string> = {
      'receptor_1.pdb': fs.readFileSync(path.join(dir, 'receptor_1.pdb'), 'utf8'),
      'ligands_1.sdf': sdfText,
    };
    vi.stubGlobal('fetch', async (input: string | URL | Request) => {
      const url = String(typeof input === 'object' && 'url' in input ? input.url : input);
      const name = Object.keys(files).find((k) => url.endsWith(k));
      return name ? new Response(files[name]) : new Response('not found', { status: 404 });
    });
    const v = await headlessViewer();
    const poses = await fetchPoses('http://stories.test/receptor_1.pdb', 'http://stories.test/ligands_1.sdf');
    expect(poses.poses).toHaveLength(4);
    await showPose(v, poses, 1);
    const structures = v.plugin.managers.structure.hierarchy.current.structures;
    expect(structures).toHaveLength(2);
    // The camera's target, the ligand, is the pose (the receptor has no ligand).
    const ligand = await targetsOf(v).resolve('ligand');
    expect(ligand.ok && ligand.parts.reduce((n, p) => n + p.selected.atoms, 0)).toBe(49);
    await expect(showPose(v, poses, 9)).rejects.toThrow(/no pose 10/);
  });
});
