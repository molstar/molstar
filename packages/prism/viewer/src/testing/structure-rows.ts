/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { MockRowSpec } from './scene.js';

const Regions: [string, string][] = [
  ['selection-box', 'Around focus'],
  ['camera-target', 'Camera target'],
  ['box', 'Fixed box'],
  ['cell', 'Whole structure'],
];

/** The volume streaming extension: density from the PDBe/RCSB servers, for X-ray and EM entries. */
export const VolumeStreamingRow: MockRowSpec = {
  id: 'volume-streaming',
  label: 'Volume streaming',
  icon: 'volume',
  when: (info) => info.method === 'x-ray' || info.method === 'em',
  params: () => ({
    enabled: PD.Boolean(false, { label: 'Enabled', isEssential: true }),
    view: PD.Select('selection-box', Regions, { label: 'Region', isEssential: true }),
    radius: PD.Numeric(
      5,
      { min: 0, max: 20, step: 0.5 },
      { label: 'Radius (Å)', isEssential: true, hideIf: (p: any) => p.view !== 'selection-box' },
    ),
    detail: PD.Select(
      '3',
      [
        ['0', 'Lowest'],
        ['3', 'Medium'],
        ['6', 'Highest'],
      ],
      { label: 'Detail', isEssential: true },
    ),
  }),
  describe: (v, info) =>
    `${info.method === 'em' ? 'EMDB' : 'X-ray'} · ${v.enabled ? (Regions.find(([r]) => r === v.view)?.[1] ?? '').toLowerCase() : 'off'}`,
  children: (_, info, make, id) =>
    info.method === 'em'
      ? [make.volume(`${id}-em`, 'EM', { relative: 3 }, 0x808080)]
      : [
          make.volume(`${id}-2fofc`, '2Fo-Fc', { relative: 1.5 }, 0x3362b2),
          make.volume(`${id}-fofc-pos`, 'Fo-Fc(+ve)', { relative: 3 }, 0x33bb33),
          make.volume(`${id}-fofc-neg`, 'Fo-Fc(−ve)', { relative: -3 }, 0xbb3333),
        ],
};

/** The RCSB assembly symmetry extension: the assembly's symmetry (annotated by RCSB PDB), its axes, coloring by it. */
export const AssemblySymmetryRow: MockRowSpec = {
  id: 'assembly-symmetry',
  label: 'Assembly symmetry',
  icon: 'orbit',
  when: () => true,
  params: (info) => {
    const sym = info.symmetry ?? 'C1';
    return {
      enabled: PD.Boolean(false, { label: 'Enabled', isEssential: true }),
      symmetry: PD.Select(
        'global',
        [
          ['global', `Global · ${sym}`],
          ['pseudo', `Pseudo · ${sym}`],
        ],
        { isEssential: true },
      ),
      axes: PD.Boolean(true, { label: 'Show axes', isEssential: true }),
      color: PD.Boolean(false, { label: 'Color by symmetry', isEssential: true }),
    };
  },
  describe: (v, info) => (v.enabled ? `${info.symmetry ?? 'C1'} · ${info.symmetry ? 'global' : 'no symmetry'}` : 'off'),
  children: (v, info, make, id) =>
    v.axes && info.symmetry ? [make.shape(`${id}-axes`, 'Symmetry axes', info.symmetry, 'orbit')] : [],
};

export const DefaultRows: readonly MockRowSpec[] = [VolumeStreamingRow, AssemblySymmetryRow];
