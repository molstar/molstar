/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { ColorTheme } from '@molstar/graphics/theme/color';
import { UniformColorThemeProvider } from '@molstar/graphics/theme/color/uniform';
import { SizeTheme } from '@molstar/graphics/theme/size';
import { PhysicalSizeThemeProvider } from '@molstar/graphics/theme/size/physical';
import { UniformSizeThemeProvider } from '@molstar/graphics/theme/size/uniform';
import type { Structure } from '@molstar/model/model/structure';
import { atomStructure } from '../visual/scenes/structure.js';

export { exampleStructure } from '../visual/scenes/structure.js';

/** Volume per heavy atom in a protein, Å³ */
const AtomVolume = 11;

/**
 * A compact, protein-like ball of `count` carbon atoms: a jittered cubic lattice (spacing of the protein heavy-atom
 * density) cut to the `count` points nearest to the center, in four chains split along z
 */
export async function globuleStructure(count: number): Promise<Structure> {
  const spacing = Math.cbrt(AtomVolume);
  const radius = Math.cbrt((3 * count * AtomVolume) / (4 * Math.PI)) + 2 * spacing;
  const n = Math.ceil((2 * radius) / spacing);
  const points: [number, number, number, number][] = [];
  let seed = 1;
  const jitter = () => {
    seed = (seed * 1664525 + 1013904223) >>> 0;
    return (seed / 4294967296 - 0.5) * 0.5;
  };
  for (let i = 0; i < n; ++i) {
    for (let j = 0; j < n; ++j) {
      for (let k = 0; k < n; ++k) {
        const x = (i + 0.5) * spacing - radius + jitter();
        const y = (j + 0.5) * spacing - radius + jitter();
        const z = (k + 0.5) * spacing - radius + jitter();
        points.push([x, y, z, x * x + y * y + z * z]);
      }
    }
  }
  points.sort((a, b) => a[3] - b[3]);
  const atoms = points.slice(0, count).sort((a, b) => a[2] - b[2]);
  const x = new Float32Array(count),
    y = new Float32Array(count),
    z = new Float32Array(count);
  const chainOf = new Int32Array(count);
  for (let i = 0; i < count; ++i) {
    [x[i], y[i], z[i]] = atoms[i];
    chainOf[i] = Math.floor((i * 4) / count);
  }
  return atomStructure('globule', x, y, z, ['A', 'B', 'C', 'D'], chainOf);
}

let registries: { colorThemeRegistry: ColorTheme.Registry; sizeThemeRegistry: SizeTheme.Registry } | undefined;
/** Uniform color, uniform and physical size: the themes the benchmarked representations use */
export function getThemeRegistries() {
  if (!registries) {
    const colorThemeRegistry = ColorTheme.createRegistry();
    colorThemeRegistry.add(UniformColorThemeProvider);
    const sizeThemeRegistry = SizeTheme.createRegistry();
    sizeThemeRegistry.add(PhysicalSizeThemeProvider);
    sizeThemeRegistry.add(UniformSizeThemeProvider);
    registries = { colorThemeRegistry, sizeThemeRegistry };
  }
  return registries;
}

export function getPhysicalSizeTheme(structure: Structure) {
  return PhysicalSizeThemeProvider.factory(
    { structure },
    PD.getDefaultValues(PhysicalSizeThemeProvider.getParams({ structure })),
  );
}
