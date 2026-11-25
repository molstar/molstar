/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as fs from 'node:fs';
import * as path from 'node:path';
import { fileURLToPath } from 'node:url';
import { describe, expect, it } from 'vitest';

const Source = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');

/**
 * The folders of graphics by layer, lowest first: a folder imports only from lower layers (and itself). `gl` draws,
 * themes map locations to colors and sizes, geometries take a theme to make their values, representations build
 * geometries, the canvas and the props draw and compute representations, `util` holds the top-level helpers (gizmos
 * and their visual sets)
 */
const Layers: { [folder: string]: number } = {
  gl: 0,
  theme: 1,
  geo: 2,
  formats: 3,
  repr: 3,
  canvas3d: 4,
  props: 4,
  util: 5,
};

const ImportPattern = /(?:import|export)\s+[^'";]*?from\s+'([^']+)'/g;

function getSourceFiles(dir: string, files: string[] = []) {
  for (const entry of fs.readdirSync(dir, { withFileTypes: true })) {
    const p = path.join(dir, entry.name);
    if (entry.isDirectory()) {
      if (entry.name !== '_test' && entry.name !== 'node_modules') getSourceFiles(p, files);
    } else if (/\.tsx?$/.test(entry.name)) {
      files.push(p);
    }
  }
  return files;
}

/** The top-level folder of graphics `specifier` (from `file`) points into, undefined outside graphics */
function getTargetFolder(file: string, specifier: string) {
  let target: string;
  if (specifier.startsWith('@molstar/graphics/')) target = specifier.slice('@molstar/graphics/'.length);
  else if (specifier.startsWith('.')) target = path.relative(Source, path.resolve(path.dirname(file), specifier));
  else return;
  return target.split(/[\\/]/)[0];
}

describe('graphics layers', () => {
  it('every folder imports only from lower layers', () => {
    const violations: string[] = [];
    for (const folder of Object.keys(Layers)) {
      for (const file of getSourceFiles(path.join(Source, folder))) {
        const text = fs.readFileSync(file, 'utf8');
        for (const [, specifier] of text.matchAll(ImportPattern)) {
          const target = getTargetFolder(file, specifier);
          if (!target || target === folder || !(target in Layers)) continue;
          if (Layers[target] >= Layers[folder]) {
            violations.push(`${path.relative(Source, file).split(path.sep).join('/')}: ${specifier}`);
          }
        }
      }
    }
    expect(violations).toEqual([]);
  });

  it('knows every top-level folder', () => {
    const folders = fs
      .readdirSync(Source, { withFileTypes: true })
      .filter((e) => e.isDirectory() && e.name !== '_test' && e.name !== 'node_modules')
      .map((e) => e.name);
    expect(folders.filter((f) => !(f in Layers))).toEqual([]);
  });
});
