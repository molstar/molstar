/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readFileSync } from 'node:fs';
import { Sphere3D } from '@molstar/core/math/geometry';
import type { GraphicsRenderObject } from '@molstar/graphics/gl/render-object';
import { DefaultPrismViewerSpec } from '@molstar/prism-viewer/default-spec';
import { createPrismViewerHeadless } from '@molstar/prism-viewer/headless';
import type { PrismViewer } from '@molstar/prism-viewer/viewer';
import { afterEach, describe, expect, it } from 'vitest';
import { GeoExportExtension } from '../prism-viewer.js';

// Geometry export in the Prism viewer on a headless viewer (parity.md §6), which has no canvas: the operation says so.
// The export itself is checked over a stand-in canvas that hands out the render objects the headless plugin built for
// 1CBS (geometry needs no WebGL to read; only texture meshes do). No network.

const cbs = () =>
  new File([readFileSync(new URL('../../../../data/examples/1cbs_updated.cif', import.meta.url))], '1cbs.cif');

const viewers: PrismViewer[] = [];
afterEach(async () => {
  for (const v of viewers.splice(0)) await v.dispose();
});

async function viewer(files: File[], o: { canvas?: boolean; extension?: boolean } = {}) {
  const base = DefaultPrismViewerSpec();
  const v = await createPrismViewerHeadless({
    spec: { ...base, extensions: [...(base.extensions ?? []), ...(o.extension === false ? [] : [GeoExportExtension])] },
    ...(o.canvas ? { modules: ['canvas', 'files', 'network'] as const } : {}),
  });
  viewers.push(v);
  for (const f of files) {
    const r = await v.api.load.files({ files: [f] });
    if (!r.ok) throw new Error(r.message);
  }
  return v;
}

/** The render objects of the scene's representations, as a canvas would hold them. */
function renderObjects(v: PrismViewer): GraphicsRenderObject[] {
  return v.plugin.managers.structure.hierarchy.current.structures.flatMap((s) =>
    s.components.flatMap((c) => c.representations.flatMap((r) => [...(r.cell.obj?.data.repr.renderObjects ?? [])])),
  );
}

/** A stand-in canvas: the render objects, the visible bounding sphere; no WebGL. */
function standInCanvas(v: PrismViewer, objects: GraphicsRenderObject[]) {
  const sphere = v.plugin.managers.structure.hierarchy.current.structures[0].cell.obj!.data.boundary.sphere;
  Object.defineProperty(v.plugin, 'canvas3d', {
    configurable: true,
    value: {
      getRenderObjects: () => objects,
      boundingSphereVisible: Sphere3D.clone(sphere),
      webgl: undefined,
      dispose: () => {},
    },
  });
}

describe('the geo-export extension in the Prism viewer', () => {
  it('fails `unavailable` without a canvas, as headless', async () => {
    const v = await viewer([cbs()]);
    const r = await v.api['geo-export']!.export({ format: 'glb' });
    expect(r.ok).toBe(false);
    if (!r.ok) {
      expect(r.code).toBe('unavailable');
      expect(r.message).toMatch(/needs/);
    }
    // Said as the reason the command and the tool cannot run.
    expect(v.api.describe('geo-export.export')?.available).not.toBe(true);
  });

  it('writes the shown geometry as GLB, STL, OBJ and USDZ, named for the entry', async () => {
    const v = await viewer([cbs()], { canvas: true });
    const objects = renderObjects(v);
    expect(objects.length).toBeGreaterThan(0);
    standInCanvas(v, objects);
    const steps = v.managers.history.steps.value.length;

    const glb = await v.api['geo-export']!.export({ format: 'glb' });
    expect(glb.ok, glb.ok ? '' : glb.message).toBe(true);
    if (!glb.ok) return;
    expect(glb.value!.filename).toBe('1CBS.glb');
    expect(glb.did).toMatch(/^Exported the geometry as 1CBS\.glb \(/);
    const head = new Uint8Array(await glb.value!.blob.arrayBuffer()).slice(0, 4);
    expect(new TextDecoder().decode(head)).toBe('glTF');
    expect(v.managers.history.steps.value.length).toBe(steps);

    const stl = await v.api['geo-export']!.export({ format: 'stl', name: 'print me' });
    expect(stl.ok && stl.value!.filename).toBe('print_me.stl');
    expect(stl.ok && stl.value!.blob.size).toBeGreaterThan(84);

    const obj = await v.api['geo-export']!.export({ format: 'obj' });
    // OBJ comes with its materials: a ZIP.
    expect(obj.ok && obj.value!.filename).toBe('1CBS.zip');

    const usdz = await v.api['geo-export']!.export({ format: 'usdz' });
    expect(usdz.ok && usdz.value!.filename).toBe('1CBS.usdz');
  });

  it('says when nothing is shown', async () => {
    const v = await viewer([cbs()], { canvas: true });
    standInCanvas(
      v,
      renderObjects(v).map((o) => ({ ...o, state: { ...o.state, visible: false } }) as GraphicsRenderObject),
    );
    const r = await v.api['geo-export']!.export({});
    expect(!r.ok && r.code).toBe('not-found');
  });

  it('brings its commands and its skill, and asks the agent every time (no undo restores a file)', async () => {
    const v = await viewer([]);
    expect(v.managers.commands.get('geo-export.export')).toBeDefined();
    expect(v.managers.commands.get('geo-export.open')).toBeDefined();
    const op = v.api['geo-export']!.export.definition;
    expect(op.undoable).toBe(false);
    expect(op.history).toBeUndefined();
    expect(op.modules).toEqual(['canvas']);
    expect(GeoExportExtension.skills?.[0].tools.map((t) => t.id)).toEqual(['export_geometry']);
    const bare = await viewer([], { extension: false });
    expect(bare.api['geo-export']).toBeUndefined();
  });
});
