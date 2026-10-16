/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Box3D } from '@molstar/core/math/geometry';
import { Task } from '@molstar/core/task';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { type ApiContext, ApiArea, ApiResult } from '@molstar/prism-viewer/api/area';
import { operationCommand } from '@molstar/prism-viewer/commands/operations';
import type { ViewerExtension, ViewerSkill } from '@molstar/prism-viewer/spec';
import { downloadBlob } from '@molstar/prism-viewer/ui/download';

// Geometry export in the Prism viewer (parity.md §6: Export › Geometry; README decision 27, commands first):
// the `geo-export` area with `geo-export.export`, which writes what the canvas shows (its render objects, as Mol*'s
// geometry exporters read them) as GLB, STL, OBJ or USDZ and gives the file's bytes; saving it is the page's
// (`ui/download.ts`), in the command and the Export card's Geometry section. The geometry is the canvas's, so it needs
// one (`modules: ['canvas']`: headless it fails `unavailable`). The exporters are a chunk of their own, imported when a
// format is first written. Like model-export, `undoable: false` without `history: false`: no history step, and the
// agent's policy asks every time (agent.md §4).

declare module '@molstar/prism-viewer/api' {
  interface ExtensionApis {
    'geo-export': GeoExportApi;
  }
}

export type GeometryFormat = 'glb' | 'stl' | 'obj' | 'usdz';

export const GeometryExportParams = {
  format: PD.Select<GeometryFormat>(
    'glb',
    [
      ['glb', 'GLB (glTF binary)'],
      ['stl', 'STL'],
      ['obj', 'OBJ'],
      ['usdz', 'USDZ (AR)'],
    ],
    { description: 'GLB keeps colors for most 3D tools; STL is shape only (3D printing); OBJ; USDZ for AR' },
  ),
  name: PD.Text('', { description: 'The file name, without its extension (empty: from the loaded entries)' }),
};

/** What `geo-export.export` gives: the file's bytes and name. */
export interface ExportedGeometry {
  blob: Blob;
  filename: string;
}

const clean = (s: string) =>
  s
    .trim()
    .replace(/[\\/:*?"<>|\s]+/g, '_')
    .replace(/^_+|_+$/g, '');

function sizeOf(bytes: number): string {
  return bytes < 1024 * 1024
    ? `${Math.max(1, Math.round(bytes / 1024))} KB`
    : `${(bytes / (1024 * 1024)).toFixed(1)} MB`;
}

/** The loaded entries' ids, as v5's export named the file. */
function entriesName(ctx: ApiContext): string {
  const ids = new Set<string>();
  for (const m of ctx.plugin.managers.structure.hierarchy.current.models) {
    const id = m.cell.obj?.data.entryId;
    if (id) ids.add(id.toUpperCase());
  }
  return [...ids].join('-') || 'molstar-model';
}

/** The format's exporter, from its own chunk. */
async function exporterFor(format: GeometryFormat, name: string, box: Box3D, radius: number) {
  switch (format) {
    case 'glb':
      return new (await import('./glb-exporter.js')).GlbExporter(box);
    case 'stl':
      return new (await import('./stl-exporter.js')).StlExporter(box);
    case 'obj':
      return new (await import('./obj-exporter.js')).ObjExporter(name, box);
    case 'usdz':
      return new (await import('./usdz-exporter.js')).UsdzExporter(box, radius);
  }
}

export class GeoExportApi extends ApiArea {
  readonly id = 'geo-export';

  readonly export = this.operation({
    id: 'geo-export.export',
    title: 'Export geometry',
    params: GeometryExportParams,
    modules: ['canvas'],
    undoable: false,
    run: async ({ format, name }, ctx): Promise<ApiResult<ExportedGeometry>> => {
      const canvas3d = ctx.plugin.canvas3d;
      if (!canvas3d) {
        return ApiResult.fail('unavailable', 'There is no canvas to export geometry from.', 'Use a viewer with one.');
      }
      // Hidden objects are skipped by the exporters too.
      const objects = canvas3d.getRenderObjects().filter((o) => o.state.visible);
      if (!objects.length) {
        return ApiResult.fail('not-found', 'Nothing is shown to export.', 'Load a structure, or show what is hidden.');
      }
      const base = clean(name.replace(/\.(glb|stl|obj|usdz|zip)$/i, '')) || entriesName(ctx);
      try {
        const sphere = canvas3d.boundingSphereVisible;
        const exporter = await exporterFor(format, base, Box3D.fromSphere3D(Box3D(), sphere), sphere.radius);
        const blob = await ctx.plugin.runTask(
          Task.create('Export geometry', async (rt) => {
            for (let i = 0, il = objects.length; i < il; i++) {
              if (ctx.signal.aborted) return undefined;
              await rt.update({
                message: `Exporting object ${i + 1}/${il}`,
                isIndeterminate: false,
                current: i,
                max: il,
              });
              await exporter.add(objects[i], canvas3d.webgl, rt);
            }
            return exporter.getBlob(rt);
          }),
          { useOverlay: true },
        );
        if (!blob || ctx.signal.aborted) return ApiResult.fail('aborted', 'Stopped.', 'Nothing was exported.');
        const filename = `${base}.${exporter.fileExtension}`;
        return ApiResult.ok(
          `Exported the geometry as ${filename} (${sizeOf(blob.size)}).`,
          {},
          {
            value: { blob, filename },
          },
        );
      } catch (e) {
        return ApiResult.fail('failed', `Could not export the geometry: ${e instanceof Error ? e.message : String(e)}`);
      }
    },
  });
}

export const GeoExportSkill: ViewerSkill = {
  id: 'export-geometry',
  title: 'Export geometry',
  summary: 'the 3D scene as GLB, STL, OBJ or USDZ (3D printing, AR, other 3D tools)',
  guidance: [
    'Skill export-geometry: export_geometry writes what is shown (its surfaces, cartoons, …) as a 3D file for the person.',
    '- format glb (colors, most 3D tools), stl (shape only, 3D printing), obj, usdz (AR on Apple devices).',
    '- Only what is visible is written: show or hide parts first.',
  ].join('\n'),
  tools: [
    {
      id: 'export_geometry',
      operation: 'geo-export.export',
      description: 'Writes the shown 3D scene as a GLB, STL, OBJ or USDZ file.',
    },
  ],
};

export const SectionId = 'geo-export.section.export';

export const GeoExportExtension: ViewerExtension = {
  id: 'geo-export',
  title: 'Geometry export',
  api: { 'geo-export': (ctx) => new GeoExportApi(ctx) },
  commands: (api) => {
    const x = api['geo-export'];
    if (!x) return [];
    return [
      {
        ...operationCommand(api, x.export, {
          category: 'Export',
          icon: 'box',
          title: 'Export geometry…',
          alias: ['glb', 'gltf', 'stl', 'obj', 'usdz', '3d print', 'mesh', 'export geometry'],
          action: 'Export',
          view: { title: 'Export › Geometry', card: 'core.overlay.export', section: SectionId },
        }),
        // The operation gives the file's bytes; the command saves it through the page.
        async run(args) {
          const r = await api['geo-export']!.export(args as never);
          if (!r.ok) throw new Error(r.hint ? `${r.message} ${r.hint}` : r.message);
          if (r.value) downloadBlob(r.value.blob, r.value.filename);
          return { message: r.did, ops: `geo-export.export(${JSON.stringify(args)})` };
        },
      },
      {
        id: 'geo-export.open',
        title: 'Show the Geometry export',
        category: 'Export',
        icon: 'box',
        alias: ['geometry export', '3d export'],
        // Opens the card: a view, not a change.
        async run() {
          const r = await api.view.open({ card: 'core.overlay.export', section: SectionId });
          if (ApiResult.isFailure(r)) throw new Error(r.message);
          return { message: 'Opened the Export card.', ops: `view.open({ card: "export", section: "${SectionId}" })` };
        },
      },
    ];
  },
  skills: [GeoExportSkill],
  ui: async (ctx) => (await import('./prism-viewer-ui.js')).geoExportUI(ctx),
};
