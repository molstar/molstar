/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Task } from '@molstar/core/task';
import { utf8ByteCount, utf8Write } from '@molstar/core/util/utf8';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { zip } from '@molstar/core/util/zip/zip';
import { type Structure, StructureElement, to_mmCIF, Unit } from '@molstar/model/model/structure';
import type { StructureRef } from '@molstar/plugin/state/manager/structure/hierarchy-state';
import { type ApiContext, ApiArea, ApiResult } from '@molstar/prism-viewer/api/area';
import { EveryStructureParam } from '@molstar/prism-viewer/api/params';
import { isFailure, labelOf, someStructures } from '@molstar/prism-viewer/api/scope';
import { operationCommand } from '@molstar/prism-viewer/commands/operations';
import type { ViewerExtension, ViewerSkill } from '@molstar/prism-viewer/spec';
import { downloadBlob } from '@molstar/prism-viewer/ui/download';
import { ModelExport } from './export.js';
import { to_pdb } from './pdb.js';

// Model export in the Prism viewer (plan 5.2, parity.md §6: Export › Models; README decision 27, commands first): the
// `model-export` area with `model-export.export`, which writes the loaded structures (or one, or only what is
// selected) as mmCIF, BinaryCIF or PDB and gives the file's bytes; saving it is the page's (`ui/download.ts`), in the
// command and the Export card's Models section. Exporting changes nothing in the scene, but a file written for the
// person is not something to do again unasked: `undoable: false` without `history: false`, so it is no history step and
// the agent's policy asks every time (agent.md §4).

declare module '@molstar/prism-viewer/api' {
  interface ExtensionApis {
    'model-export': ModelExportApi;
  }
}

export type ModelFormat = 'cif' | 'bcif' | 'pdb';

export const ModelExportParams = {
  structure: EveryStructureParam(),
  format: PD.Select<ModelFormat>(
    'cif',
    [
      ['cif', 'mmCIF'],
      ['bcif', 'BinaryCIF'],
      ['pdb', 'PDB'],
    ],
    { description: 'mmCIF keeps everything; BinaryCIF is the same, smaller; PDB only the coordinates' },
  ),
  selection: PD.Boolean(false, { description: 'Only the selected atoms' }),
  name: PD.Text('', {
    description: 'The file name, without its extension (empty: from the structure; a ZIP when there are several)',
  }),
};

/** What `model-export.export` gives: the file (one structure's, or a ZIP of several) and the files in it. */
export interface ExportedModels {
  blob: Blob;
  filename: string;
  files: readonly string[];
}

const Extensions: Record<ModelFormat, string> = { cif: 'cif', bcif: 'bcif', pdb: 'pdb' };

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

function bytesOf(data: string | Uint8Array<ArrayBuffer>): Uint8Array<ArrayBuffer> {
  if (typeof data !== 'string') return data;
  const bytes = new Uint8Array(utf8ByteCount(data));
  utf8Write(bytes, 0, data);
  return bytes;
}

/** What one structure exports as: the scene's (moved by a superposition), or its selected atoms; else why not. */
function exportedStructure(ctx: ApiContext, s: StructureRef, selection: boolean): Structure | string {
  const label = labelOf(ctx, s);
  const shown = s.transform?.cell.obj?.data ?? s.cell.obj?.data;
  if (!shown) return `${label} is not built yet.`;
  if (shown.models.length > 1) return `${label} holds several models; one at a time can be exported.`;
  if (shown.units.some((u) => !Unit.isAtomic(u))) return `${label} is a coarse model; only atomic ones export.`;
  if (!selection) return shown;
  const own = s.cell.obj?.data;
  const loci = own && ctx.plugin.managers.structure.selection.getLoci(own);
  if (!loci || !StructureElement.Loci.is(loci) || StructureElement.Loci.isEmpty(loci)) {
    return `Nothing is selected in ${label}.`;
  }
  return StructureElement.Loci.toStructure(StructureElement.Loci.remap(loci, shown));
}

function encode(name: string, s: Structure, format: ModelFormat, whole: boolean, warnings: Set<string>) {
  if (format === 'pdb') {
    const pdb = to_pdb(s);
    for (const w of pdb.warnings) warnings.add(w);
    return pdb.text;
  }
  // A whole structure keeps every category of its file; a selection only what Mol* writes of the atoms it keeps.
  return to_mmCIF(name, s, format === 'bcif', { copyAllCategories: whole });
}

export class ModelExportApi extends ApiArea {
  readonly id = 'model-export';

  readonly export = this.operation({
    id: 'model-export.export',
    title: 'Export models',
    params: ModelExportParams,
    undoable: false,
    when: (ctx) => (ctx.sources.empty.value ? 'Nothing is loaded.' : true),
    run: async ({ structure, format, selection, name }, ctx): Promise<ApiResult<ExportedModels>> => {
      const list = someStructures(ctx, structure || undefined);
      if (isFailure(list)) return list;

      const skipped: string[] = [];
      const parts: { name: string; structure: Structure }[] = [];
      for (const s of list) {
        const x = exportedStructure(ctx, s, selection);
        if (typeof x === 'string') skipped.push(x);
        else
          parts.push({ name: clean(ModelExport.getStructureName(x) || labelOf(ctx, s)) || 'structure', structure: x });
      }
      if (!parts.length) {
        const reason = skipped.join(' ');
        return selection
          ? ApiResult.fail('not-found', reason, 'Select atoms first, or export the whole structure.')
          : ApiResult.fail('invalid', reason || 'There is nothing to export.');
      }

      const ext = Extensions[format];
      const taken = new Map<string, number>();
      const warnings = new Set<string>();
      try {
        const files = await ctx.plugin.runTask(
          Task.create('Export models', async (rt) => {
            const out: [string, string | Uint8Array<ArrayBuffer>][] = [];
            for (const p of parts) {
              if (ctx.signal.aborted) break;
              const n = taken.get(p.name) ?? 0;
              taken.set(p.name, n + 1);
              const file = n ? `${p.name}_${n + 1}.${ext}` : `${p.name}.${ext}`;
              await rt.update({ message: `Exporting ${file}…`, isIndeterminate: true });
              out.push([file, encode(p.name, p.structure, format, !selection, warnings)]);
            }
            return out;
          }),
          { useOverlay: true },
        );
        if (ctx.signal.aborted) return ApiResult.fail('aborted', 'Stopped.', 'Nothing was exported.');

        const base = clean(name.replace(/\.(cif|bcif|pdb|zip)$/i, ''));
        let blob: Blob;
        let filename: string;
        if (files.length === 1) {
          blob = new Blob([files[0][1]]);
          filename = base ? `${base}.${ext}` : files[0][0];
        } else {
          const data: Record<string, Uint8Array<ArrayBuffer>> = {};
          for (const [f, d] of files) data[f] = bytesOf(d);
          const buffer = await ctx.plugin.runTask(Task.create('Compress models', (rt) => zip(rt, data)));
          blob = new Blob([new Uint8Array(buffer as ArrayBuffer)], { type: 'application/zip' });
          filename = `${base || 'structures'}.zip`;
        }
        const what = selection ? 'the selection of ' : '';
        const of = files.map((f) => f[0]).join(', ');
        const did =
          files.length === 1
            ? `Exported ${what}${parts[0].name} as ${filename} (${sizeOf(blob.size)}).`
            : `Exported ${what}${files.length} structures as ${filename} (${of}; ${sizeOf(blob.size)}).`;
        const all = [...skipped, ...warnings];
        return ApiResult.ok<ExportedModels>(
          did,
          {},
          {
            value: { blob, filename, files: files.map((f) => f[0]) },
            ...(all.length ? { warnings: all } : {}),
          },
        );
      } catch (e) {
        const big = format === 'cif' && parts.some((p) => p.structure.elementCount > 2_000_000);
        return ApiResult.fail(
          'failed',
          `Could not export: ${e instanceof Error ? e.message : String(e)}`,
          big ? 'The structure may be too big for text mmCIF: try BinaryCIF.' : undefined,
        );
      }
    },
  });
}

export const ModelExportSkill: ViewerSkill = {
  id: 'export-models',
  title: 'Export models',
  summary: 'structures as mmCIF, BinaryCIF or PDB files',
  guidance: [
    'Skill export-models: export_models writes structures as a file for the person to save.',
    '- format cif (mmCIF, everything), bcif (the same, binary), pdb (coordinates only).',
    '- selection: true writes only the selected atoms; select them first.',
    '- structure empty: every loaded one (several go into a ZIP).',
  ].join('\n'),
  tools: [
    {
      id: 'export_models',
      operation: 'model-export.export',
      description: 'Writes structures (or the selected atoms) as an mmCIF, BinaryCIF or PDB file.',
    },
  ],
};

export const SectionId = 'model-export.section.export';

export const ModelExportExtension: ViewerExtension = {
  id: 'model-export',
  title: 'Model export',
  api: { 'model-export': (ctx) => new ModelExportApi(ctx) },
  commands: (api) => {
    const x = api['model-export'];
    if (!x) return [];
    return [
      {
        ...operationCommand(api, x.export, {
          category: 'Export',
          icon: 'file',
          title: 'Export models…',
          alias: ['export models', 'mmcif', 'bcif', 'pdb', 'save coordinates', 'export structure'],
          action: 'Export',
          view: { title: 'Export › Models', card: 'core.overlay.export', section: SectionId },
        }),
        // The operation gives the file's bytes; the command saves it through the page.
        async run(args) {
          const r = await api['model-export']!.export(args as never);
          if (!r.ok) throw new Error(r.hint ? `${r.message} ${r.hint}` : r.message);
          if (r.value) downloadBlob(r.value.blob, r.value.filename);
          return { message: r.did, ops: `model-export.export(${JSON.stringify(args)})` };
        },
      },
      {
        id: 'model-export.open',
        title: 'Show the Models export',
        category: 'Export',
        icon: 'file',
        alias: ['model export', 'coordinates export'],
        // Opens the card: a view, not a change.
        async run() {
          const r = await api.view.open({ card: 'core.overlay.export', section: SectionId });
          if (ApiResult.isFailure(r)) throw new Error(r.message);
          return { message: 'Opened the Export card.', ops: `view.open({ card: "export", section: "${SectionId}" })` };
        },
      },
    ];
  },
  skills: [ModelExportSkill],
  ui: async (ctx) => (await import('./prism-viewer-ui.js')).modelExportUI(ctx),
};
