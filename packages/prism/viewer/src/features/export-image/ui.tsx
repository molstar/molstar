/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { bind, Button, domView, Icon, observableState, onActive, onCleanup } from '@molstar/prism-ui';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import { fovNormalizedCameraPosition } from '@molstar/core/util/camera';
import type { CameraValue } from '../../api/camera.js';
import type { ViewerExtensionUI, ViewerUIContext } from '../../spec.js';
import { ImageExport } from '../../views/export/image-export.js';
import { downloadBlob } from '../../ui/download.js';
import { createImageExportModel } from './model.js';

// The Image section of the Export card (parity.md §4's export row): the existing view over a model built on
// `export.preview`, `export.image` and `export.image-copy`; the model is made when the section is shown and released
// with it. It previews only while the section is active (folded: no subscription, no render).
// The Viewport section (v5's screenshot panel's camera info and state, parity.md): the camera in numbers while it is
// shown (`export.camera`, bound: folded, it follows nothing), copying it as MolViewSpec's camera, and saving or loading
// the viewport state (`export.viewport`, `load.viewport`).

const Rows: readonly [string, (c: CameraValue) => string][] = [
  ['Position', (c) => c.position.join(', ')],
  ['Target', (c) => c.target.join(', ')],
  ['Up', (c) => c.up.join(', ')],
  ['Distance', (c) => `${c.distance} Å`],
  ['Field of view', (c) => `${c.fov}° (${c.mode})`],
];

const round2 = (v: number) => Math.round(v * 100) / 100;

/** MolViewSpec's camera: the position normalized by the field of view (as v5's "Copy MVS JSON"). */
function mvsCamera(c: CameraValue): string {
  const target = Vec3.ofArray(c.target);
  const position = fovNormalizedCameraPosition(target, Vec3.ofArray(c.position), c.mode, (c.fov * Math.PI) / 180);
  const list = (v: ArrayLike<number>) => `[${Array.from(v, round2).join(', ')}]`;
  return `{\n  "position": ${list(position)},\n  "target": ${list(c.target)},\n  "up": ${list(c.up)}\n}`;
}

function viewportSection(ctx: ViewerUIContext) {
  const camera = observableState<CameraValue | undefined>(undefined);
  const busy = observableState(false);
  const say = (r: { ok: boolean; did?: string; message?: string; hint?: string }) =>
    ctx.managers.notifications.post(
      r.ok ? (r.did ?? '') : `${r.message}${r.hint ? ` ${r.hint}` : ''}`,
      r.ok ? 'success' : 'error',
    );
  const save = async () => {
    if (busy.value || !ctx.api.export) return;
    busy.set(true);
    try {
      const r = await ctx.api.export.viewport({});
      if (r.ok && r.value) downloadBlob(r.value.blob, r.value.filename);
      say(r);
    } finally {
      busy.set(false);
    }
  };
  const file = (
    <input
      type="file"
      accept=".json,.molj"
      hidden
      onChange={(e) => {
        const picked = e.currentTarget.files?.[0];
        e.currentTarget.value = '';
        if (picked && ctx.api.load.viewport) void ctx.api.load.viewport({ file: picked }).then(say);
      }}
    />
  ) as HTMLInputElement;
  const copy = () => {
    const c = camera.value;
    if (!c || typeof navigator === 'undefined' || !navigator.clipboard) return;
    void navigator.clipboard.writeText(mvsCamera(c)).then(
      () => ctx.managers.notifications.post('Copied the camera as MolViewSpec JSON.', 'success'),
      () => ctx.managers.notifications.post('Could not copy to the clipboard.', 'error'),
    );
  };
  if (ctx.api.export) bind(ctx.api.export.camera, (c) => camera.set(c));
  const text = 'font-size: var(--msp-font-size-sm); min-width: 0';
  return (
    <div style="display: grid; gap: var(--msp-space-3); min-width: 0">
      <div
        style={`display: grid; grid-template-columns: auto 1fr; gap: var(--msp-space-1) var(--msp-space-3); ${text}`}
      >
        {Rows.map(([label, value]) => (
          <>
            <span style="color: var(--msp-color-text-muted)">{label}</span>
            <span
              style="font-family: var(--msp-font-mono); white-space: nowrap; overflow: hidden"
              ref={(el) => bind(camera.asReadonly, (c) => (el.textContent = c ? value(c) : '—'))}
            />
          </>
        ))}
      </div>
      <div style="display: flex; flex-wrap: wrap; gap: var(--msp-space-2)">
        <Button size="sm" icon={<Icon name="copy" />} title="Copy the camera as MolViewSpec JSON" onClick={copy}>
          Copy MVS camera
        </Button>
        <Button
          size="sm"
          icon={<Icon name="download" />}
          title="Save the camera and the canvas settings as JSON"
          disabled={busy.asReadonly}
          onClick={() => void save()}
        >
          Save
        </Button>
        <Button
          size="sm"
          icon={<Icon name="folder-open" />}
          title="Load a saved viewport state"
          disabled={!ctx.api.load.viewport}
          onClick={() => file.click()}
        >
          Load
        </Button>
        {file}
      </div>
    </div>
  );
}

export function exportImageUI(ctx: ViewerUIContext): ViewerExtensionUI {
  return {
    sections: [
      {
        id: 'export-image.section.image',
        slot: 'export',
        title: 'Image',
        icon: 'camera',
        order: 0,
        keywords: ['screenshot', 'png', 'picture', 'copy'],
        render: domView(() => {
          const model = createImageExportModel({
            api: ctx.api,
            sources: ctx.sources,
            notify: (text, kind) => ctx.managers.notifications.post(text, kind),
            save: downloadBlob,
          });
          onActive((active) => model.setActive(active));
          onCleanup(() => model.dispose());
          return <ImageExport model={model} />;
        }),
      },
      {
        id: 'export-image.section.viewport',
        slot: 'export',
        title: 'Viewport',
        icon: 'move-3d',
        order: 0.5,
        keywords: ['camera', 'viewport state', 'position', 'target', 'mvs', 'molviewspec'],
        render: domView(() => viewportSection(ctx)),
      },
    ],
  };
}
