/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import {
  bind,
  Button,
  domView,
  Icon,
  observableState,
  onActive,
  onCleanup,
  type ShellViewContext,
  Show,
  ViewHost,
} from '@molstar/prism-ui';
import { distinctUntilChanged, map } from 'rxjs';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import { fovNormalizedCameraPosition } from '@molstar/core/util/camera';
import type { CameraValue } from '../../api/camera.js';
import type { ExportFrame, ViewerExtensionUI, ViewerUIContext } from '../../spec.js';
import { ImageExport } from '../../views/export/image-export.js';
import { downloadBlob } from '../../ui/download.js';
import { createImageExportModel, frameValues } from './model.js';

// The Image section of the Export card (parity.md §4's export row): the existing view over a model built on
// `export.preview`, `export.image` and `export.image-copy`; the model is made when the section is shown and released
// with it. It previews only while the section is active (folded: no subscription, no render).
// The Camera section (v5's screenshot panel's Camera, parity.md): the camera in numbers while it is shown
// (`export.camera`, bound: folded, it follows nothing) and copying it as MolViewSpec's camera. v5's State group of that
// panel is Export › Session & state.

const round2 = (v: number) => Math.round(v * 100) / 100;
const list = (v: ArrayLike<number>) => Array.from(v, round2).join(', ');

const Rows: readonly [string, (c: CameraValue) => string][] = [
  ['Position', (c) => list(c.position)],
  ['FoV norm. position', (c) => list(fovNormalized(c))],
  ['Target', (c) => list(c.target)],
  ['Direction', (c) => list(direction(c))],
  ['Up', (c) => list(c.up)],
  ['Distance', (c) => `${c.distance} Å`],
  ['Radius', (c) => `${round2(c.radius)} Å`],
  ['Field of view', (c) => `${c.fov}° (${c.mode})`],
];

/** The position normalized by the field of view (MolViewSpec's, v5's "FoV Norm. Pos."). */
function fovNormalized(c: CameraValue): Vec3 {
  return fovNormalizedCameraPosition(Vec3.ofArray(c.target), Vec3.ofArray(c.position), c.mode, (c.fov * Math.PI) / 180);
}

/** Where the camera looks: a unit vector from its position to its target. */
function direction(c: CameraValue): Vec3 {
  return Vec3.normalize(Vec3(), Vec3.sub(Vec3(), Vec3.ofArray(c.target), Vec3.ofArray(c.position)));
}

/** MolViewSpec's camera: the position normalized by the field of view (as v5's "Copy MVS JSON"). */
function mvsCamera(c: CameraValue): string {
  const json = (v: ArrayLike<number>) => `[${list(v)}]`;
  return `{\n  "position": ${json(fovNormalized(c))},\n  "target": ${json(c.target)},\n  "up": ${json(c.up)}\n}`;
}

function cameraSection(ctx: ViewerUIContext) {
  const camera = observableState<CameraValue | undefined>(undefined);
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
        keywords: ['screenshot', 'png', 'picture', 'copy', 'animation', 'video', 'movie'],
        render: domView<ShellViewContext>(({ ctx: view }) => {
          const model = createImageExportModel({
            api: ctx.api,
            sources: ctx.sources,
            notify: (text, kind) => ctx.managers.notifications.post(text, kind),
            save: downloadBlob,
          });
          onActive((active) => model.setActive(active));
          onCleanup(() => model.dispose());
          // Opened at one of them (`export-image.section.image/<its id>`: the Export animation command).
          if (view.focus) model.form.set('mode', view.focus);
          // What else renders the view (an animation, `export.picture`): the switch's choices, each with its controls,
          // which get the shared frame (the resolution, the axes) as their subject.
          const pictures = view.shell.sections('export.picture');
          const frame: ExportFrame = { get: () => frameValues(model.form.values.value) };
          return (
            <ImageExport
              model={model}
              kinds={pictures.pipe(map((list) => list.map((p) => ({ id: p.id, title: p.title, icon: p.icon }))))}
              // Followed: a lazy extension's placeholder is replaced by its view once it loads (same id).
              kind={(id) => (
                <Show
                  when={pictures.pipe(
                    map((list) => list.find((x) => x.id === id)),
                    distinctUntilChanged(),
                  )}
                >
                  {(p) => <ViewHost renderer={p.render} ctx={{ shell: view.shell, subject: frame }} />}
                </Show>
              )}
            />
          );
        }),
      },
      {
        id: 'export-image.section.camera',
        slot: 'export',
        title: 'Camera',
        icon: 'move-3d',
        order: 0.5,
        keywords: ['camera', 'position', 'target', 'direction', 'mvs', 'molviewspec'],
        render: domView(() => cameraSection(ctx)),
      },
    ],
  };
}
