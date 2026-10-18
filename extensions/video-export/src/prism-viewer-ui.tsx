/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { ActionFooter, Button, domView, Icon, ParamsForm, type ShellViewContext, Show } from '@molstar/prism-ui';
import { ParamsFormModel } from '@molstar/prism-ui/forms/params-form';
import { observableState } from '@molstar/prism-ui/reactive/observable-state';
import type { ExportFrame, ViewerExtensionUI, ViewerUIContext } from '@molstar/prism-viewer/spec';
import { downloadBlob } from '@molstar/prism-viewer/ui/download';
import { VideoOperations, VideoSection } from './prism-viewer-stub.js';

// Animation in the Export card's Image section (parity.md §6, `export.picture`): v5's Export Animation panel as a form
// over `video.export`'s own params (the animation, its length, frame rate, format, quality), Render, and Cancel while it
// renders; the frame (the resolution, the axes) is the section's, shared with the image (its subject, `ExportFrame`).
// It acts only through the person's API; the file it gets back is saved through the page. The animation's own params
// (spin speed, axis) stay with the operation's `params`, for commands and the agent.

const FormFields = ['animation', 'duration', 'fps', 'format', 'quality', 'name'] as const;

/** The operation's params as the form shows them (the frame is the section's). */
function formParams(ctx: ViewerUIContext): PD.Params {
  const params = ctx.api.describe(VideoOperations.export.id)?.params ?? {};
  const out: PD.Params = {};
  for (const k of FormFields) if (params[k]) out[k] = params[k];
  return out;
}

export function videoExportUI(ctx: ViewerUIContext): ViewerExtensionUI {
  return {
    sections: [
      {
        id: VideoSection,
        slot: 'export.picture',
        title: 'Animation',
        icon: 'clapperboard',
        order: 2,
        keywords: ['video', 'movie', 'mp4', 'webm', 'spin', 'rock', 'record'],
        render: domView<ShellViewContext>(({ ctx: view }) => {
          // The frame: the Image section's (without it, the viewport).
          const frame = view.subject as ExportFrame | undefined;
          // Made when the section shows: the animations registered by then are the choices.
          const form = new ParamsFormModel(formParams(ctx));
          const busy = observableState(false);
          const render = async () => {
            const video = ctx.api.video;
            if (busy.value || !video) return;
            busy.set(true);
            try {
              const r = await video.export({ ...form.values.value, ...frame?.get() } as never);
              if (!r.ok || !r.value) return;
              downloadBlob(r.value.blob, r.value.filename);
              ctx.managers.notifications.post(r.did, 'success');
            } finally {
              busy.set(false);
            }
          };
          const available = ctx.api.describe(VideoOperations.export.id)?.supported;
          return (
            <div style="display: grid; gap: var(--msp-space-3); min-width: 0">
              <ParamsForm form={form} />
              <Show when={busy.asReadonly}>
                {() => (
                  <Button size="sm" icon={<Icon name="x" />} onClick={() => void ctx.api.video?.cancel()}>
                    Cancel
                  </Button>
                )}
              </Show>
              <ActionFooter
                icon="clapperboard"
                label="Render"
                busy={busy.asReadonly}
                info={available === true ? '' : (available ?? '')}
                onReset={() => form.reset()}
                onApply={() => void render()}
              />
            </div>
          );
        }),
      },
    ],
  };
}
