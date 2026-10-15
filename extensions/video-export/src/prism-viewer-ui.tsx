/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { ActionFooter, Button, domView, Icon, ParamsForm, Show } from '@molstar/prism-ui';
import { ParamsFormModel } from '@molstar/prism-ui/forms/params-form';
import { observableState } from '@molstar/prism-ui/reactive/observable-state';
import type { ViewerExtensionUI, ViewerUIContext } from '@molstar/prism-viewer/spec';
import { downloadBlob } from '@molstar/prism-viewer/ui/download';
import { VideoOperations, VideoSection } from './prism-viewer-stub.js';

// The Animation section of the Export card (parity.md §6): v5's Export Animation panel as a form over `video.export`'s
// own params (the animation, its length, frame rate, size, format, quality), Render, and Cancel while it renders. It
// acts only through the person's API; the file it gets back is saved through the page. The animation's own params
// (spin speed, axis) stay with the operation's `params`, for commands and the agent.

const FormFields = ['animation', 'duration', 'fps', 'size', 'width', 'height', 'format', 'quality', 'name'] as const;

/** The operation's params as the form shows them: a custom size's sides only for a custom size. */
function formParams(ctx: ViewerUIContext): PD.Params {
  const params = ctx.api.describe(VideoOperations.export.id)?.params ?? {};
  const out: PD.Params = {};
  for (const k of FormFields) {
    const p = params[k];
    if (!p) continue;
    out[k] = k === 'width' || k === 'height' ? { ...p, hideIf: (v: { size?: string }) => v.size !== 'custom' } : p;
  }
  return out;
}

export function videoExportUI(ctx: ViewerUIContext): ViewerExtensionUI {
  return {
    sections: [
      {
        id: VideoSection,
        slot: 'export',
        title: 'Animation',
        icon: 'clapperboard',
        order: 2,
        keywords: ['video', 'movie', 'mp4', 'webm', 'spin', 'rock', 'record'],
        render: domView(() => {
          // Made when the section shows: the animations registered by then are the choices.
          const form = new ParamsFormModel(formParams(ctx));
          const busy = observableState(false);
          const render = async () => {
            const video = ctx.api.video;
            if (busy.value || !video) return;
            busy.set(true);
            try {
              const r = await video.export(form.values.value as never);
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
