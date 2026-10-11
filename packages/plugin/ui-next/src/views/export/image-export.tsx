/** @jsxImportSource @molstar/plugin-ui-next/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { map, of } from 'rxjs';
import { Button } from '../../components/button/button.js';
import { Icon } from '../../components/icon/icon.js';
import { ParamsForm } from '../../components/params/params-form.js';
import { Show } from '../../dom/control-flow.js';
import type { ImageExportModel } from '../../model/image-export.js';
import { classes } from './_generated/image-export.module.css.js';

/**
 * Image export: the preview at the output's shape with the crop on it (auto-crop toggled in its corner), the screenshot
 * settings as a form, and the size with Copy and Download at the bottom right.
 */
export function ImageExport(props: { model: ImageExportModel }) {
  const { model } = props;
  const busy = model.busy ?? of(false);
  const pct = (v: number) => `${(v * 100).toFixed(2)}%`;
  return (
    <div class={classes.export}>
      <div class={classes.stage}>
        <div
          class={classes.preview}
          style={{
            'aspect-ratio': model.aspect.pipe(map((a) => `${a}`)),
            background: model.preview.pipe(map((p) => p ?? 'var(--msp-color-control)')),
          }}
        >
          <div
            class={classes.crop}
            data-auto={model.autoCrop.asReadonly.pipe(map((a) => (a ? '' : undefined)))}
            style={{
              left: model.crop.pipe(map((c) => pct(c.x))),
              top: model.crop.pipe(map((c) => pct(c.y))),
              width: model.crop.pipe(map((c) => pct(c.width))),
              height: model.crop.pipe(map((c) => pct(c.height))),
            }}
          />
          <button
            type="button"
            class={classes.autoCrop}
            aria-pressed={model.autoCrop.asReadonly}
            title="Crop to the visible content"
            onClick={() => model.autoCrop.set(!model.autoCrop.value)}
          >
            <Icon name="crop" size={13} />
            <span>Auto-crop</span>
          </button>
        </div>
      </div>
      <ParamsForm form={model.form} />
      <div class={classes.footer}>
        <span class={classes.size}>{model.size}</span>
        <Button size="sm" variant="ghost" icon={<Icon name="copy" />} disabled={busy} onClick={() => model.copy()}>
          Copy
        </Button>
        <Button
          size="sm"
          variant="primary"
          disabled={busy}
          icon={
            <Show when={busy} fallback={() => <Icon name="download" />}>
              {() => <Icon name="loader-circle" class={classes.spin} />}
            </Show>
          }
          onClick={() => model.download()}
        >
          Download
        </Button>
      </div>
    </div>
  );
}
