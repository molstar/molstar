/** @jsxImportSource @molstar/plugin-prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { map } from 'rxjs';
import { ActionFooter } from '../../components/action-footer/action-footer.js';
import { ParamsForm } from '../../components/params/params-form.js';
import type { ImageExportModel } from '../../model/image-export.js';
import { classes } from './_generated/image-export.module.css.js';

/**
 * Image export: the preview at the output's shape with the crop on it (auto-crop toggled in its corner), the screenshot
 * settings as a form, and the footer: set defaults, the size, Copy and Download.
 */
export function ImageExport(props: { model: ImageExportModel }) {
  const { model } = props;
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
            <span>Auto-crop</span>
          </button>
        </div>
      </div>
      <ParamsForm form={model.form} />
      <ActionFooter
        icon="download"
        label="Download"
        info={model.size}
        busy={model.busy}
        secondary={{ icon: 'copy', label: 'Copy', onClick: () => model.copy() }}
        onReset={() => model.form.reset()}
        onApply={() => model.download()}
      />
    </div>
  );
}
