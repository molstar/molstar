/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { combineLatest, map, type Observable, of } from 'rxjs';
import { ActionFooter, onCleanup } from '@molstar/prism-ui';
import { ParamsForm } from '@molstar/prism-ui';
import type { CropBox, ImageExportModel } from '../../model/image-export.js';
import { classes } from './_generated/image-export.module.css.js';

/**
 * Image export: the preview at the output's shape with the crop on it (auto-crop toggled in its corner), the screenshot
 * settings as a form, and the footer: set defaults, the size, Copy and Download.
 */
export function ImageExport(props: { model: ImageExportModel }) {
  const { model } = props;
  // Whether the preview matches the settings (the first one after opening may not yet).
  const settled = model.settled ?? of(true);
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
          <Crop model={model} settled={settled} />
          <button
            type="button"
            class={classes.autoCrop}
            aria-pressed={combineLatest([model.autoCrop.asReadonly, settled]).pipe(map(([a, s]) => a && s))}
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

/** The smallest crop, as a fraction of the image's side. */
const MinSide = 0.05;
/** An arrow key moves (or with ⇧ resizes) the box by this fraction of the image. */
const KeyStep = 0.01;
const Edges = ['n', 's', 'e', 'w', 'ne', 'nw', 'se', 'sw'] as const;
type Edge = (typeof Edges)[number] | 'move';

const pct = (v: number) => `${(v * 100).toFixed(2)}%`;
const clamp = (v: number, lo: number, hi: number) => Math.min(hi, Math.max(lo, v));

/** The box after a drag of `dx`, `dy` (fractions of the image) on an edge, a corner or the box itself. */
function dragged(start: CropBox, edge: Edge, dx: number, dy: number): CropBox {
  if (edge === 'move') {
    return {
      x: clamp(start.x + dx, 0, 1 - start.width),
      y: clamp(start.y + dy, 0, 1 - start.height),
      width: start.width,
      height: start.height,
    };
  }
  let left = start.x;
  let top = start.y;
  let right = start.x + start.width;
  let bottom = start.y + start.height;
  if (edge.includes('w')) left = clamp(left + dx, 0, right - MinSide);
  if (edge.includes('e')) right = clamp(right + dx, left + MinSide, 1);
  if (edge.includes('n')) top = clamp(top + dy, 0, bottom - MinSide);
  if (edge.includes('s')) bottom = clamp(bottom + dy, top + MinSide, 1);
  return { x: left, y: top, width: right - left, height: bottom - top };
}

/**
 * The crop box on the preview (v5's screenshot crop): drag it to move it, its edges and corners to resize it; arrow
 * keys move it and ⇧ arrow keys resize it. A drag turns auto-crop off (`setCrop`). The preview's size is read once
 * when a drag starts and the box is set at most once per frame (AGENTS.md, snappy).
 */
function Crop(props: { model: ImageExportModel; settled: Observable<boolean> }) {
  const { model } = props;
  let drag: { edge: Edge; width: number; height: number; start: CropBox; x0: number; y0: number } | undefined;
  let x = 0;
  let y = 0;
  let frame = 0;
  const apply = () => {
    frame = 0;
    if (drag) model.setCrop(dragged(drag.start, drag.edge, (x - drag.x0) / drag.width, (y - drag.y0) / drag.height));
  };
  const end = (e: PointerEvent) => {
    if (!drag) return;
    if (frame) cancelAnimationFrame(frame);
    apply();
    drag = undefined;
    const el = e.currentTarget as HTMLElement;
    delete el.dataset.dragging;
    if (el.hasPointerCapture?.(e.pointerId)) el.releasePointerCapture(e.pointerId);
  };
  onCleanup(() => {
    if (frame) cancelAnimationFrame(frame);
    drag = undefined;
  });
  return (
    <div
      class={classes.crop}
      tabIndex={0}
      role="group"
      aria-label="Crop: drag to move, drag an edge or a corner to resize; arrow keys move it, ⇧ arrow keys resize it"
      // Until the preview matches the settings, an automatic crop is not drawn (a stale preview would show the wrong
      // one); a crop the person set is theirs and stays.
      hidden={combineLatest([props.settled, model.autoCrop.asReadonly]).pipe(map(([s, a]) => a && !s))}
      data-auto={model.autoCrop.asReadonly.pipe(map((a) => (a ? '' : undefined)))}
      style={{
        left: model.crop.pipe(map((c) => pct(c.x))),
        top: model.crop.pipe(map((c) => pct(c.y))),
        width: model.crop.pipe(map((c) => pct(c.width))),
        height: model.crop.pipe(map((c) => pct(c.height))),
      }}
      onPointerDown={(e: PointerEvent) => {
        if (e.button !== 0 || drag) return;
        const el = e.currentTarget as HTMLElement;
        const preview = el.parentElement?.getBoundingClientRect();
        if (!preview?.width || !preview.height) return;
        const edge = ((e.target as HTMLElement).dataset.edge as Edge | undefined) ?? 'move';
        drag = {
          edge,
          width: preview.width,
          height: preview.height,
          start: { ...model.crop.value },
          x0: e.clientX,
          y0: e.clientY,
        };
        x = e.clientX;
        y = e.clientY;
        el.dataset.dragging = '';
        el.setPointerCapture?.(e.pointerId);
        e.preventDefault();
      }}
      onPointerMove={(e: PointerEvent) => {
        if (!drag) return;
        x = e.clientX;
        y = e.clientY;
        if (!frame) frame = requestAnimationFrame(apply);
      }}
      onPointerUp={end}
      onPointerCancel={end}
      onKeyDown={(e: KeyboardEvent) => {
        const step = { ArrowLeft: [-1, 0], ArrowRight: [1, 0], ArrowUp: [0, -1], ArrowDown: [0, 1] }[e.key];
        if (!step) return;
        e.preventDefault();
        model.setCrop(dragged(model.crop.value, e.shiftKey ? 'se' : 'move', step[0] * KeyStep, step[1] * KeyStep));
      }}
    >
      {Edges.map((edge) => (
        <span class={classes.handle} data-edge={edge} aria-hidden="true" />
      ))}
    </div>
  );
}
