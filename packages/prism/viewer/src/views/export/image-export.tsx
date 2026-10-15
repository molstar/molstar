/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { combineLatest, distinctUntilChanged, map, type Observable, of } from 'rxjs';
import { ActionFooter, Choices, onCleanup } from '@molstar/prism-ui';
import { ParamsForm } from '@molstar/prism-ui';
import type { IconName } from '@molstar/prism-ui/components/icon/icon';
import type { Child } from '@molstar/prism-ui/dom/jsx-types';
import { Show } from '@molstar/prism-ui/dom/control-flow';
import type { CropBox, ImageExportModel } from '../../model/image-export.js';
import { classes } from './_generated/image-export.module.css.js';

/** Besides an image, what else renders the view (an animation): a choice of the switch, and its controls. */
export interface PictureKind {
  readonly id: string;
  readonly title: string;
  readonly icon?: IconName;
}

/**
 * Image export: the preview at the output's shape with the crop on it (auto-crop toggled in its corner), the screenshot
 * settings as a form, and the footer: set defaults, the size, Copy and Download. With `kinds` (an animation), a switch
 * at the top: Image, or one of them, whose controls (`kind`) take the footer's place; the resolution and the axes stay
 * shared, the crop, the format and the transparency are the image's (the form's `mode`).
 */
export function ImageExport(props: {
  model: ImageExportModel;
  kinds?: Observable<readonly PictureKind[]>;
  kind?: (id: string) => Child;
}) {
  const { model } = props;
  // Whether the preview matches the settings (the first one after opening may not yet).
  const settled = model.settled ?? of(true);
  const mode = model.form.values.pipe(
    map((v) => (v.mode as string | undefined) || 'image'),
    distinctUntilChanged(),
  );
  const still = mode.pipe(map((m) => m === 'image'));
  const kinds = props.kinds ?? of([]);
  return (
    <div class={classes.export}>
      <Show when={kinds.pipe(map((k) => (k.length ? k : undefined)))}>
        {(k) => (
          <Choices
            choices={[
              { id: 'image', label: 'Image', icon: 'image' },
              ...k.map((x) => ({ id: x.id, label: x.title, ...(x.icon ? { icon: x.icon } : {}) })),
            ]}
            active={mode}
            onPick={(id) => model.form.set('mode', id)}
          />
        )}
      </Show>
      <div class={classes.stage}>
        <div
          class={classes.preview}
          style={{
            'aspect-ratio': model.aspect.pipe(map((a) => `${a}`)),
            background: model.preview.pipe(map((p) => p ?? 'var(--msp-color-control)')),
          }}
        >
          <Show when={still}>
            {() => (
              <>
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
              </>
            )}
          </Show>
        </div>
      </div>
      <ParamsForm form={model.form} />
      <Show when={mode}>
        {(m) =>
          m === 'image' || !props.kind ? (
            <ActionFooter
              icon="download"
              label="Download"
              info={model.size}
              busy={model.busy}
              secondary={{ icon: 'copy', label: 'Copy', onClick: () => model.copy() }}
              onReset={() => {
                model.form.reset();
                model.form.set('mode', m);
              }}
              onApply={() => model.download()}
            />
          ) : (
            props.kind(m)
          )
        }
      </Show>
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
