/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ActionFooter, Choices } from '@molstar/prism-ui';
import { ParamsForm } from '@molstar/prism-ui';
import type { IconName } from '@molstar/prism-ui/components/icon/icon';
import type { Child } from '@molstar/prism-ui/dom/jsx-types';
import { Show } from '@molstar/prism-ui/dom/control-flow';
import { effect, frameState } from '@molstar/prism-ui/dom/scope';
import { drag } from '@molstar/prism-ui/dom/scoped';
import type { ReadonlyState } from '@molstar/prism-ui/reactive/readonly-state';
import { derive, State } from '@molstar/prism-ui/reactive/state';
import type { CropBox, ImageExportModel } from '../../model/image-export.js';
import { classes } from './_generated/image-export.module.css.js';

/** What a model without a settled flag, or a view without kinds, gets: states that never change. */
const Settled: ReadonlyState<boolean> = new State(true).asReadonly;
const NoKinds: ReadonlyState<readonly PictureKind[]> = new State<readonly PictureKind[]>([]).asReadonly;

/** Besides an image, what else renders the view (an animation): a choice of the switch, and its controls. */
export interface PictureKind {
  readonly id: string;
  readonly title: string;
  readonly icon?: IconName;
}

/**
 * Image export: the preview at the output's shape with the crop on it (auto-crop toggled in its corner), the screenshot
 * settings as a form, and the footer: set defaults, the size, Copy and Download. With `kinds` (an animation), a switch
 * at the top: Image, or one of them, whose controls (`kind`) take the footer's place; the view (the resolution, the
 * crop, the axes) stays shared, the format and the transparency are the image's (the form's `mode`).
 */
export function ImageExport(props: {
  model: ImageExportModel;
  kinds?: ReadonlyState<readonly PictureKind[]>;
  kind?: (id: string) => Child;
}) {
  const { model } = props;
  // Whether the preview matches the settings (the first one after opening may not yet).
  const settled = model.settled ?? Settled;
  const mode = derive(model.form.values, (v) => (v.mode as string | undefined) || 'image');
  const kinds = props.kinds ?? NoKinds;
  return (
    <div class={classes.export}>
      <Show when={derive(kinds, (k) => (k.length ? k : undefined))}>
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
            'aspect-ratio': derive(model.aspect, (a) => `${a}`),
            background: derive(model.preview, (p) => p ?? 'var(--msp-color-control)'),
          }}
        >
          <Crop model={model} settled={settled} />
          <button
            type="button"
            class={classes.autoCrop}
            aria-pressed={derive([model.autoCrop.asReadonly, settled], (a, s) => a && s)}
            title="Crop to the visible content"
            onClick={() => model.autoCrop.set(!model.autoCrop.value)}
          >
            <span>Auto-crop</span>
          </button>
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
function Crop(props: { model: ImageExportModel; settled: ReadonlyState<boolean> }) {
  const { model } = props;
  let x = 0;
  let y = 0;
  let held: { edge: Edge; width: number; height: number; start: CropBox; x0: number; y0: number } | undefined;
  // The pointer's moves, counted: a frame state, so the box is set once a frame (its value is only a tick).
  const moves = frameState(0);
  let n = 0;
  const apply = () => {
    if (held) model.setCrop(dragged(held.start, held.edge, (x - held.x0) / held.width, (y - held.y0) / held.height));
  };
  effect(moves, apply);
  const el = (
    <div
      class={classes.crop}
      tabIndex={0}
      role="group"
      aria-label="Crop: drag to move, drag an edge or a corner to resize; arrow keys move it, ⇧ arrow keys resize it"
      // Until the preview matches the settings, an automatic crop is not drawn (a stale preview would show the wrong
      // one); a crop the person set is theirs and stays.
      hidden={derive([props.settled, model.autoCrop.asReadonly], (s, a) => a && !s)}
      data-auto={derive(model.autoCrop.asReadonly, (a) => (a ? '' : undefined))}
      style={{
        left: derive(model.crop, (c) => pct(c.x)),
        top: derive(model.crop, (c) => pct(c.y)),
        width: derive(model.crop, (c) => pct(c.width)),
        height: derive(model.crop, (c) => pct(c.height)),
      }}
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
  ) as HTMLElement;
  // The drag's listeners attached once (`drag`, dom/scoped.ts: every way it ends ends it); the preview's size read at
  // the press.
  drag(el, {
    start: (e) => {
      if (e.button !== 0) return undefined;
      const preview = el.parentElement?.getBoundingClientRect();
      if (!preview?.width || !preview.height) return undefined;
      const edge = ((e.target as HTMLElement).dataset.edge as Edge | undefined) ?? 'move';
      held = {
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
      e.preventDefault();
      return held;
    },
    move: (e) => {
      x = e.clientX;
      y = e.clientY;
      moves.offer(++n);
    },
    end: (_s, _how, e) => {
      // Ended by the pointer: where it is, now (a step still waiting for its frame is then a no-op); by the scope: not.
      if (e) apply();
      held = undefined;
      delete el.dataset.dragging;
    },
  });
  return el;
}
