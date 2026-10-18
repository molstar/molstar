/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Icon } from '@molstar/prism-ui/components/icon/icon';
import { bind, ParamsForm } from '@molstar/prism-ui';
import { For, Show } from '@molstar/prism-ui/dom/control-flow';
import { observableState } from '@molstar/prism-ui/reactive/observable-state';
import type { ReadonlyState } from '@molstar/prism-ui/reactive/readonly-state';
import { derive } from '@molstar/prism-ui/reactive/state';
import type { MeasurementKind, MeasurementPoint, MeasurementsModel } from '../../model/measurements.js';
import type { SceneNode } from '../../model/scene.js';
import { classes } from './_generated/measurements-add.module.css.js';
import { reorderRows } from '@molstar/prism-ui/components/reorder/drag-to-reorder';

/** How many of the latest selections a measurement can use. */
const MaxPoints = 4;

/**
 * Adding measurements (the Selection tab's Measure section, after Mol*'s measurement controls): the points (the latest
 * selections, in order, the first four used; reorder, remove, a click focuses one), what can be measured from them
 * (each enabled once there are enough points), and the options.
 */
export function MeasurementsAdd(props: { model: MeasurementsModel }) {
  const { model } = props;
  const used = derive(model.points, (p) => p.slice(0, MaxPoints));
  const older = derive(model.points, (p) => Math.max(0, p.length - MaxPoints));
  const options = observableState(false);
  let list: HTMLElement | undefined;
  const add = (
    <div class={classes.add}>
      <p class={classes.caption}>Points · the latest selections</p>
      <div class={classes.points} ref={(el: HTMLElement) => (list = el)}>
        <For
          each={used}
          key={(p) => p.id}
          fallback={() => (
            <p class={classes.hint}>
              Select atoms or residues (Select mode, Shift adds): the latest ones become the points to measure.
            </p>
          )}
        >
          {(point, index) => <PointRow model={model} point={point} index={index} />}
        </For>
        <Show when={derive(older, (n) => n > 0)}>
          {() => <p class={classes.more}>{derive(older, (n) => `+${n} older, not used`)}</p>}
        </Show>
      </div>
      <p class={classes.caption}>Measure</p>
      <div class={classes.kinds}>
        {model.kinds.map((kind) => (
          <KindRow model={model} kind={kind} />
        ))}
      </div>
      {/* The measurements made: adjust them here too (the same as their rows in the Structure tab). */}
      <Show when={derive(model.added, (a) => a.length > 0)}>
        {() => (
          <>
            <p class={classes.caption}>Added</p>
            <div class={classes.added}>
              <For each={model.added} key={(n) => n.id}>
                {(node) => <AddedRow model={model} node={node} />}
              </For>
            </div>
          </>
        )}
      </Show>
      <button
        type="button"
        class={classes.optionsToggle}
        aria-expanded={options}
        onClick={() => options.set(!options.value)}
      >
        <Icon name="chevron-right" size={14} class={classes.chevron} />
        Options
      </button>
      <Show when={options}>
        {() => {
          // The form follows the manager while it shows (a hidden section holds no subscription).
          if (model.optionValues) bind(model.optionValues, (v) => model.loadOptions?.(v));
          return <ParamsForm form={model.options} />;
        }}
      </Show>
    </div>
  );
  if (list) reorderRows(list, `.${classes.grip}`, (id, index) => model.movePointTo(id, index));
  return add;
}

function PointRow(props: { model: MeasurementsModel; point: MeasurementPoint; index: ReadonlyState<number> }) {
  const { model, point } = props;
  return (
    <div
      class={classes.point}
      data-id={point.id}
      onPointerEnter={() => model.highlightPoint?.(point.id)}
      onPointerLeave={() => model.highlightPoint?.(undefined)}
    >
      {/* The position; on hover, the handle that drags the point (the arrows stay for the keyboard). */}
      <span class={classes.slot}>
        <span class={classes.index}>{derive(props.index, (i) => `${i + 1}`)}</span>
        <span class={classes.grip} aria-hidden="true" data-tip="Drag to reorder">
          <Icon name="grip-vertical" size={13} />
        </span>
      </span>
      {/* The short form; the full label (where the pick is) in the tooltip. */}
      <button
        type="button"
        class={classes.pointLabel}
        title={point.label}
        aria-label={`Focus ${point.label}`}
        onClick={() => model.focusPoint(point.id)}
      >
        {point.short ?? point.label}
      </button>
      <span class={classes.pointActions}>
        <button
          type="button"
          class={classes.action}
          aria-label="Move up"
          title="Move up"
          onClick={() => model.movePoint(point.id, -1)}
        >
          <Icon name="arrow-up" size={13} />
        </button>
        <button
          type="button"
          class={classes.action}
          aria-label="Move down"
          title="Move down"
          onClick={() => model.movePoint(point.id, 1)}
        >
          <Icon name="arrow-down" size={13} />
        </button>
        <button
          type="button"
          class={classes.action}
          aria-label="Remove"
          title="Remove"
          onClick={() => model.removePoint(point.id)}
        >
          <Icon name="x" size={13} />
        </button>
      </span>
    </div>
  );
}

/** A measurement: its icon and name, and how many points it takes; enabled once there are enough. */
function KindRow(props: { model: MeasurementsModel; kind: MeasurementKind }) {
  const { model, kind } = props;
  const ready = derive(model.points, (p) => (kind.points === 0 ? p.length > 0 : p.length >= kind.points));
  const hint =
    kind.points === 0
      ? derive(ready, (r) => (r ? 'of the selection' : 'needs a selection'))
      : derive(ready, (r) =>
          r ? (kind.points === 1 ? 'point 1' : `points 1–${kind.points}`) : `needs ${kind.points} points`,
        );
  return (
    <button type="button" class={classes.kind} disabled={derive(ready, (r) => !r)} onClick={() => model.add(kind.id)}>
      <span class={classes.kindTitle}>{kind.title}</span>
      <span class={classes.kindHint}>{hint}</span>
    </button>
  );
}

/** A measurement made: its icon, label and value; show or hide, its style (inline under it), remove. */
function AddedRow(props: { model: MeasurementsModel; node: SceneNode }) {
  const { model, node } = props;
  const open = observableState(false);
  const style = model.styleOf(node);
  return (
    <div class={classes.addedItem}>
      <div class={classes.addedRow} data-hidden={derive(node.isHidden, (h) => (h ? '' : undefined))}>
        <span class={classes.pointLabel}>{node.label}</span>
        <span class={classes.addedValue}>{node.description}</span>
        <span class={classes.pointActions}>
          <button
            type="button"
            class={classes.addedAction}
            title={derive(node.isHidden, (h) => (h ? 'Show' : 'Hide'))}
            onClick={() => node.toggleVisibility()}
          >
            <Show when={node.isHidden} fallback={() => <Icon name="eye" size={14} />}>
              {() => <Icon name="eye-off" size={14} />}
            </Show>
          </button>
          {style ? (
            <button
              type="button"
              class={classes.addedAction}
              title="Style"
              aria-label="Style"
              aria-expanded={open}
              onClick={() => open.set(!open.value)}
            >
              <Icon name="sliders-horizontal" size={14} />
            </button>
          ) : null}
          <button
            type="button"
            class={classes.addedAction}
            title="Remove"
            aria-label="Remove"
            onClick={() => node.remove()}
          >
            <Icon name="x" size={13} />
          </button>
        </span>
      </div>
      {style ? (
        <Show when={open} reveal>
          {() => (
            <div class={classes.style}>
              <ParamsForm form={style} />
            </div>
          )}
        </Show>
      ) : null}
    </div>
  );
}
