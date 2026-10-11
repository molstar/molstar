/** @jsxImportSource @molstar/plugin-ui-next/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { distinctUntilChanged, map } from 'rxjs';
import { Icon } from '../../components/icon/icon.js';
import { ParamsForm } from '../../components/params/params-form.js';
import { dragToReorder } from '../../components/reorder/drag-to-reorder.js';
import { For, Show } from '../../dom/control-flow.js';
import type { MeasurementKind, MeasurementPoint, MeasurementsModel } from '../../model/measurements.js';
import type { SceneNode } from '../../model/scene.js';
import { observableState } from '../../reactive/observable-state.js';
import type { ReadonlyState } from '../../reactive/readonly-state.js';
import { classes } from './_generated/measurements-add.module.css.js';

/** How many of the latest selections a measurement can use. */
const MaxPoints = 4;

/**
 * Adding measurements (the Selection tab's Measure section, after Mol*'s measurement controls): the points (the latest
 * selections, in order, the first four used; reorder, remove, a click focuses one), what can be measured from them
 * (each enabled once there are enough points), and the options.
 */
export function MeasurementsAdd(props: { model: MeasurementsModel }) {
  const { model } = props;
  const used = model.points.pipe(map((p) => p.slice(0, MaxPoints)));
  const older = model.points.pipe(
    map((p) => Math.max(0, p.length - MaxPoints)),
    distinctUntilChanged(),
  );
  const options = observableState(false);
  let list: HTMLElement | undefined;
  return (
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
          {(point, index) => (
            <PointRow
              model={model}
              point={point}
              index={index}
              onGrab={(e) => dragToReorder(e, list, point.id, (i) => model.movePointTo(point.id, i))}
            />
          )}
        </For>
        <Show when={older.pipe(map((n) => n > 0))}>
          {() => <p class={classes.more}>{older.pipe(map((n) => `+${n} older, not used`))}</p>}
        </Show>
      </div>
      <p class={classes.caption}>Measure</p>
      <div class={classes.kinds}>
        {model.kinds.map((kind) => (
          <KindRow model={model} kind={kind} />
        ))}
      </div>
      {/* The measurements made: adjust them here too (the same as their rows in the Structure tab). */}
      <Show when={model.added.pipe(map((a) => a.length > 0))}>
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
      <Show when={options}>{() => <ParamsForm form={model.options} />}</Show>
    </div>
  );
}

function PointRow(props: {
  model: MeasurementsModel;
  point: MeasurementPoint;
  index: ReadonlyState<number>;
  onGrab: (e: PointerEvent) => void;
}) {
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
        <span class={classes.index}>{props.index.pipe(map((i) => `${i + 1}`))}</span>
        <span class={classes.grip} aria-hidden="true" data-tip="Drag to reorder" onPointerDown={props.onGrab}>
          <Icon name="grip-vertical" size={13} />
        </span>
      </span>
      <button type="button" class={classes.pointLabel} title="Focus" onClick={() => model.focusPoint(point.id)}>
        {point.label}
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
  const ready = model.points.pipe(
    map((p) => (kind.points === 0 ? p.length > 0 : p.length >= kind.points)),
    distinctUntilChanged(),
  );
  const hint =
    kind.points === 0
      ? ready.pipe(map((r) => (r ? 'of the selection' : 'needs a selection')))
      : ready.pipe(
          map((r) => (r ? (kind.points === 1 ? 'point 1' : `points 1–${kind.points}`) : `needs ${kind.points} points`)),
        );
  return (
    <button type="button" class={classes.kind} disabled={ready.pipe(map((r) => !r))} onClick={() => model.add(kind.id)}>
      <Icon name={kind.icon} size={15} />
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
      <div class={classes.addedRow} data-hidden={node.isHidden.pipe(map((h) => (h ? '' : undefined)))}>
        <Icon name={node.icon} size={15} class={classes.addedIcon} />
        <span class={classes.pointLabel}>{node.label}</span>
        <span class={classes.addedValue}>{node.description}</span>
        <span class={classes.pointActions}>
          <button
            type="button"
            class={classes.addedAction}
            title={node.isHidden.pipe(map((h) => (h ? 'Show' : 'Hide')))}
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
        <Show when={open}>
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
