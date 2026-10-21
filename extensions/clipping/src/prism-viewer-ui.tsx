/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Vec3 } from '@molstar/core/math/linear-algebra';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Button, domView, Icon, ParamsForm } from '@molstar/prism-ui';
import type { IconName } from '@molstar/prism-ui/components/icon/icon';
import { For, Show } from '@molstar/prism-ui/dom/control-flow';
import type { ReadonlyState } from '@molstar/prism-ui/core';
import { effect, onActive, onCleanup } from '@molstar/prism-ui/dom/scope';
import { ParamsFormModel } from '@molstar/prism-ui/forms/params-form';
import { observableState } from '@molstar/prism-ui/reactive/observable-state';
import { derive } from '@molstar/prism-ui/reactive/state';
import type { ViewerExtensionUI, ViewerUIContext } from '@molstar/prism-viewer/spec';
import {
  ConeOpening,
  type Cut,
  type CutShape,
  CutShapeLabels,
  CutShapes,
  coneAngle,
  coneScale,
  LocalUp,
  rotated,
} from './cuts.js';
import type { ClippingApi } from './operations.js';
import { ClippingSectionId } from './prism-viewer.js';

// The Clipping section of the Selection tab (_spec/README.md §3): a cut through what the selection is drawn as, by
// shape, narrowed to some representation types; the cuts, each flipped, edited or removed. It acts only through the
// person's API (`clipping.*`) and reads the cuts the area follows.

export const ShapeIcons: Readonly<Record<CutShape, IconName>> = {
  plane: 'square',
  sphere: 'circle',
  cube: 'box',
  cylinder: 'cylinder',
  infiniteCone: 'cone',
};

export function clippingUI(ctx: ViewerUIContext): ViewerExtensionUI {
  return {
    sections: [
      {
        id: ClippingSectionId,
        slot: 'selection',
        title: 'Clipping',
        icon: 'scissors',
        order: 0.5,
        keywords: ['clip', 'cut', 'slice', 'cross section', 'plane', 'inside'],
        render: domView(() => {
          const clipping = ctx.api.clipping;
          return clipping ? <ClippingSection ctx={ctx} clipping={clipping} /> : <></>;
        }),
      },
    ],
    nodeActions: [
      {
        id: 'clipping.node-action.clip-this',
        title: 'Clip this…',
        icon: 'scissors',
        order: 10,
        kinds: ['component', 'representation'],
        when: (node) => !!ctx.api.clipping?.canClipThis(node.id),
        run: (node) => {
          const clipping = ctx.api.clipping;
          if (!clipping) return;
          void clipping.clipThis(ctx.api, node.id).then((ok) => {
            if (ok) void ctx.api.view.open({ tab: 'selection', section: ClippingSectionId });
          });
        },
      },
    ],
    selectionActions: [
      {
        id: 'clipping.selection-action.clip',
        title: 'Clip…',
        hint: 'Cut through what it is drawn as',
        icon: 'scissors',
        order: 10,
        run: () => void ctx.api.view.open({ tab: 'selection', section: ClippingSectionId }),
      },
    ],
  };
}

const row = 'display: flex; align-items: center; gap: var(--msp-space-2); min-width: 0';
const caption = 'margin: 0; font-size: var(--msp-font-size-xs); color: var(--msp-color-text-muted)';

function ClippingSection(props: { ctx: ViewerUIContext; clipping: ClippingApi }) {
  const { ctx, clipping } = props;
  const selected = derive(ctx.sources.interaction.selection, (s) => s !== undefined);
  // The types it cuts: every type the selection is drawn as but those unchecked (a new type comes checked)
  const types = clipping.selectionTypes;
  const filter = clipping.filter;
  const chosen = derive([types, filter.asReadonly], (ts, f) =>
    ts.filter((t) => (f.only ? f.only.includes(t) : !f.excluded.includes(t))),
  );
  // Set by hand, the filter is the types left out (a type drawn later comes checked)
  const exclude = (excluded: readonly string[]) => filter.set({ excluded });
  const toggle = (t: string) => {
    const on = chosen.value.includes(t);
    exclude(types.value.filter((x) => (x === t ? on : !chosen.value.includes(x))));
  };
  const add = (shape: CutShape) => {
    // all chosen: no filter, so a type drawn later is cut too
    const only = chosen.value.length === types.value.length ? [] : chosen.value;
    void ctx.api.run('clipping.add', { shape, types: only }).then((r) => {
      // a cut made here comes with its handles
      const cut = (r as { value?: Cut }).value;
      if (r.ok && cut) clipping.gizmos.setHandled(cut.index, true);
    });
  };
  // The handled cuts' gizmos show while the section does: closed (or its tab), they go; opened, they come back
  onActive((active) => (active ? clipping.showGizmos(ctx.api) : clipping.gizmos.hide()));
  onCleanup(() => clipping.gizmos.hide());
  const canAdd = derive([selected, chosen], (on, c) => on && c.length > 0);
  const numbers = derive(clipping.cuts, (cuts) => cuts.map((c) => c.index), sameNumbers);
  return (
    <div style="display: grid; gap: var(--msp-space-3); min-width: 0">
      <p style={caption}>
        {derive(selected, (on) =>
          on ? 'Cut through what the selection is drawn as:' : 'Select something to cut through what it is drawn as.',
        )}
      </p>
      {/* As a representation's Type: a strip of icons, their names as tooltips (a shape is plain from its icon) */}
      <div
        style="display: flex; padding: 2px; border-radius: var(--msp-radius-md); background: var(--msp-color-control); min-width: 0"
        role="group"
        aria-label="Cut with"
      >
        {CutShapes.map((shape) => (
          <span style="flex: 1 1 0; min-width: 0; display: grid">
            {/* the icon as content, not `icon`: an icon-only button is a fixed square, this one fills its slot */}
            <Button
              size="sm"
              variant="ghost"
              block
              title={`Cut with a ${CutShapeLabels[shape].toLowerCase()}`}
              disabled={derive(canAdd, (on) => !on)}
              onClick={() => add(shape)}
            >
              <Icon name={ShapeIcons[shape]} label={CutShapeLabels[shape]} />
            </Button>
          </span>
        ))}
      </div>
      <Show when={derive(types, (t) => t.length > 1)}>
        {() => (
          <div style="display: flex; flex-wrap: wrap; align-items: center; gap: var(--msp-space-1); min-width: 0">
            <span style={`${caption}; margin-right: var(--msp-space-1)`}>Through</span>
            <div style="display: contents" role="group" aria-label="Cut through these representations">
              <For each={types}>
                {(t) => {
                  const on = derive(chosen, (c) => c.includes(t));
                  return (
                    <Button
                      size="xs"
                      variant="ghost"
                      pressed={on}
                      icon={
                        <span style={derive(on, (v) => (v ? 'display: flex' : 'display: flex; visibility: hidden'))}>
                          <Icon name="check" size={14} />
                        </span>
                      }
                      title={derive(on, (v) => (v ? `Leave the ${t} whole` : `Cut through the ${t} too`))}
                      onClick={() => toggle(t)}
                    >
                      <span style={derive(on, (v) => (v ? '' : 'color: var(--msp-color-text-muted)'))}>{t}</span>
                    </Button>
                  );
                }}
              </For>
            </div>
            <span style="flex: 1" />
            <Button size="xs" variant="ghost" title="Cut through every type" onClick={() => exclude([])}>
              All
            </Button>
            <Button size="xs" variant="ghost" title="Uncheck every type" onClick={() => exclude(types.value)}>
              None
            </Button>
          </div>
        )}
      </Show>
      <Show when={derive(numbers, (n) => n.length > 0)}>
        {() => (
          <div style="display: grid; gap: var(--msp-space-1)">
            <div style={row}>
              <p style={`${caption}; flex: 1`}>Cuts</p>
              <Button
                size="xs"
                variant="ghost"
                title="Handles in 3D on every cut: drag to move, turn or size it (Esc cancels)"
                onClick={() => clipping.gizmos.handleAll(true)}
              >
                All
              </Button>
              <Button size="xs" variant="ghost" title="No handles" onClick={() => clipping.gizmos.handleAll(false)}>
                None
              </Button>
            </div>
            <For each={numbers}>
              {(index) => (
                <CutRow
                  ctx={ctx}
                  clipping={clipping}
                  cut={derive(clipping.cuts, (cs) => cs.find((c) => c.index === index))}
                />
              )}
            </For>
          </div>
        )}
      </Show>
    </div>
  );
}

function sameNumbers(a: readonly number[], b: readonly number[]) {
  return a.length === b.length && a.every((x, i) => x === b[i]);
}

/** A cut: its shape and what it is on; flip, edit (its numbers), remove */
function CutRow(props: { ctx: ViewerUIContext; clipping: ClippingApi; cut: ReadonlyState<Cut | undefined> }) {
  const { ctx, clipping, cut } = props;
  const handled = derive([clipping.gizmos.handled, cut], (h, c) => !!c && h.includes(c.index));
  const editing = observableState(false);
  const run = (op: string, values: Record<string, unknown>) => {
    const c = cut.value;
    return c ? ctx.api.run(op, { cut: c.index, ...values }) : undefined;
  };
  return (
    // keyed by its shape (a string, stable): a new `Cut` each state change must not render the row (and a form being
    // edited) again
    <Show when={derive(cut, (c) => c?.object.type)}>
      {(type) => {
        const shape = type as CutShape;
        return (
          <div style="display: grid; gap: var(--msp-space-2); min-width: 0">
            <div style={row}>
              <Icon name={ShapeIcons[shape]} size={16} />
              <span style="flex: 1; min-width: 0; display: grid">
                <span style="overflow: hidden; text-overflow: ellipsis; white-space: nowrap">
                  {derive(cut, (c) => (c ? `${CutShapeLabels[shape]} ${c.index}` : ''))}
                </span>
                <span style={`${caption}; overflow: hidden; text-overflow: ellipsis; white-space: nowrap`}>
                  {derive(cut, (c) => (c ? c.targets.map((t) => t.label).join(', ') : ''))}
                </span>
              </span>
              <Button
                size="xs"
                variant="ghost"
                icon={<Icon name="move-3d" />}
                title={derive(handled, (on) =>
                  on ? 'Hide its handles' : 'Handles in 3D: drag to move, turn or size it (Esc cancels)',
                )}
                pressed={handled}
                onClick={() => {
                  const c = cut.value;
                  if (c) clipping.gizmos.setHandled(c.index, !handled.value);
                }}
              />
              <Button
                size="xs"
                variant="ghost"
                icon={<Icon name="flip-vertical-2" />}
                title={shape === 'plane' ? 'Cut away the other side' : 'Cut away the outside instead'}
                pressed={derive(cut, (c) => !!c?.object.invert)}
                onClick={() => run('clipping.update', { invert: !cut.value?.object.invert })}
              />
              <Button
                size="xs"
                variant="ghost"
                icon={<Icon name="sliders-horizontal" />}
                title="Edit its numbers"
                pressed={editing.asReadonly}
                onClick={() => editing.set(!editing.value)}
              />
              <Button
                size="xs"
                variant="ghost"
                icon={<Icon name="trash-2" />}
                title="Remove the cut"
                onClick={() => run('clipping.remove', {})}
              />
            </div>
            <Show when={editing.asReadonly}>
              {() => (
                <>
                  <CutTargets clipping={clipping} cut={cut} run={run} />
                  <CutForm cut={cut} shape={shape} run={run} />
                </>
              )}
            </Show>
          </div>
        );
      }}
    </Show>
  );
}

/** What a cut is on, and the representations near them it could be: each checked or not (the last stays checked) */
function CutTargets(props: {
  clipping: ClippingApi;
  cut: ReadonlyState<Cut | undefined>;
  run(op: string, values: Record<string, unknown>): void;
}) {
  const candidates = derive(props.cut, (c) => (c ? props.clipping.candidates(c) : []), sameTargets);
  return (
    <div style="display: grid; gap: 2px; min-width: 0" role="group" aria-label="The representations it cuts">
      <p style={caption}>Cuts through</p>
      <For each={candidates}>
        {(t) => {
          const on = derive(props.cut, (c) => !!c?.targets.some((x) => x.ref === t.ref));
          const last = derive(props.cut, (c) => !!c && c.targets.length === 1 && c.targets[0].ref === t.ref);
          return (
            <Button
              size="xs"
              variant="ghost"
              block
              pressed={on}
              disabled={last}
              icon={
                <span style={derive(on, (v) => (v ? 'display: flex' : 'display: flex; visibility: hidden'))}>
                  <Icon name="check" size={14} />
                </span>
              }
              title={derive([on, last], (v, l) =>
                l
                  ? 'The cut is only on this: remove the cut instead'
                  : v
                    ? `Leave ${t.label} whole`
                    : `Cut ${t.label} too`,
              )}
              onClick={() => props.run('clipping.apply', { representation: t.ref, on: !on.value })}
            >
              <span
                style={derive(
                  on,
                  (v) =>
                    `flex: 1; text-align: start; min-width: 0; overflow: hidden; text-overflow: ellipsis; white-space: nowrap${v ? '' : '; color: var(--msp-color-text-muted)'}`,
                )}
              >
                {t.label}
              </span>
            </Button>
          );
        }}
      </For>
    </div>
  );
}

function sameTargets(a: readonly { ref: string }[], b: readonly { ref: string }[]) {
  return a.length === b.length && a.every((t, i) => t.ref === b[i].ref);
}

/** A cut's numbers by its shape: where, which way, how large (Å, degrees) */
function CutForm(props: {
  cut: ReadonlyState<Cut | undefined>;
  shape: CutShape;
  run(op: string, values: Record<string, unknown>): Promise<unknown> | void;
}) {
  const { cut, shape } = props;
  const params = formParams(shape);
  // An edit applies at once, one update in flight: fields changed meanwhile go in the next (the latest values win).
  // Only the changed fields are sent: the others as shown are rounded, and a box's up alone would drop its twist.
  let running = false;
  let latest: Record<string, unknown> | undefined;
  const changed = new Set<string>();
  const flush = async () => {
    running = true;
    while (latest) {
      const values = latest;
      const fields = new Set(changed);
      latest = undefined;
      changed.clear();
      await props.run('clipping.update', updateOf(shape, values, fields));
    }
    running = false;
  };
  const form = new ParamsFormModel(params, cut.value ? formValues(cut.value) : undefined, {
    onChange: (values, name) => {
      latest = values;
      changed.add(name);
      if (!running) void flush();
    },
  });
  // The cut as it is now (a gizmo dragged, an undo): the form follows
  effect(cut, (c) => {
    if (c && !running) form.load(formValues(c) as never);
  });
  return <ParamsForm form={form} />;
}

const Position = PD.Vec3(Vec3(), {}, { description: 'Its center (Å); a plane goes through it, a cone starts at it' });
const Direction = (description: string) => PD.Vec3(Vec3.create(0, 1, 0), {}, { description });
const Length = (label: string) => PD.Numeric(10, { min: 0.1, max: 2000, step: 0.1 }, { label });

function formParams(shape: CutShape): PD.Params {
  switch (shape) {
    case 'plane':
      return { position: Position, normal: Direction('It cuts away the side this points to') };
    case 'sphere':
      return { position: Position, diameter: Length('Diameter (Å)') };
    case 'cube':
      return { position: Position, up: Direction('Its height axis'), size: PD.Vec3(Vec3.create(10, 10, 10)) };
    case 'cylinder':
      return {
        position: Position,
        axis: Direction('Its axis'),
        diameter: Length('Diameter (Å)'),
        height: Length('Height (Å)'),
      };
    case 'infiniteCone':
      return {
        position: Position,
        axis: Direction('Which way it opens'),
        angle: PD.Numeric(25, { min: 1, max: 89, step: 1 }, { label: 'Half-angle (°)' }),
      };
  }
}

/** Shown to 2 decimals (Å), directions to 3: the fields fit, and nothing finer is meant */
const round = (v: ArrayLike<number>, digits: number) =>
  Vec3.create(...([0, 1, 2].map((i) => +v[i].toFixed(digits)) as [number, number, number]));

function formValues(cut: Cut): Record<string, unknown> {
  const raw = cut.object;
  const o = { ...raw, position: round(raw.position, 2), scale: round(raw.scale, 2) };
  const up = round(rotated(raw.rotation, LocalUp), 3);
  switch (o.type as CutShape) {
    case 'plane':
      return { position: o.position, normal: up };
    case 'sphere':
      return { position: o.position, diameter: o.scale[0] };
    case 'cube':
      return { position: o.position, up, size: o.scale };
    case 'cylinder':
      return { position: o.position, axis: up, diameter: o.scale[0], height: o.scale[1] };
    case 'infiniteCone':
      return {
        position: o.position,
        axis: round(rotated(raw.rotation, ConeOpening), 3),
        angle: +coneAngle(raw.scale).toFixed(1),
      };
  }
  return {};
}

/** The changed fields of the form as `clipping.update`'s values */
function updateOf(shape: CutShape, v: Record<string, any>, fields: ReadonlySet<string>): Record<string, unknown> {
  const out: Record<string, unknown> = {};
  if (fields.has('position')) out.position = [...v.position];
  const direction = { plane: 'normal', cube: 'up', cylinder: 'axis', infiniteCone: 'axis', sphere: '' }[shape];
  if (direction && fields.has(direction)) out.normal = [...v[direction]];
  switch (shape) {
    case 'sphere':
      if (fields.has('diameter')) out.scale = [v.diameter, v.diameter, v.diameter];
      break;
    case 'cube':
      if (fields.has('size')) out.scale = [...v.size];
      break;
    case 'cylinder':
      if (fields.has('diameter') || fields.has('height')) out.scale = [v.diameter, v.height, v.diameter];
      break;
    case 'infiniteCone':
      if (fields.has('angle')) out.scale = [...coneScale(v.angle)];
      break;
  }
  return out;
}
