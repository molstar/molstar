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
import { effect, requireScope, runInScope } from '@molstar/prism-ui/dom/scope';
import { timeout } from '@molstar/prism-ui/dom/scoped';
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
  // Narrowed to these types (none pressed: every type the selection is drawn as)
  const only = observableState<readonly string[]>([]);
  const types = clipping.selectionTypes;
  const toggle = (t: string) =>
    only.set(only.value.includes(t) ? only.value.filter((x) => x !== t) : [...only.value, t]);
  const add = (shape: CutShape) => {
    const present = only.value.filter((t) => types.value.includes(t));
    void ctx.api.run('clipping.add', { shape, types: present });
  };
  const numbers = derive(clipping.cuts, (cuts) => cuts.map((c) => c.index), sameNumbers);
  return (
    <div style="display: grid; gap: var(--msp-space-3); min-width: 0">
      <p style={caption}>
        {derive(selected, (on) =>
          on ? 'Cut through what the selection is drawn as:' : 'Select something to cut through what it is drawn as.',
        )}
      </p>
      <div style={`${row}; flex-wrap: wrap`} role="group" aria-label="Cut with">
        {CutShapes.map((shape) => (
          <Button
            size="sm"
            icon={<Icon name={ShapeIcons[shape]} />}
            title={`Cut with a ${CutShapeLabels[shape].toLowerCase()}`}
            disabled={derive(selected, (on) => !on)}
            onClick={() => add(shape)}
          >
            {CutShapeLabels[shape]}
          </Button>
        ))}
      </div>
      <Show when={derive(types, (t) => t.length > 1)}>
        {() => (
          <div style={`${row}; flex-wrap: wrap`} role="group" aria-label="Only these representations">
            <span style={caption}>Only</span>
            <For each={types}>
              {(t) => (
                <Button
                  size="xs"
                  pressed={derive(only.asReadonly, (o) => o.includes(t))}
                  title={`Cut only the ${t} (none pressed: all)`}
                  onClick={() => toggle(t)}
                >
                  {t}
                </Button>
              )}
            </For>
          </div>
        )}
      </Show>
      <Show when={derive(numbers, (n) => n.length > 0)}>
        {() => (
          <div style="display: grid; gap: var(--msp-space-1)">
            <p style={caption}>Cuts</p>
            <For each={numbers}>
              {(index) => <CutRow ctx={ctx} cut={derive(clipping.cuts, (cs) => cs.find((c) => c.index === index))} />}
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
function CutRow(props: { ctx: ViewerUIContext; cut: ReadonlyState<Cut | undefined> }) {
  const { ctx, cut } = props;
  const editing = observableState(false);
  const run = (op: string, values: Record<string, unknown>) => {
    const c = cut.value;
    if (c) void ctx.api.run(op, { cut: c.index, ...values });
  };
  return (
    <Show when={cut}>
      {(c0) => {
        const shape = c0.object.type as CutShape;
        return (
          <div style="display: grid; gap: var(--msp-space-2); min-width: 0">
            <div style={row}>
              <Icon name={ShapeIcons[shape]} size={16} />
              <span style="flex: 1; min-width: 0; display: grid">
                <span>{derive(cut, (c) => (c ? `Cut ${c.index} · ${CutShapeLabels[shape]}` : ''))}</span>
                <span style={`${caption}; overflow: hidden; text-overflow: ellipsis; white-space: nowrap`}>
                  {derive(cut, (c) => (c ? c.targets.map((t) => t.label).join(', ') : ''))}
                </span>
              </span>
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
            <Show when={editing.asReadonly}>{() => <CutForm cut={cut} shape={shape} run={run} />}</Show>
          </div>
        );
      }}
    </Show>
  );
}

/** A cut's numbers by its shape: where, which way, how large (Å, degrees) */
function CutForm(props: {
  cut: ReadonlyState<Cut | undefined>;
  shape: CutShape;
  run(op: string, values: Record<string, unknown>): void;
}) {
  const { cut, shape } = props;
  const params = formParams(shape);
  const owner = requireScope('CutForm');
  let pending: (() => void) | undefined;
  const form = new ParamsFormModel(params, cut.value ? formValues(cut.value) : undefined, {
    // Typing and sliders change values at their rate: one update a moment after the last
    onChange: (values) => {
      pending?.();
      pending = runInScope(owner, () =>
        timeout(() => {
          pending = undefined;
          props.run('clipping.update', updateOf(shape, values));
        }, 150),
      );
    },
  });
  // The cut as it is now (a gizmo dragged, an undo): the form follows
  effect(cut, (c) => {
    if (c && !pending) form.load(formValues(c) as never);
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

function formValues(cut: Cut): Record<string, unknown> {
  const o = cut.object;
  const up = rotated(o.rotation, LocalUp);
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
      return { position: o.position, axis: rotated(o.rotation, ConeOpening), angle: coneAngle(o.scale) };
  }
  return {};
}

/** The form's values as `clipping.update`'s */
function updateOf(shape: CutShape, v: Record<string, any>): Record<string, unknown> {
  const position = [...v.position];
  switch (shape) {
    case 'plane':
      return { position, normal: [...v.normal] };
    case 'sphere':
      return { position, scale: [v.diameter, v.diameter, v.diameter] };
    case 'cube':
      return { position, normal: [...v.up], scale: [...v.size] };
    case 'cylinder':
      return { position, normal: [...v.axis], scale: [v.diameter, v.height, v.diameter] };
    case 'infiniteCone':
      return { position, normal: [...v.axis], scale: [...coneScale(v.angle)] };
  }
}
