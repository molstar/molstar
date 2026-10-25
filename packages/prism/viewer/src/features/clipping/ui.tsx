/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Ludovic Autin <ludovic.autin@gmail.com>
 */

import { Vec3 } from '@molstar/core/math/linear-algebra';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Button, domView, Icon, ParamsForm, TextInput } from '@molstar/prism-ui';
import type { IconName } from '@molstar/prism-ui/components/icon/icon';
import { For, Show } from '@molstar/prism-ui/dom/control-flow';
import type { ReadonlyState } from '@molstar/prism-ui/core';
import { effect, onActive, onCleanup } from '@molstar/prism-ui/dom/scope';
import { ParamsFormModel } from '@molstar/prism-ui/forms/params-form';
import { observableState } from '@molstar/prism-ui/reactive/observable-state';
import { derive } from '@molstar/prism-ui/reactive/state';
import type { ViewerExtensionUI, ViewerUIContext } from '../../spec.js';
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
} from '../../sources/clipping/cuts.js';
import type { ClippingApi, ClipSettings } from './operations.js';
import { classes } from './_generated/clipping.module.css.js';
import { ClippingSectionId } from './index.js';

// The Clipping section of the Selection tab (clipping.md §3): a cut through what the selection is drawn as, by
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
        kinds: [
          'component',
          'representation',
          'volume',
          'volume-representation',
          'particles',
          'particles-representation',
        ],
        when: (node) => ctx.api.clipping.canClipThis(node.id),
        run: (node) => {
          const clipping = ctx.api.clipping;
          void clipping.clipThis(ctx.api, node.id).then((done) => {
            if (!done) return;
            // a cut made here (a volume's, a shape's): its handles on, its editor the one open
            if (done.cut) {
              clipping.gizmos.handleOnly(done.cut.index);
            }
            void ctx.api.view.open({ tab: 'selection', section: ClippingSectionId });
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
  // One cut is active at a time: its editor open, its handles shown (so a cut drawn in the scene is never one that can
  // not be dragged); the gizmo controller keeps it, following it when the cuts are renumbered
  const active = clipping.gizmos.handled;
  const setOpen = (index: number, on: boolean) => clipping.gizmos.handleOnly(on ? index : undefined);
  // + Add: a plane through what the selection is drawn as; its shape is picked on the cut
  const add = () => {
    // all chosen: no filter, so a type drawn later is cut too
    const only = chosen.value.length === types.value.length ? [] : chosen.value;
    // nothing selected (a map alone, say): everything drawn, narrowed on the cut afterwards
    const values = selected.value
      ? { shape: 'plane', types: only, ...(clipping.keepNew.value ? { keepWhole: true } : {}) }
      : { shape: 'plane', everything: true };
    void ctx.api.run('clipping.add', values).then((r) => {
      // a new cut, or the one like it already there (Add pressed again): either is what is worked on now
      const cut = (r as { value?: Cut }).value;
      if (!r.ok || !cut) return;
      // the new cut (or the one like it) is the active one: its editor open, its handles shown
      clipping.gizmos.handleOnly(cut.index);
    });
  };
  // The handled cuts' gizmos show while the section does: closed (or its tab), they go; opened, they come back
  onActive((active) => (active ? clipping.showGizmos(ctx.api) : clipping.gizmos.hide()));
  onCleanup(() => clipping.gizmos.hide());
  const canAdd = derive([selected, chosen], (on, c) => !on || c.length > 0);
  const numbers = derive(clipping.cuts, (cuts) => cuts.map((c) => c.index), sameNumbers);
  const any = derive(numbers, (n) => n.length > 0);
  return (
    <div style="display: grid; gap: var(--msp-space-2); min-width: 0">
      <div style={row}>
        <p style={`${caption}; flex: 1`}>Cuts</p>
        <Button
          size="xs"
          variant="ghost"
          icon={<Icon name="plus" />}
          title={derive(selected, (on) =>
            on
              ? 'A cut through what the selection is drawn as (pick its shape on it)'
              : 'A cut through everything drawn (narrow it on the cut)',
          )}
          disabled={derive(canAdd, (on) => !on)}
          onClick={add}
        >
          Add
        </Button>
      </div>
      <Show when={derive([any, selected], (a, s) => !a && !s)}>
        {() => <p style={caption}>Add cuts what the selection is drawn as; with nothing selected, everything drawn.</p>}
      </Show>
      <Show when={derive(types, (t) => t.length > 1)}>
        {() => (
          <div style="display: flex; flex-wrap: wrap; align-items: center; gap: var(--msp-space-1); min-width: 0">
            <span style={`${caption}; margin-right: var(--msp-space-1)`}>Through</span>
            <div style="display: contents" role="group" aria-label="New cuts through these representation types">
              <For each={types}>
                {(t) => {
                  const on = derive(chosen, (c) => c.includes(t));
                  return (
                    <Button
                      size="xs"
                      variant="ghost"
                      pressed={on}
                      icon={<Check on={on} />}
                      title={derive(on, (v) => (v ? `Leave the ${t} whole` : `Cut through the ${t} too`))}
                      onClick={() => toggle(t)}
                    >
                      <Muted on={on}>{t}</Muted>
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
      <Show when={selected}>
        {() => {
          const on = clipping.keepNew.asReadonly;
          return (
            <Button
              size="xs"
              variant="ghost"
              block
              pressed={on}
              icon={<Check on={on} />}
              title="A cut added keeps the selected atoms whole: it goes around them"
              onClick={() => clipping.keepNew.set(!on.value)}
            >
              <Muted on={on}>Keep the selection whole</Muted>
            </Button>
          );
        }}
      </Show>
      <For each={numbers}>
        {(index) => (
          <CutRow
            ctx={ctx}
            clipping={clipping}
            cut={derive(clipping.cuts, (cs) => cs.find((c) => c.index === index))}
            open={derive(active, (a) => a.includes(index))}
            setOpen={(on) => setOpen(index, on)}
          />
        )}
      </For>
      <Fold title="Settings" meta={derive(clipping.cutSettings, settingsMeta)}>
        {() => <CutSettings ctx={ctx} clipping={clipping} />}
      </Fold>
    </div>
  );
}

/** "· whole copies" while cuts leave out whole instances; nothing for the usual pixel cuts */
const settingsMeta = (s: ClipSettings) => (s.variant === 'instance' ? '· whole copies' : '');

/** How every cut cuts (clipping.settings): the variant, the solid caps; following them as they change elsewhere */
function CutSettings(props: { ctx: ViewerUIContext; clipping: ClippingApi }) {
  const { ctx, clipping } = props;
  // the operation's own params (each optional there: absent, kept), labeled as a form's rows
  const op = clipping.settings.definition.params as PD.Params;
  const params = {
    variant: { ...op.variant, isOptional: false, label: 'Cut by' },
    solidSurface: { ...op.solidSurface, isOptional: false, label: 'Solid surface' },
  } as PD.Params;
  const form = new ParamsFormModel(params, clipping.cutSettings.value as never, {
    onChange: (values, name) => void ctx.api.run('clipping.settings', { [name]: values[name] }),
  });
  effect(clipping.cutSettings, (s) => {
    const v = form.values.value;
    if (v.variant !== s.variant || v.solidSurface !== s.solidSurface) form.load(s as never);
  });
  return <ParamsForm form={form} />;
}

function sameNumbers(a: readonly number[], b: readonly number[]) {
  return a.length === b.length && a.every((x, i) => x === b[i]);
}

/** A check mark, its place kept when off (the labels line up) */
function Check(props: { on: ReadonlyState<boolean>; partly?: ReadonlyState<boolean> }) {
  const partly = props.partly;
  return (
    <span style={derive(props.on, (v) => (v ? 'display: flex' : 'display: flex; visibility: hidden'))}>
      {partly ? (
        <Show when={partly} fallback={() => <Icon name="check" size={14} />}>
          {() => <Icon name="minus" size={14} />}
        </Show>
      ) : (
        <Icon name="check" size={14} />
      )}
    </span>
  );
}

/** A label, muted while its check is off */
function Muted(props: { on: ReadonlyState<boolean>; children?: unknown }) {
  return (
    <span
      class={classes.fade}
      style={derive(
        props.on,
        (v) => `flex: 1; text-align: start; min-width: 0${v ? '' : '; color: var(--msp-color-text-muted)'}`,
      )}
    >
      {props.children as never}
    </span>
  );
}

/** A fold in a cut's editor: its title and a count, the chevron on the right; closed at first */
function Fold(props: { title: string; meta?: ReadonlyState<string>; children: () => unknown }) {
  const open = observableState(false);
  return (
    <div style="display: grid; gap: 2px; min-width: 0">
      <button
        type="button"
        class={classes.fold}
        aria-expanded={derive(open.asReadonly, (o) => String(o))}
        onClick={() => open.set(!open.value)}
      >
        <span class={`${classes.foldTitle} ${classes.fade}`}>
          {props.title}
          {props.meta ? <> {props.meta}</> : null}
        </span>
        <Icon name="chevron-right" size={14} class={classes.foldChevron} />
      </button>
      <Show when={open.asReadonly}>{() => props.children() as never}</Show>
    </div>
  );
}

/** The five shapes as a strip of icons (as a representation's Type), the current one pressed */
function ShapeStrip(props: { current: CutShape; pick(shape: CutShape): void }) {
  return (
    <div
      style="display: flex; padding: 2px; border-radius: var(--msp-radius-md); background: var(--msp-color-control); min-width: 0"
      role="radiogroup"
      aria-label="Shape"
    >
      {CutShapes.map((shape) => (
        <span style="flex: 1 1 0; min-width: 0; display: grid">
          {/* the icon as content, not `icon`: an icon-only button is a fixed square, this one fills its slot */}
          <Button
            size="sm"
            variant="ghost"
            block
            pressed={shape === props.current}
            title={
              shape === props.current
                ? `A ${CutShapeLabels[shape].toLowerCase()}`
                : `Make it a ${CutShapeLabels[shape].toLowerCase()}`
            }
            onClick={() => shape !== props.current && props.pick(shape)}
          >
            <Icon name={ShapeIcons[shape]} label={CutShapeLabels[shape]} />
          </Button>
        </span>
      ))}
    </div>
  );
}

/** A cut: its shape and what it is on; handles, flip, edit (its shape, what it cuts, its numbers), remove */
function CutRow(props: {
  ctx: ViewerUIContext;
  clipping: ClippingApi;
  cut: ReadonlyState<Cut | undefined>;
  open: ReadonlyState<boolean>;
  setOpen(on: boolean): void;
}) {
  const { ctx, clipping, cut } = props;
  const run = (op: string, values: Record<string, unknown>) => {
    const c = cut.value;
    // by its tag where it has one: a run waits its turn (a double click's second, behind its first), and a number
    // shifts when a cut before it goes
    return c ? ctx.api.run(op, { ...(c.object.tag ? { tag: c.object.tag } : { cut: c.index }), ...values }) : undefined;
  };
  return (
    // keyed by its shape (a string, stable): a new `Cut` each state change must not render the row (and a form being
    // edited) again; another shape is another form
    <Show when={derive(cut, (c) => c?.object.type)}>
      {(type) => {
        const shape = type as CutShape;
        return (
          <div style="display: grid; gap: var(--msp-space-2); min-width: 0">
            <div class={classes.cutRow}>
              {/* its icon and name open and close its editor */}
              <button
                type="button"
                class={classes.cutLabel}
                aria-expanded={derive(props.open, (o) => String(o))}
                title={derive(props.open, (o) =>
                  o ? 'Close its editor' : 'Its shape, what it cuts, where, how drawn',
                )}
                onClick={() => props.setOpen(!props.open.value)}
              >
                <span class={classes.cutIcon}>
                  <Icon name={ShapeIcons[shape]} size={18} />
                </span>
                <span class={classes.cutText}>
                  <span class={classes.fade}>
                    {derive(cut, (c) => (c ? `${CutShapeLabels[shape]} ${c.index}` : ''))}
                  </span>
                  <span class={`${classes.caption} ${classes.fade}`}>
                    {derive(cut, (c) => (c ? targetsSummary(c) : ''))}
                  </span>
                </span>
              </button>
              <Button
                size="xs"
                variant="ghost"
                icon={<Icon name="scan" />}
                title={derive(cut, (c) =>
                  c?.drawing ? 'Drawn in the scene (saved with it): hide it' : 'Draw it in the scene (saved with it)',
                )}
                pressed={derive(cut, (c) => !!c?.drawing)}
                onClick={() => run('clipping.show', { on: !cut.value?.drawing })}
              />
              <Button
                size="xs"
                variant="ghost"
                icon={<Icon name="trash-2" />}
                title="Remove the cut"
                onClick={() => run('clipping.remove', {})}
              />
            </div>
            <Show when={props.open}>
              {() => (
                <>
                  <ShapeStrip current={shape} pick={(next) => run('clipping.update', { shape: next })} />
                  <Button
                    size="xs"
                    variant="ghost"
                    block
                    pressed={derive(cut, (c) => !!c?.object.invert)}
                    icon={<Icon name="flip-vertical-2" />}
                    title={shape === 'plane' ? 'Which side it cuts away' : 'Its inside cut away, or its outside'}
                    onClick={() => run('clipping.update', { invert: !cut.value?.object.invert })}
                  >
                    <span style="flex: 1; text-align: start">
                      {shape === 'plane' ? 'Cut away the other side' : 'Cut away the outside instead'}
                    </span>
                  </Button>
                  <CutTargets clipping={clipping} cut={cut} run={run} />
                  <Fold title="Kept whole" meta={derive(cut, (c) => (c?.kept ? `(${c.kept.label})` : '(nothing)'))}>
                    {() => <KeptWhole ctx={ctx} cut={cut} run={run} />}
                  </Fold>
                  <Fold title="Position and size">{() => <CutForm cut={cut} shape={shape} run={run} />}</Fold>
                  <Show when={derive(cut, (c) => !!c?.drawing)}>
                    {() => (
                      <Fold title="In the scene">{() => <DrawingForm clipping={clipping} cut={cut} run={run} />}</Fold>
                    )}
                  </Show>
                </>
              )}
            </Show>
          </div>
        );
      }}
    </Show>
  );
}

/** "Cartoon of Polymer, Ball & Stick of Ligand", or "12 representations" past a few */
function targetsSummary(cut: Cut): string {
  return cut.targets.length > 3 ? `${cut.targets.length} representations` : cut.targets.map((t) => t.label).join(', ');
}

/** Past this many representations, the list is searched and shows this many until asked for all */
const TargetsShown = 8;

/**
 * What a cut is on, and the representations near them it could be: by type (each on or off for all of its kind, one
 * step), then one by one, those it is on first; past a few, a search and "Show all". Turning off the last says to
 * remove the cut instead (the operation refuses).
 */
function CutTargets(props: {
  clipping: ClippingApi;
  cut: ReadonlyState<Cut | undefined>;
  run(op: string, values: Record<string, unknown>): Promise<unknown> | void;
}) {
  const candidates = derive(props.cut, (c) => (c ? props.clipping.candidates(c) : []), sameTargets);
  const on = derive(props.cut, (c) => new Set(c?.targets.map((t) => t.ref)));
  const byType = derive([candidates, on], (cs, o) => {
    const types = new Map<string, { refs: string[]; on: number }>();
    for (const t of cs) {
      const e = types.get(t.type) ?? { refs: [], on: 0 };
      e.refs.push(t.ref);
      if (o.has(t.ref)) e.on++;
      types.set(t.type, e);
    }
    return [...types].map(([type, e]) => ({ type, ...e }));
  });
  const query = observableState('');
  const all = observableState(false);
  const many = derive(candidates, (cs) => cs.length > TargetsShown);
  const shown = derive(
    [candidates, on, query.asReadonly, all.asReadonly],
    (cs, o, q, a) => {
      const needle = q.trim().toLowerCase();
      const list = [...cs]
        .filter((t) => !needle || t.label.toLowerCase().includes(needle))
        .sort((x, y) => Number(o.has(y.ref)) - Number(o.has(x.ref)));
      return a || needle ? list : list.slice(0, TargetsShown);
    },
    sameTargets,
  );
  const apply = (refs: readonly string[], value: boolean) =>
    props.run('clipping.apply', { representations: refs, on: value });
  const count = derive([candidates, on], (cs, o) => `(${cs.filter((t) => o.has(t.ref)).length}/${cs.length})`);
  return (
    <Fold title="Cuts through" meta={count}>
      {() => (
        <div style="display: grid; gap: 2px; min-width: 0" role="group" aria-label="The representations it cuts">
          <Show when={derive([byType, candidates], (t, cs) => t.length > 1 && cs.length > t.length)}>
            {() => (
              <div style="display: flex; flex-wrap: wrap; gap: var(--msp-space-1); min-width: 0" aria-label="By type">
                <For each={derive(byType, (t) => t.map((x) => x.type), sameNumbersOrStrings)}>
                  {(type) => {
                    const entry = derive(byType, (t) => t.find((x) => x.type === type));
                    const any = derive(entry, (e) => !!e && e.on > 0);
                    const partly = derive(entry, (e) => !!e && e.on > 0 && e.on < e.refs.length);
                    return (
                      <Button
                        size="xs"
                        variant="ghost"
                        pressed={any}
                        icon={<Check on={any} partly={partly} />}
                        title={derive(entry, (e) =>
                          e && e.on === e.refs.length ? `Leave every ${type} whole` : `Cut every ${type}`,
                        )}
                        onClick={() => {
                          const e = entry.value;
                          if (e) apply(e.refs, e.on < e.refs.length);
                        }}
                      >
                        <Muted on={any}>{derive(entry, (e) => (e ? `${type} ${e.on}/${e.refs.length}` : type))}</Muted>
                      </Button>
                    );
                  }}
                </For>
              </div>
            )}
          </Show>
          <Show when={many}>
            {() => (
              <TextInput
                size="sm"
                search
                placeholder="Find a representation"
                aria-label="Find a representation"
                value={query}
                onInput={(v) => query.set(v)}
              />
            )}
          </Show>
          <For each={shown} by={(t) => t.ref}>
            {(t) => {
              const checked = derive(on, (o) => o.has(t.ref));
              return (
                <Button
                  size="xs"
                  variant="ghost"
                  block
                  pressed={checked}
                  icon={<Check on={checked} />}
                  title={derive(checked, (v) => (v ? `Leave ${t.label} whole` : `Cut ${t.label} too`))}
                  onClick={() => apply([t.ref], !checked.value)}
                >
                  <Muted on={checked}>{t.label}</Muted>
                </Button>
              );
            }}
          </For>
          <Show when={derive([many, query.asReadonly], (m, q) => m && !q.trim())}>
            {() => (
              <Button size="xs" variant="ghost" onClick={() => all.set(!all.value)}>
                {derive([all.asReadonly, candidates], (a, cs) => (a ? 'Show fewer' : `Show all ${cs.length}`))}
              </Button>
            )}
          </Show>
        </div>
      )}
    </Fold>
  );
}

function sameTargets(a: readonly { ref: string }[], b: readonly { ref: string }[]) {
  return a.length === b.length && a.every((t, i) => t.ref === b[i].ref);
}

function sameNumbersOrStrings<T extends string | number>(a: readonly T[], b: readonly T[]) {
  return a.length === b.length && a.every((x, i) => x === b[i]);
}

/** What a cut keeps whole: the atoms (in words), cut again; the selection added */
function KeptWhole(props: {
  ctx: ViewerUIContext;
  cut: ReadonlyState<Cut | undefined>;
  run(op: string, values: Record<string, unknown>): Promise<unknown> | void;
}) {
  const { cut } = props;
  const selected = derive(props.ctx.sources.interaction.selection, (s) => s !== undefined);
  return (
    <div style="display: grid; gap: 2px; min-width: 0">
      <Show when={derive(cut, (c) => c?.kept?.label)}>
        {(label) => (
          <div class={classes.cutRow}>
            <span class={`${classes.fade}`} style="flex: 1; min-width: 0">
              {label}
              <span class={classes.caption}> {derive(cut, (c) => (c?.kept ? `· ${c.kept.count} atoms` : ''))}</span>
            </span>
            <Button
              size="xs"
              variant="ghost"
              title="Cut these atoms again"
              onClick={() => props.run('clipping.release', {})}
            >
              Release
            </Button>
          </div>
        )}
      </Show>
      <Button
        size="xs"
        variant="ghost"
        block
        icon={<Icon name="plus" />}
        disabled={derive(selected, (s) => !s)}
        title={derive(selected, (s) =>
          s ? 'The selected atoms left whole by this cut' : 'Select the atoms to keep whole',
        )}
        onClick={() => props.run('clipping.keep', {})}
      >
        <span style="flex: 1; text-align: start">Keep the selection whole</span>
      </Button>
    </div>
  );
}

/** How a drawn cut looks in the scene (its walls, its edges): each change applied at once, the latest wins */
function DrawingForm(props: {
  clipping: ClippingApi;
  cut: ReadonlyState<Cut | undefined>;
  run(op: string, values: Record<string, unknown>): Promise<unknown> | void;
}) {
  const { clipping, cut } = props;
  let running = false;
  let latest: Record<string, unknown> | undefined;
  const flush = async () => {
    running = true;
    while (latest) {
      const options = latest;
      latest = undefined;
      await props.run('clipping.show', { on: true, options });
    }
    running = false;
    // as the drawing is now: an edit refused (or nothing changed) leaves no change to follow
    const c = cut.value;
    if (c) form.load(clipping.drawingValues(c) as never);
  };
  const form = new ParamsFormModel(clipping.drawingParams, cut.value ? clipping.drawingValues(cut.value) : undefined, {
    onChange: (values) => {
      latest = values;
      if (!running) void flush();
    },
  });
  effect(cut, (c) => {
    if (c && !running) form.load(clipping.drawingValues(c) as never);
  });
  return <ParamsForm form={form} />;
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
  // Only the changed fields are sent: the others as shown are rounded (a box's up turns it from where it points, its
  // twist kept: clipping.update).
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
    // as the cut is now: an edit refused (a zero direction) leaves no change to follow
    const c = cut.value;
    if (c) form.load(formValues(c) as never);
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
