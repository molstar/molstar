/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Button, ParamsForm, Popover } from '@molstar/prism-ui';
import { Icon } from '@molstar/prism-ui/components/icon/icon';
import { For, Show } from '@molstar/prism-ui/dom/control-flow';
import { effect } from '@molstar/prism-ui/dom/scope';
import { setAttr } from '@molstar/prism-ui/dom/writes';
import { observableState } from '@molstar/prism-ui/reactive/observable-state';
import type { ReadonlyState } from '@molstar/prism-ui/reactive/readonly-state';
import { derive } from '@molstar/prism-ui/reactive/state';
import type { LookOverrideSummary, NodeLookModel, SceneLookModel } from '../../model/look.js';
import { classes } from './_generated/look-popover.module.css.js';

// The look's popovers (api.md §2, "Look"): the scene defaults from the Structure panel's sliders, a node's own values
// from its inspector's sliders. Both are the same form; under it, what keeps values of its own, each with the way back.

/** The names in a tooltip: a few, then how many more. */
const namesTip = (names: readonly string[]) =>
  names.length <= 6 ? names.join('\n') : `${names.slice(0, 6).join('\n')}\n…and ${names.length - 6} more`;

/**
 * The sliders button: opens a look popover; a dot at its corner while `marked` (overrides are in place), said in its
 * name too.
 */
export function LookButton(props: {
  title: string;
  /** Overrides are in place: what the dot says. */
  marked: ReadonlyState<string | undefined>;
  class?: string;
  popover: (anchor: HTMLElement, close: (refocus?: boolean) => void) => unknown;
}) {
  const open = observableState(false);
  let anchor: HTMLElement | undefined;
  // its dialog's state, said to assistive technology (the press state is the button's own)
  effect(open.asReadonly, (o) => anchor && setAttr(anchor, 'aria-expanded', o ? 'true' : 'false'));
  return (
    <span
      class={classes.button}
      ref={(el: HTMLElement) => {
        anchor = el.querySelector('button') ?? undefined;
        if (anchor) setAttr(anchor, 'aria-haspopup', 'dialog');
        if (anchor) setAttr(anchor, 'aria-expanded', open.value ? 'true' : 'false');
      }}
    >
      <Button
        size="sm"
        variant="ghost"
        class={props.class}
        icon={<Icon name="sliders-horizontal" />}
        title={derive(props.marked, (m) => (m ? `${props.title} (${m})` : props.title))}
        pressed={open.asReadonly}
        onClick={() => open.set(!open.value)}
      />
      <Show when={props.marked}>{() => <span class={classes.dot} aria-hidden="true" />}</Show>
      <Show when={open.asReadonly}>
        {() =>
          props.popover(anchor!, (refocus) => {
            open.set(false);
            if (refocus) anchor?.focus();
          }) as never
        }
      </Show>
    </span>
  );
}

/** "2 representations keep their own" for the button's name and dot. */
export function sceneOverridesText(list: readonly LookOverrideSummary[]): string | undefined {
  const n = new Set(list.flatMap((o) => o.names)).size;
  if (!n) return undefined;
  return n === 1 ? 'a representation has its own look' : `${n} representations have their own look`;
}

/** The scene defaults: the form, then the values some representations keep as their own, each with Apply to all. */
export function SceneLookPopover(props: {
  model: SceneLookModel;
  anchor: HTMLElement;
  onClose(refocus?: boolean): void;
}) {
  const { model } = props;
  return (
    <Popover anchor={props.anchor} onClose={props.onClose} title="Scene defaults">
      <ParamsForm form={model.form} />
      <Show when={derive(model.overrides, (o) => o.length > 0)}>
        {() => (
          <section class={classes.overrides} aria-label="Kept as their own">
            <span class={classes.overridesTitle}>
              <Icon name="sliders-horizontal" size={12} class={classes.mark} />
              Kept as their own
            </span>
            <For each={model.overrides} by={(o) => o.key}>
              {(o) => (
                <div class={classes.override} data-tip={namesTip(o.names)}>
                  <span class={classes.overrideLabel}>
                    {o.label}{' '}
                    <span class={classes.overrideCount}>
                      {o.names.length === 1 ? '· 1 representation' : `· ${o.names.length} representations`}
                    </span>
                  </span>
                  <Button
                    size="xs"
                    variant="ghost"
                    title={`Use the default ${o.label.toLowerCase()} everywhere`}
                    onClick={() => model.applyToAll(o.key)}
                  >
                    Apply to all
                  </Button>
                </div>
              )}
            </For>
          </section>
        )}
      </Show>
    </Popover>
  );
}

/** A node's look: the form (its representations' values), then its own values, each with Use default. */
export function NodeLookPopover(props: {
  model: NodeLookModel;
  title: ReadonlyState<string>;
  anchor: HTMLElement;
  onClose(refocus?: boolean): void;
}) {
  const { model } = props;
  return (
    <Popover anchor={props.anchor} onClose={props.onClose} title={props.title}>
      <ParamsForm form={model.form} />
      <Show when={derive(model.own, (o) => o.length > 0)}>
        {() => (
          <section class={classes.overrides} aria-label="Its own values">
            <span class={classes.overridesTitle}>
              <Icon name="sliders-horizontal" size={12} class={classes.mark} />
              Not the scene defaults
            </span>
            <For each={model.own} by={(o) => o.key}>
              {(o) => (
                <div class={classes.override} data-tip={namesTip(o.names)}>
                  <span class={classes.overrideLabel}>
                    {o.label}
                    {o.names.length > 1 ? (
                      <span class={classes.overrideCount}>{` · ${o.names.length} representations`}</span>
                    ) : null}
                  </span>
                  <Button
                    size="xs"
                    variant="ghost"
                    title={`Use the scene's ${o.label.toLowerCase()}`}
                    onClick={() => model.reset(o.key)}
                  >
                    Use default
                  </Button>
                </div>
              )}
            </For>
            <Show when={derive(model.own, (o) => o.length > 1)}>
              {() => (
                <Button
                  size="xs"
                  variant="ghost"
                  icon={<Icon name="rotate-ccw" size={12} />}
                  onClick={() => model.reset()}
                >
                  Use all the scene defaults
                </Button>
              )}
            </Show>
          </section>
        )}
      </Show>
    </Popover>
  );
}
