/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { map, of } from 'rxjs';
import { Icon, type IconName } from '@molstar/prism-ui/components/icon/icon';
import type { CanvasMenuItem } from '@molstar/prism-ui';
import { Show } from '@molstar/prism-ui/dom/control-flow';
import { onCleanup } from '@molstar/prism-ui/dom/scope';
import { deepActiveElement, focusFirstItem, menuKeys } from '@molstar/prism-ui/dom/focus';
import type { ReadonlyState } from '@molstar/prism-ui/reactive/readonly-state';
import { classes } from './_generated/mode-switch.module.css.js';

import { derive } from '@molstar/prism-ui/reactive/state';
/** The picking level's menu (one per shadow root or page: the viewer's shell is in one): Select controls it. */
const LevelMenuId = 'msp-pick-level';

/** The levels' menu next to `el` (in its shadow root or document). */
function levelMenu(el: Element): HTMLElement | null {
  return (el.getRootNode() as Document | ShadowRoot).getElementById?.(LevelMenuId) ?? null;
}

/** What a click on the canvas does: interact (focus, rotate, measure) or select. */
export type CanvasMode = 'interact' | 'select';

/** A picking level (Mol*'s granularity): its value, label, icon, and an optional group (multi-instance levels). */
export interface Granularity {
  value: string;
  label: string;
  icon: IconName;
  group?: string;
}

export interface ModeSwitchProps {
  mode: ReadonlyState<CanvasMode>;
  setMode(mode: CanvasMode): void;
  /** What a click selects (residue, chain, …), offered below the switch while selecting. */
  granularity?: ReadonlyState<string>;
  granularities?: readonly Granularity[];
  setGranularity?(value: string): void;
}

/**
 * The canvas mode, at the top of the canvas stack: two toggle buttons, Interact and Select (the titles explain the
 * difference; Select takes the accent and carries the picking level as a small badge). `PickLevel` changes the level:
 * it shows beside Select while the pointer is on Select or on it.
 */
export function ModeSwitch(props: ModeSwitchProps) {
  const selecting = derive(props.mode, (m) => m === 'select');
  // The picking level as a badge on Select: what a click picks, at a glance.
  const levels = props.granularities;
  const current =
    levels && props.granularity
      ? derive(props.granularity, (g) => levels.find((l) => l.value === g) ?? levels[0])
      : undefined;
  return (
    <div class={classes.switch} role="group" aria-label="Canvas mode" data-canvas-mode={props.mode}>
      <button
        type="button"
        class={classes.segment}
        aria-pressed={selecting.pipe(map((s) => !s))}
        data-tip="Interact: drag to rotate, click to focus"
        data-tip-side="top"
        aria-label="Interact"
        onClick={() => props.setMode('interact')}
      >
        <Icon name="rotate-3d" size={20} />
      </button>
      <button
        type="button"
        class={`${classes.segment} ${classes.select}`}
        aria-pressed={selecting}
        data-tip={current ? undefined : 'Select: click to select, Shift+click to extend'}
        data-tip-side="top"
        aria-label="Select"
        aria-haspopup={current ? 'menu' : undefined}
        aria-controls={current ? LevelMenuId : undefined}
        aria-expanded={current ? false : undefined}
        onClick={() => props.setMode('select')}
        onKeyDown={
          current
            ? (e) => {
                // ↑ or ↓ goes into the picking levels (the menu button pattern), at the current one.
                if (e.key !== 'ArrowUp' && e.key !== 'ArrowDown') return;
                const menu = levelMenu(e.currentTarget);
                if (menu && focusFirstItem(menu, true)) e.preventDefault();
              }
            : undefined
        }
      >
        <Icon name="mouse-pointer-click" size={20} />
        {current ? <Show when={current}>{(l) => <Icon name={l.icon} size={10} class={classes.badge} />}</Show> : null}
      </button>
    </div>
  );
}

/**
 * The mode switch and the picking level as rows of the collapsed toolbar's menu (ui shell.md §2): Interact and Select
 * as a choice, and the levels as a submenu (picking one turns select mode on, as in `PickLevel`).
 */
export function modeMenuRows(props: ModeSwitchProps): readonly CanvasMenuItem[] {
  const selecting = derive(props.mode, (m) => m === 'select');
  const rows: CanvasMenuItem[] = [
    {
      id: 'mode-interact',
      label: 'Interact',
      icon: 'rotate-3d',
      group: 'mode',
      choice: true,
      pressed: selecting.pipe(map((s) => !s)),
      run: () => props.setMode('interact'),
    },
    {
      id: 'mode-select',
      label: 'Select',
      icon: 'mouse-pointer-click',
      group: 'mode',
      choice: true,
      pressed: selecting,
      run: () => props.setMode('select'),
    },
  ];
  const { granularity, granularities, setGranularity } = props;
  if (granularity && granularities && setGranularity) {
    rows.push({
      id: 'mode-level',
      label: 'What a click selects',
      group: 'mode',
      run: () => {},
      items: of(
        granularities.map(
          (l): CanvasMenuItem => ({
            id: `level-${l.value}`,
            label: l.label,
            icon: l.icon,
            group: l.group ?? '',
            choice: true,
            pressed: derive(granularity, (g) => g === l.value),
            run: () => {
              setGranularity(l.value);
              props.setMode('select');
            },
          }),
        ),
      ),
    });
  }
  return rows;
}

/**
 * The picking levels, above the canvas toolbar over Select (outside the toolbar's frame, so it is not clipped): the list
 * itself, one click from the toolbar: the current level marked; picking one sets it and turns select mode on. Shown while
 * the pointer is on Select or on it, or Select or a row has the keyboard focus (a click's focus does not keep it); on touch screens (no hover) a tap on
 * Select opens it. It opens upward. From the keyboard: ↑ or ↓ on Select goes into it (at the current level), the
 * arrows move between levels, Escape or a pick goes back to Select (it is outside the toolbar in the DOM, so Tab does
 * not reach it).
 */
export function PickLevel(props: ModeSwitchProps) {
  const levels = props.granularities;
  const granularity = props.granularity;
  if (!levels || !granularity) return null;
  // On a touch screen it shows while Select has the focus, and a tap on the canvas (not focusable) does not take the
  // focus away: a tap outside the list and Select lets it go, so the list closes.
  let level: HTMLElement | undefined;
  const onDown = (e: PointerEvent) => {
    const active = deepActiveElement();
    if (!level || !(active instanceof HTMLElement)) return;
    const menu = level.querySelector(`#${LevelMenuId}`);
    const select = menu && opener(menu);
    const holds = (el: Element | null | undefined) => !!el && (el === active || el.contains(active));
    if (!holds(menu) && !holds(select)) return;
    const path = e.composedPath();
    if ((menu && path.includes(menu)) || (select && path.includes(select))) return;
    active.blur();
  };
  document.addEventListener('pointerdown', onDown, true);
  onCleanup(() => document.removeEventListener('pointerdown', onDown, true));
  return (
    <div class={classes.level} ref={(el) => (level = el)}>
      <div
        class={classes.levels}
        id={LevelMenuId}
        role="menu"
        aria-label="What a click selects"
        onKeyDown={(e) => {
          menuKeys(e);
          if (e.key !== 'Escape' && e.key !== 'Tab') return;
          if (e.key === 'Escape') {
            e.preventDefault();
            e.stopPropagation();
          }
          if (e.key === 'Tab') e.preventDefault();
          opener(e.currentTarget)?.focus();
        }}
        onFocusIn={(e) => opener(e.currentTarget)?.setAttribute('aria-expanded', 'true')}
        onFocusOut={(e) => {
          const menu = e.currentTarget;
          if (menu.contains(e.relatedTarget as Node | null)) return;
          opener(menu)?.setAttribute('aria-expanded', 'false');
        }}
      >
        {levels.map((l, i) => (
          <>
            {l.group && l.group !== levels[i - 1]?.group ? <div class={classes.group}>{l.group}</div> : null}
            <button
              type="button"
              role="menuitemradio"
              class={classes.levelRow}
              aria-checked={derive(granularity, (g) => g === l.value)}
              onClick={(e) => {
                const menu = e.currentTarget.closest<HTMLElement>(`#${LevelMenuId}`);
                const inside = !!menu && menu.contains(deepActiveElement() ?? null);
                props.setGranularity?.(l.value);
                props.setMode('select');
                // Picked from the keyboard: back to Select (the list hides once the focus leaves it). Picked with the
                // pointer: the row lets the focus go, so the list closes as the pointer leaves it.
                if (inside && e.detail === 0) opener(menu)?.focus();
                else if (e.detail > 0) e.currentTarget.blur();
              }}
            >
              <Icon name={l.icon} size={15} />
              {l.label}
            </button>
          </>
        ))}
      </div>
    </div>
  );
}

/** The Select button that controls the levels' menu. */
function opener(menu: Element): HTMLElement | null {
  return (menu.getRootNode() as Document | ShadowRoot).querySelector?.(`[aria-controls="${LevelMenuId}"]`) ?? null;
}
