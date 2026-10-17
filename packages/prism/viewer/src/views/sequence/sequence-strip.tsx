/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type Observable } from 'rxjs';
import { Button } from '@molstar/prism-ui';
import { Icon } from '@molstar/prism-ui/components/icon/icon';
import { ScrollArea } from '@molstar/prism-ui';
import { OptionList, type SelectOption } from '@molstar/prism-ui';
import { Show } from '@molstar/prism-ui/dom/control-flow';
import { observableState } from '@molstar/prism-ui/reactive/observable-state';
import type { ReadonlyState } from '@molstar/prism-ui/reactive/readonly-state';
import { classes } from './_generated/sequence-strip.module.css.js';

import { derive } from '@molstar/prism-ui/reactive/state';
/**
 * One of the sequence view's choices (Mol*'s sequence controls): the structure, the mode (Chain, Polymers,
 * Everything), the entity and the chain (in Chain mode), the instance (an assembly's operators, when there are several).
 */
export interface SequenceSelector {
  id: string;
  /** What it chooses ("Structure", "Mode"): its tooltip. */
  title: string;
  options: readonly SelectOption[];
  value: ReadonlyState<string>;
  set(value: string): void;
  /** The button's text for a value, shorter than the option ("1TQN" for "1TQN | Crystal structure of …"). */
  display?: (value: string) => string;
  /** Shown only while true (entity and chain in Chain mode; the instance with several). */
  visible?: Observable<boolean>;
  /** Set in the text color and bold (the structure); the others quieter. */
  strong?: boolean;
}

export interface SequenceChain {
  /** The structure it belongs to ("1TQN"). */
  structure: string;
  /** The chain ("A") and its entity description. */
  name: string;
  info?: string;
  /** One-letter codes. */
  sequence: string;
  /** The number of the first residue. */
  start: number;
}

/**
 * A sequence strip for the top dock (a mock of the sequence view). Residues are rendered once; hover and selection
 * are applied to the existing nodes, so thousands of residues cost no subscriptions or re-rendering.
 */
export function SequenceStrip(props: {
  chain: SequenceChain;
  onHover?: (seqNumber: number | undefined) => void;
  onToggle?: (seqNumber: number, selected: boolean) => void;
  /** Collapses the strip (to the dock's slim bar). */
  onCollapse?: () => void;
  /** The view's choices; without them, the chain's structure and name as plain selectors. */
  selectors?: readonly SequenceSelector[];
}) {
  const { chain } = props;
  const fragment = document.createDocumentFragment();
  for (let i = 0; i < chain.sequence.length; i++) {
    const span = document.createElement('span');
    span.className = classes.residue;
    span.textContent = chain.sequence[i];
    const number = chain.start + i;
    span.dataset.seq = String(number);
    if (number % 10 === 0) span.dataset.number = String(number);
    fragment.appendChild(span);
  }
  // One listener per strip (delegation), not per residue.
  const residues = (
    <ScrollArea
      class={classes.residues}
      contentClass={classes.residuesContent}
      onPointerOver={(e) => {
        const seq = (e.target as HTMLElement).dataset?.seq;
        props.onHover?.(seq === undefined ? undefined : Number(seq));
      }}
      onPointerLeave={() => props.onHover?.(undefined)}
      onClick={(e) => {
        const target = e.target as HTMLElement;
        if (target.dataset?.seq === undefined) return;
        const selected = target.classList.toggle(classes.selected);
        props.onToggle?.(Number(target.dataset.seq), selected);
      }}
    >
      {fragment}
    </ScrollArea>
  );
  // The current structure and chain as selectors on top (as in the classic UI), so the residues get the full width.
  return (
    <div class={classes.strip}>
      <div class={classes.header}>
        {props.selectors ? (
          props.selectors.map((selector) =>
            selector.visible ? (
              <Show when={selector.visible}>{() => <HeaderSelect selector={selector} />}</Show>
            ) : (
              <HeaderSelect selector={selector} />
            ),
          )
        ) : (
          <>
            <button type="button" class={classes.selector} title="Structure">
              {chain.structure} <Icon name="chevron-down" size={12} />
            </button>
            <button type="button" class={classes.selector} title="Chain">
              {chain.name}
              {chain.info ? <span class={classes.info}>{chain.info}</span> : null}
              <Icon name="chevron-down" size={12} />
            </button>
          </>
        )}
        <span class={classes.spacer} />
        {props.onCollapse ? (
          <Button
            variant="ghost"
            size="sm"
            icon={<Icon name="chevron-up" />}
            title="Hide the sequence"
            onClick={() => props.onCollapse?.()}
          />
        ) : null}
      </div>
      {residues}
    </div>
  );
}

/** A choice in the header: its value as quiet text and a chevron; the options open in a small list under it. */
function HeaderSelect(props: { selector: SequenceSelector }) {
  const { selector } = props;
  const open = observableState(false);
  // Fixed to the window under the button: the strip's card clips what overflows it.
  const at = observableState({ top: 0, left: 0 });
  const text = derive(
    selector.value,
    (v) => selector.display?.(v) ?? selector.options.find((o) => o[0] === v)?.[1] ?? v,
  );
  const tip = derive(selector.value, (v) => `${selector.title}: ${selector.options.find((o) => o[0] === v)?.[1] ?? v}`);
  return (
    <span class={classes.choice}>
      <button
        type="button"
        class={classes.selector}
        data-strong={selector.strong ? '' : undefined}
        title={tip}
        aria-label={selector.title}
        aria-haspopup="listbox"
        aria-expanded={open}
        onClick={(e) => {
          const r = (e.currentTarget as HTMLElement).getBoundingClientRect();
          at.set({ top: r.bottom + 2, left: r.left });
          open.set(!open.value);
        }}
      >
        <span class={classes.value}>{text}</span>
        <Icon name="chevron-down" size={12} class={classes.chevron} />
      </button>
      <Show when={open}>
        {() => (
          <div
            class={classes.popover}
            style={{
              top: derive(at.asReadonly, (p) => `${p.top}px`),
              left: derive(at.asReadonly, (p) => `${p.left}px`),
            }}
          >
            <OptionList
              options={selector.options}
              value={selector.value}
              onPick={(v) => {
                open.set(false);
                selector.set(v);
              }}
              onClose={() => open.set(false)}
              aria-label={selector.title}
            />
          </div>
        )}
      </Show>
    </span>
  );
}
