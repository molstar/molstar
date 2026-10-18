/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { formatChord, type ReadonlyState, ScrollArea } from '@molstar/prism-ui';
import { watch } from '@molstar/prism-ui/dom/scope';
import { setText } from '@molstar/prism-ui/dom/writes';
import type { AgentModel } from '../model.js';
import { AISetup } from './ai-setup.js';
import { classes } from './_generated/assistant-tab.module.css.js';

/** A way into the Inspector from the Assistant tab: a view and what it holds ("3 turns", "14 tools, 3 skills"). */
export interface InspectLink {
  id: string;
  label: string;
  detail: string | ReadonlyState<string>;
}

/** The Inspect block: the Inspector's views, each opening the page there (the Inspector itself is a separate entry). */
export interface AssistantInspect {
  links: readonly InspectLink[];
  open(id: string): void;
}

/** One way into the Inspector: its name, and what it holds, which follows the Inspector's data. */
function InspectRow(props: { link: InspectLink; onOpen(): void }) {
  const { detail } = props.link;
  const text = document.createTextNode('');
  if (typeof detail === 'string') text.data = detail;
  else watch([detail], () => setText(text, detail.value));
  return (
    <button type="button" class={classes.link} onClick={props.onOpen}>
      <span>{props.link.label}</span>
      <span class={classes.detail}>{text}</span>
    </button>
  );
}

/**
 * The Assistant tab, at the rail's foot by Settings (only with the agent module included; _spec/README.md, "The
 * Assistant tab"): the model it uses and the keys, then (with the Inspector) a way into its trace, tools and limits;
 * asking happens in the command menu.
 */
export function AssistantTab(props: { agent: AgentModel; inspect?: AssistantInspect }) {
  const { inspect } = props;
  return (
    <ScrollArea class={classes.tab} contentClass={classes.content}>
      <AISetup settings={props.agent.settings} />
      {inspect ? (
        <section class={classes.inspect} aria-label="Inspect">
          <p class={classes.caption}>Inspect</p>
          {inspect.links.map((l) => (
            <InspectRow link={l} onOpen={() => inspect.open(l.id)} />
          ))}
        </section>
      ) : null}
      <p class={classes.hint}>
        Ask from the command menu ({formatChord({ key: 'k', mod: true })}): type what you want, then{' '}
        <strong>Ask AI assistant</strong>. Its work shows over the canvas and in History, and undoes in one step.
      </p>
    </ScrollArea>
  );
}
