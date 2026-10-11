/** @jsxImportSource @molstar/plugin-ui-next/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ScrollArea } from '../../components/scroll-area/scroll-area.js';
import type { AgentModel } from '../../model/agent.js';
import { AISetup } from './ai-setup.js';
import { classes } from './_generated/assistant-tab.module.css.js';

/**
 * The Assistant tab, at the rail's foot by Settings (only with the agent module included, design/agent.md §6): the
 * model it uses and the keys; asking happens in the command menu.
 */
export function AssistantTab(props: { agent: AgentModel }) {
  return (
    <ScrollArea class={classes.tab} contentClass={classes.content}>
      <AISetup settings={props.agent.settings} />
      <p class={classes.hint}>
        Ask from the command menu (⌘K): type what you want, then <strong>Ask the assistant</strong>. Its work shows over
        the canvas and in History, and undoes in one step.
      </p>
    </ScrollArea>
  );
}
