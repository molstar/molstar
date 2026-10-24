/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Choices } from '@molstar/prism-ui';
import { Icon } from '@molstar/prism-ui/components/icon/icon';
import { For, Show } from '@molstar/prism-ui/dom/control-flow';
import { observableState } from '@molstar/prism-ui/reactive/observable-state';
import { derive } from '@molstar/prism-ui/reactive/state';
import type { SequenceChainsEntity, SequenceChainsModel } from '../../model/scene.js';
import { classes } from './_generated/sequence-chains.module.css.js';

// The chains a structure's sequence shows, in its Sequence node's inspector (sequence.md §8.1): the presets as a row of
// choices (every chain, none, the polymers; the one in effect marked), then the entities, each one line: a box that
// shows or leaves out its chains (mixed while some are shown), a dot for its kind, its name and how many of its chains
// are shown. An entity of several chains opens to its chains as toggles that wrap. A mesoscale structure (hundreds of
// entities and chains): the first entities, then Show all; an entity's chains only once it is opened.

const PageSize = 8;

const Kinds: Readonly<Record<SequenceChainsEntity['kind'], { label: string; class: string }>> = {
  polymer: { label: 'Polymer', class: classes.polymer },
  ligand: { label: 'Ligand', class: classes.ligand },
  branched: { label: 'Carbohydrate', class: classes.branched },
};

export function SequenceChains(props: { model: SequenceChainsModel }) {
  const { model } = props;
  const all = observableState(model.entities.length <= PageSize + 1);
  const listed = derive(all.asReadonly, (a) => (a ? model.entities : model.entities.slice(0, PageSize)));
  return (
    <section class={classes.chains} aria-label="Chains">
      <Choices
        label="Chains"
        choices={[
          { id: 'all', label: 'All', title: 'Every chain' },
          { id: 'none', label: 'None', title: 'No chain: then show them one by one' },
          { id: 'polymer', label: 'Polymers', title: 'The polymer chains, no ligands' },
        ]}
        active={model.preset}
        onPick={(id) => model.apply(id as 'all' | 'none' | 'polymer')}
      />
      <div class={classes.list} role="list" aria-label="Entities">
        <For each={listed} by={(e) => e.id}>
          {(e) => <EntityRow model={model} entity={e} />}
        </For>
      </div>
      <Show when={derive(all.asReadonly, (a) => !a)}>
        {() => (
          <button type="button" class={classes.more} onClick={() => all.set(true)}>
            {`Show all ${model.entities.length}`}
            <span class={classes.moreCount}>{`${model.entities.length - PageSize} more`}</span>
          </button>
        )}
      </Show>
    </section>
  );
}

function EntityRow(props: { model: SequenceChainsModel; entity: SequenceChainsEntity }) {
  const { model, entity } = props;
  const kind = Kinds[entity.kind];
  const several = entity.chains.length > 1;
  const open = observableState(false);
  const count = derive(model.shown, (s) => entity.chains.filter((c) => s.has(c.name)).length);
  const checked = derive(count, (n) => (n === 0 ? 'false' : n === entity.chains.length ? 'true' : 'mixed'));
  const meta = derive(count, (n) =>
    several ? `${n} of ${entity.chains.length} chains` : `${kind.label} · chain ${entity.chains[0].label}`,
  );
  return (
    <div role="listitem" class={classes.entity}>
      <div class={classes.row}>
        <button
          type="button"
          role="checkbox"
          class={classes.box}
          aria-checked={checked}
          aria-label={`Show ${entity.name}`}
          onClick={() => model.toggleEntity(entity.id)}
        >
          <Show when={derive(checked, (c) => c !== 'false')}>
            {() => <Icon name={checked.value === 'mixed' ? 'minus' : 'check'} size={12} />}
          </Show>
        </button>
        <span class={`${classes.dot} ${kind.class}`} data-tip={kind.label} aria-label={kind.label} />
        <button
          type="button"
          class={classes.text}
          aria-expanded={several ? open.asReadonly : undefined}
          onClick={() => (several ? open.set(!open.value) : model.toggleEntity(entity.id))}
        >
          <span class={classes.name} data-tip-overflow="">
            {entity.name}
          </span>
          <span class={classes.meta}>{meta}</span>
        </button>
        {several ? (
          <span class={classes.chevron} classList={{ [classes.chevronOpen]: open.asReadonly }}>
            <Icon name="chevron-right" size={14} />
          </span>
        ) : null}
      </div>
      <Show when={open.asReadonly}>
        {() => (
          <div class={classes.chips} role="group" aria-label={`Chains of ${entity.name}`}>
            {entity.chains.map((c) => (
              <button
                type="button"
                class={classes.chip}
                aria-pressed={derive(model.shown, (s) => s.has(c.name))}
                title={c.label}
                onClick={() => model.toggleChain(c.name)}
              >
                {c.label}
              </button>
            ))}
          </div>
        )}
      </Show>
    </div>
  );
}
