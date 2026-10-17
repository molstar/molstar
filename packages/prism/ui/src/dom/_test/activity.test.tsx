/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { BehaviorSubject } from 'rxjs';
import { describe, expect, it, vi } from 'vitest';
import { observableState } from '../../reactive/observable-state.js';
import { ViewHost } from '../../views/shell/view-host.js';
import { activity, For, Match, Show, Switch } from '../control-flow.js';
import type { Child } from '../jsx-types.js';
import { render } from '../render.js';
import { bind, onActive, onCleanup, runInScope, Scope } from '../scope.js';
import { domView, type ViewRenderer } from '../view.js';

import { derive } from '../../reactive/state.js';
function mount(fn: () => Child) {
  const container = document.createElement('div');
  const dispose = render(fn, container);
  return { container, dispose, html: () => container.innerHTML.replace(/<!--[^>]*-->/g, '') };
}

describe('Scope.setActive', () => {
  it('drops the binding subscription while paused and subscribes again, with the current value, on resume', () => {
    const source = new BehaviorSubject('a');
    const seen: string[] = [];
    const scope = new Scope();
    runInScope(scope, () => bind(source, (v) => seen.push(v)));
    expect(seen).toEqual(['a']);
    expect(source.observed).toBe(true);

    scope.setActive(false);
    expect(source.observed).toBe(false);
    source.next('b');
    source.next('c');
    expect(seen).toEqual(['a']);

    scope.setActive(true);
    expect(seen).toEqual(['a', 'c']);
    expect(source.observed).toBe(true);
    source.next('d');
    expect(seen).toEqual(['a', 'c', 'd']);

    scope.dispose();
    expect(source.observed).toBe(false);
  });

  it('pauses the children with it, and a child made while paused starts paused', () => {
    const source = new BehaviorSubject(1);
    const parent = new Scope();
    const early = parent.child();
    const seen: string[] = [];
    runInScope(early, () => bind(source, (v) => seen.push(`early ${v}`)));
    parent.setActive(false);
    expect(source.observed).toBe(false);
    expect(early.active).toBe(false);

    const late = parent.child();
    expect(late.active).toBe(false);
    runInScope(late, () => bind(source, (v) => seen.push(`late ${v}`)));
    expect(source.observed).toBe(false);

    // A child cannot be more active than its parent.
    late.setActive(true);
    expect(late.active).toBe(false);

    source.next(2);
    parent.setActive(true);
    expect(seen.sort()).toEqual(['early 1', 'early 2', 'late 2']);
  });

  it('a child that was paused itself stays paused when its parent resumes', () => {
    const source = new BehaviorSubject(1);
    const parent = new Scope();
    const child = parent.child();
    runInScope(child, () => bind(source, () => {}));
    child.setActive(false);
    parent.setActive(false);
    parent.setActive(true);
    expect(child.active).toBe(false);
    expect(source.observed).toBe(false);
    child.setActive(true);
    expect(source.observed).toBe(true);
  });

  it('runs onActive hooks with the state they start in and on each change, not after dispose', () => {
    const scope = new Scope();
    const calls: boolean[] = [];
    runInScope(scope, () => onActive((a) => calls.push(a)));
    scope.setActive(false);
    scope.setActive(false);
    scope.setActive(true);
    expect(calls).toEqual([true, false, true]);
    scope.dispose();
    scope.setActive(false);
    expect(calls).toEqual([true, false, true]);

    const paused = new Scope();
    paused.setActive(false);
    const late: boolean[] = [];
    runInScope(paused, () => onActive((a) => late.push(a)));
    expect(late).toEqual([false]);
  });

  it('a throwing hook does not keep the others from running', () => {
    const scope = new Scope();
    const error = vi.spyOn(console, 'error').mockImplementation(() => {});
    const ran: boolean[] = [];
    scope.onActive(() => {
      throw new Error('x');
    });
    scope.onActive((a) => ran.push(a));
    scope.setActive(false);
    expect(ran).toEqual([true, false]);
    error.mockRestore();
  });

  it('follows a stream, and a stream that is set while its subscription is being made does not leave it subscribed', () => {
    const flag = observableState(false);
    const source = new BehaviorSubject(1);
    const scope = new Scope();
    scope.follow(flag);
    runInScope(scope, () => bind(source, () => {}));
    expect(source.observed).toBe(false);
    flag.set(true);
    expect(source.observed).toBe(true);
    flag.set(false);
    expect(source.observed).toBe(false);
  });
});

describe('activity', () => {
  it('a bound text does not change while hidden and shows the latest value on show; the part is inert meanwhile', () => {
    const shown = observableState(true);
    const text = observableState('a');
    const { container } = mount(() =>
      activity(shown, () => (
        <section>
          <b>{text}</b>
        </section>
      )),
    );
    const section = container.querySelector('section') as HTMLElement;
    const b = container.querySelector('b') as HTMLElement;
    expect(section.hasAttribute('inert')).toBe(false);
    text.set('b');
    expect(b.textContent).toBe('b');

    shown.set(false);
    expect(section.hasAttribute('inert')).toBe(true);
    const mutations: MutationRecord[] = [];
    const observer = new MutationObserver((m) => mutations.push(...m));
    observer.observe(section, { subtree: true, childList: true, characterData: true, attributes: true });
    text.set('c');
    text.set('d');
    expect(b.textContent).toBe('b');
    expect(observer.takeRecords()).toEqual([]);

    shown.set(true);
    expect(b.textContent).toBe('d');
    expect(section.hasAttribute('inert')).toBe(false);
    text.set('e');
    expect(b.textContent).toBe('e');
    observer.disconnect();
  });

  it('keeps the DOM mounted: the same nodes, the scopes (and what lives in them) are not rebuilt', () => {
    const shown = observableState(true);
    const text = observableState('a');
    let renders = 0;
    let cleaned = 0;
    const { container } = mount(() =>
      activity(shown, () => (
        <div>
          {() => {
            renders++;
            onCleanup(() => cleaned++);
            return <i>{text}</i>;
          }}
        </div>
      )),
    );
    const i = container.querySelector('i');
    shown.set(false);
    text.set('b');
    shown.set(true);
    expect(container.querySelector('i')).toBe(i);
    expect(i?.textContent).toBe('b');
    expect([renders, cleaned]).toEqual([1, 0]);
  });

  it('a part that starts hidden renders its structure and binds nothing until it is shown', () => {
    const shown = observableState(false);
    const source = new BehaviorSubject('x');
    const { container } = mount(() =>
      activity(shown, () => (
        <p>
          <span>{source}</span>
        </p>
      )),
    );
    expect(container.querySelector('p')?.hasAttribute('inert')).toBe(true);
    expect(container.querySelector('span')?.textContent).toBe('');
    expect(source.observed).toBe(false);
    shown.set(true);
    expect(container.querySelector('span')?.textContent).toBe('x');
    expect(container.querySelector('p')?.hasAttribute('inert')).toBe(false);
  });

  it('a nested For catches up on show: kept rows keep their nodes, rows follow their own bindings', () => {
    const shown = observableState(true);
    const names = observableState({ a: 'A', b: 'B', c: 'C' });
    const keys = observableState(['a', 'b']);
    const { container } = mount(() =>
      activity(shown, () => (
        <ul>
          <For each={keys}>{(k) => <li>{derive(names.asReadonly, (n) => n[k as 'a'])}</li>}</For>
        </ul>
      )),
    );
    const text = () => [...container.querySelectorAll('li')].map((li) => li.textContent);
    const first = container.querySelector('li');
    expect(text()).toEqual(['A', 'B']);

    shown.set(false);
    keys.set(['a', 'b', 'c']);
    names.set({ a: 'A!', b: 'B', c: 'C' });
    expect(text()).toEqual(['A', 'B']);

    shown.set(true);
    expect(text()).toEqual(['A!', 'B', 'C']);
    expect(container.querySelector('li')).toBe(first);

    shown.set(false);
    keys.set(['c']);
    shown.set(true);
    expect(text()).toEqual(['C']);
  });

  it('a Show inside keeps its content across a pause (not rendered again) and follows what changed meanwhile', () => {
    const shown = observableState(true);
    const on = observableState(true);
    const item = observableState({ n: 1 });
    let renders = 0;
    const { container } = mount(() =>
      activity(shown, () => (
        <div>
          <Show when={on}>
            {() => {
              renders++;
              return <b>on</b>;
            }}
          </Show>
          <Show when={item}>{(v) => <i>{v.n}</i>}</Show>
        </div>
      )),
    );
    const b = container.querySelector('b');
    shown.set(false);
    shown.set(true);
    expect(renders).toBe(1);
    expect(container.querySelector('b')).toBe(b);

    shown.set(false);
    on.set(false);
    item.set({ n: 2 });
    expect(container.querySelector('b')).toBe(b);
    expect(container.querySelector('i')?.textContent).toBe('1');
    shown.set(true);
    expect(container.querySelector('b')).toBeNull();
    expect(container.querySelector('i')?.textContent).toBe('2');
  });

  it('a Show over a value derived anew on every subscription keeps its content across a pause', () => {
    const shown = observableState(true);
    const focus = observableState<string | undefined>(undefined);
    let renders = 0;
    const { container } = mount(() =>
      activity(shown, () => (
        <div>
          <Show when={derive(focus.asReadonly, (section) => ({ section }))}>
            {() => {
              renders++;
              return <b />;
            }}
          </Show>
        </div>
      )),
    );
    const b = container.querySelector('b');
    shown.set(false);
    shown.set(true);
    expect(renders).toBe(1);
    expect(container.querySelector('b')).toBe(b);
    shown.set(false);
    focus.set('measure');
    shown.set(true);
    expect(renders).toBe(2);
  });

  it('a Switch inside keeps its case across a pause', () => {
    const shown = observableState(true);
    const mode = observableState('a');
    let renders = 0;
    const { container } = mount(() =>
      activity(shown, () => (
        <div>
          <Switch>
            <Match when={derive(mode.asReadonly, (m) => m === 'a')}>
              {() => {
                renders++;
                return <b>a</b>;
              }}
            </Match>
            <Match when={derive(mode.asReadonly, (m) => m === 'b')}>{() => <i>b</i>}</Match>
          </Switch>
        </div>
      )),
    );
    shown.set(false);
    shown.set(true);
    expect(renders).toBe(1);
    shown.set(false);
    mode.set('b');
    expect(container.querySelector('i')).toBeNull();
    shown.set(true);
    expect(container.querySelector('i')?.textContent).toBe('b');
    expect(container.querySelector('b')).toBeNull();
  });

  it('is disposed with its owner while paused, leaving no subscription', () => {
    const shown = observableState(false);
    const source = new BehaviorSubject('x');
    const { dispose } = mount(() => activity(shown, () => <p>{source}</p>));
    shown.set(true);
    expect(source.observed).toBe(true);
    shown.set(false);
    dispose();
    shown.set(true);
    expect(source.observed).toBe(false);
  });

  it('nests: an outer pause holds the inner part, and the inner flag still counts after it', () => {
    const outer = observableState(true);
    const inner = observableState(true);
    const text = observableState('a');
    const { container } = mount(() =>
      activity(outer, () => (
        <div>
          {activity(inner, () => (
            <p>{text}</p>
          ))}
        </div>
      )),
    );
    const p = container.querySelector('p') as HTMLElement;
    outer.set(false);
    text.set('b');
    expect(p.textContent).toBe('a');
    inner.set(false);
    outer.set(true);
    expect(p.textContent).toBe('a');
    expect(p.hasAttribute('inert')).toBe(true);
    inner.set(true);
    expect(p.textContent).toBe('b');
  });
});

describe('views in an inactive scope', () => {
  it('a domView follows the scope it is mounted in', () => {
    const shown = observableState(true);
    const text = observableState('a');
    const view = domView<{ text: typeof text }>((props) => <em>{props.ctx.text}</em>);
    const { container } = mount(() => activity(shown, () => <ViewHost renderer={view} ctx={{ text }} />));
    const em = container.querySelector('em') as HTMLElement;
    shown.set(false);
    text.set('b');
    expect(em.textContent).toBe('a');
    shown.set(true);
    expect(em.textContent).toBe('b');
  });

  it('an imperative view pauses and resumes its own work through scope.onActive', () => {
    const shown = observableState(true);
    const events: string[] = [];
    const view: ViewRenderer<undefined> = {
      mount(host, _ctx, scope) {
        host.textContent = 'imperative';
        scope.onActive((active) => events.push(active ? 'start' : 'stop'));
        return () => events.push('dispose');
      },
    };
    const { dispose } = mount(() => activity(shown, () => <ViewHost renderer={view} ctx={undefined} />));
    shown.set(false);
    shown.set(true);
    shown.set(false);
    expect(events).toEqual(['start', 'stop', 'start', 'stop']);
    dispose();
    expect(events).toEqual(['start', 'stop', 'start', 'stop', 'dispose']);
  });
});

describe('Show closing', () => {
  it('makes the content inert while its closing animation plays, then removes it', async () => {
    let finish!: () => void;
    const finished = new Promise<void>((resolve) => {
      finish = resolve;
    });
    const animate = vi.fn(() => ({ finished, cancel: () => {} }));
    (HTMLElement.prototype as any).animate = animate;
    const on = observableState(true);
    const { container } = mount(() => (
      <Show when={on} reveal="side">
        {() => <section />}
      </Show>
    ));
    document.body.append(container);
    const section = container.querySelector('section') as HTMLElement;
    on.set(false);
    expect(animate).toHaveBeenCalled();
    expect(container.querySelector('section')).toBe(section);
    expect(section.hasAttribute('inert')).toBe(true);
    finish();
    await finished;
    await Promise.resolve();
    expect(container.querySelector('section')).toBeNull();
    container.remove();
    delete (HTMLElement.prototype as any).animate;
  });
});
