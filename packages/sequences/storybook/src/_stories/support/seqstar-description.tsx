/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The DOM UI story's case 10 description: what a host does for a `markdown` description, which
// the viewer's default `Description` does not show. The host's part goes in through the `Description` slot of the `<seq-viewer>`
// element; it reads the description from the context (`descriptionOf`), renders it with its own markdown library
// (unified: remark-parse, remark-rehype, rehype-sanitize), marks the `!` links with `@molstar/seq-core/links`'s
// `rehypeSeqLinks` and runs them with its runner over `ctx.commands.run`. The link component is the reference for hosts:
//   - pointer enter and leave, focus and blur, and click call `runner.run(href, event, id)`;
//   - a click prevents the default (the href is relative and would navigate);
//   - the id is released when it leaves the tree or changes (an effect's cleanup keyed on the id);
//   - the runner is disposed with the viewer;
//   - each rendered document takes its own `prefix`, so the same link in two descriptions lights two channels.

import { createLinkRunner, type LinkRunner } from '@molstar/seq-core/links';
import { descriptionOf } from '@molstar/seq-viewer';
import type { DomPart, PartProps } from '@molstar/seq-viewer/dom';
import { type CSSProperties, createElement, type ReactNode, useEffect, useState } from 'react';
import { createRoot } from 'react-dom/client';
import { type HastNode, markdownToHast } from '../../../../core/_test/links/markdown.js';

const LINK: CSSProperties = { color: 'var(--seq-accent)', cursor: 'pointer', textDecoration: 'underline dotted' };

/** A `!` link: the runner's events, released when the link leaves the tree. */
function SeqLink({
  href,
  id,
  runner,
  children,
}: {
  href: string;
  id: string;
  runner: LinkRunner;
  children: ReactNode;
}) {
  useEffect(() => () => runner.release(id), [id, runner]);
  const run = (event: 'enter' | 'leave') => () => void runner.run(href, event, id);
  return (
    <a
      href={href}
      style={LINK}
      onPointerEnter={run('enter')}
      onPointerLeave={run('leave')}
      onFocus={run('enter')}
      onBlur={run('leave')}
      onClick={(e) => {
        e.preventDefault();
        runner.run(href, 'click', id);
      }}
    >
      {children}
    </a>
  );
}

const VOID = new Set(['br', 'hr', 'img']);
function toReact(n: HastNode, key: number, runner: LinkRunner): ReactNode {
  if (n.type === 'text') return n.value;
  if (n.type === 'root') return n.children?.map((c, i) => toReact(c, i, runner));
  if (n.type !== 'element' || !n.tagName) return null;
  const kids = (n.children ?? []).map((c, i) => toReact(c, i, runner));
  const props = n.properties ?? {};
  if (n.tagName === 'a' && typeof props.dataSeqLink === 'string')
    return (
      <SeqLink key={key} href={String(props.href)} id={props.dataSeqLink} runner={runner}>
        {kids}
      </SeqLink>
    );
  if (n.tagName === 'a')
    return (
      <a key={key} href={String(props.href ?? '')} target="_blank" rel="noopener noreferrer" style={LINK}>
        {kids}
      </a>
    );
  return createElement(n.tagName, { key }, ...(VOID.has(n.tagName) ? [] : kids));
}

/** One rendered document: its own prefix for the link ids, so the same link in another description is another channel. */
function Markdown({ text, prefix, runner }: { text: string; prefix: string; runner: LinkRunner }) {
  const [tree, setTree] = useState<HastNode | null>(null);
  useEffect(() => {
    let live = true;
    void markdownToHast(text, prefix).then((t) => live && setTree(t));
    return () => {
      live = false;
    };
  }, [text, prefix]);
  return tree ? toReact(tree, 0, runner) : null;
}

/** The host's `Description`: markdown, interactive. A DOM part (the slot's contract) that renders React into its node. */
export const MarkdownDescription: DomPart<PartProps> = (host, ctx) => {
  const box = document.createElement('div');
  box.className = 'seqv-description';
  host.append(box);
  const root = createRoot(box);
  const runner = createLinkRunner((inv) => ctx.commands.run(inv, 'host:links'));
  const prefixes = new Map<string, string>();
  const render = () => {
    const d = descriptionOf(ctx);
    if (!d) return root.render(null);
    const prefix = prefixes.get(d.text) ?? `d${prefixes.size}-`;
    prefixes.set(d.text, prefix);
    root.render(<Markdown text={d.text} prefix={prefix} runner={runner} />);
  };
  const subs = [
    ctx.document.shown.subscribe(render),
    ctx.document.base.subscribe(render),
    ctx.snapshots.current.subscribe(render),
  ];
  return {
    dispose() {
      for (const s of subs) s.unsubscribe();
      runner.dispose();
      // Deferred: the viewer can dispose its parts while React renders it (a part swapped in a render), and a root
      // unmounted synchronously then warns ("Attempted to synchronously unmount a root while React was already rendering").
      queueMicrotask(() => {
        root.unmount();
        box.remove();
      });
    },
  };
};
