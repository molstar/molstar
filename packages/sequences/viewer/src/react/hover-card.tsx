/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The hover card's React part (seqspec:interaction#the-hover-card): the same controller and placement as the DOM part.

import { useEffect, useLayoutEffect, useMemo, useRef, useSyncExternalStore } from 'react';
import { dismissOn, placeCard } from '../dom/hover-card.js';
import { createHoverCard } from '../hover-card.js';
import { type PartProps, subscriber, useScopeMaybe } from './hooks.js';
import { NO_SURFACES } from './parts.js';

export function HoverCard({ context: ctx }: PartProps) {
  // A canvas overlay, placed at its UI's surfaces: nothing without a provider (as the DOM's without a UI).
  const scope = useScopeMaybe();
  const surfaces = scope?.surfaces ?? NO_SURFACES;
  const card = useMemo(() => createHoverCard(ctx), [ctx]);
  useEffect(() => () => card.dispose(), [card]);
  // The state's replay on subscribe is not a change.
  const subscribe = useMemo(() => subscriber([card.current]), [card]);
  const s = useSyncExternalStore(subscribe, () => card.current.value);
  const el = useRef<HTMLDivElement>(null);
  useEffect(() => dismissOn(card, () => el.current), [card]);
  useLayoutEffect(() => {
    if (s && el.current) placeCard(surfaces, el.current, s);
  }, [s, surfaces]);
  if (!s || !scope) return null;
  return (
    <div
      ref={el}
      role="dialog"
      aria-label={s.label?.title ?? 'Item'}
      className="seqv-hover-card"
      onPointerEnter={() => card.enter()}
      onPointerLeave={() => card.leave()}
    >
      {s.label?.title && <div className="seqv-tooltip-title">{s.label.title}</div>}
      {s.label?.lines.map((l) => (
        <div key={`${l.key}:${l.value}`} className="seqv-tooltip-line">
          <span className="seqv-quiet">{l.key}</span> {l.value}
        </div>
      ))}
      <div className="seqv-hover-card-actions">
        {s.actions.map((a) => (
          <button key={a.name} type="button" className="seqv-hover-card-action" onClick={() => card.run(a.name)}>
            {a.title}
          </button>
        ))}
      </div>
    </div>
  );
}
