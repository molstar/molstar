/** @jsxImportSource @molstar/plugin-ui-next/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

/**
 * A stand-in for the 3D canvas in stories: a stylised chain on a backdrop, so layouts can be judged against something
 * molecule-like without WebGL.
 */
export function CanvasPlaceholder(props: { empty?: boolean; background?: string }) {
  const background = props.background ?? 'radial-gradient(circle at 50% 45%, #ffffff 0%, #eef0f3 70%, #dfe3e8 100%)';
  if (props.empty) return <div style={{ position: 'absolute', inset: 0, background }} />;
  const beads: SVGElement[] = [];
  for (let i = 0; i < 90; i++) {
    const t = i / 90;
    const x = 120 + t * 560 + 40 * Math.sin(t * 19);
    const y = 300 + 150 * Math.sin(t * 7.3) + 30 * Math.cos(t * 23);
    const hue = 200 + 160 * t;
    beads.push((<circle cx={x} cy={y} r={14 + 4 * Math.sin(t * 31)} fill={`hsl(${hue} 55% 60%)`} />) as SVGElement);
  }
  return (
    <div style={{ position: 'absolute', inset: 0, background }}>
      <svg viewBox="0 0 800 600" preserveAspectRatio="xMidYMid meet" style="width: 100%; height: 100%">
        <defs>
          <radialGradient id="msp-placeholder-shade" cx="35%" cy="35%" r="65%">
            <stop offset="0%" stop-color="#fff" stop-opacity="0.55" />
            <stop offset="100%" stop-color="#000" stop-opacity="0.2" />
          </radialGradient>
        </defs>
        <g>{beads}</g>
        <g style="mix-blend-mode: soft-light">
          {beads.map((b) => (
            <circle
              cx={b.getAttribute('cx')!}
              cy={b.getAttribute('cy')!}
              r={b.getAttribute('r')!}
              fill="url(#msp-placeholder-shade)"
            />
          ))}
        </g>
      </svg>
    </div>
  );
}
