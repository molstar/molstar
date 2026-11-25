## Seq* (next): a viewer ships only what its spec lists

- Layer types and section layouts are registry entries (`layer-type`, `layout`): the canvas and the SVG export draw the
  composition's own, and say what they cannot draw; tracks is the one fixed layout.
- Navigators (`NavigatorsPart`), story steps and transitions (`StoryPart`), the mounted history (`ui.history`, a
  factory) and the rows rail are parts and behaviors the default spec lists. `DefaultSeqViewerSpec()` and
  `createSeqContext()` (now `@molstar/seqx-core/context/default`) list everything; `composeSeqContext` and `exportSvg`
  compose only what they are given.
- The default compositions moved off the barrels to subpaths: `render/render-svg`, `render/export-png`,
  `render/navigator`, `render/built-ins`. The seq core and viewer packages say `sideEffects: false`.
- A minimal viewer (tracks, letters, ranges) drops from 1014 to 830 kB minified.
