## Seq* (next): stories

- Stories: steps of a view and a spec, stored as patches from keyframes (`story.maxChain`), with full-state steps; the
  `seq.story` channels (the view with `top` and bands, heights, colors blended in Oklab with theme tokens resolved,
  entering and leaving data) tweened by a plan a host can play or scrub. The entry `@molstar/seqx-core/story`.
- The viewer plays steps through an animation port (the host owns the loop); input by hand or a spec change ends a
  motion; reduced motion cuts (`story.reducedMotion` for stories, `view.motion` for jumps).
- The p53 story example, with an SVG per step.
