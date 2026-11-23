# Lane XD: the extension openings (strong, design only)

**Task**: draft `drafts/extension-openings.md` (shape: [design-drafts.md](../../third-slice/briefs/design-drafts.md)),
one draft in three parts, each cut by Gate 0's answer 2 (cut order CN, then LE):
1. The extension descriptor and `ext:` layer kinds with their fallback ([extensibility.md](../../../design/extensibility.md),
   EC-235; GR's research `../../third-slice/research-genome-ring.md`: the Rondo-like ring is the fifth slice's
   showcase as an extension package, so this must be enough for it): what an extension ships (layer types, layouts,
   sizings, behaviors, themes, formats), how a spec names it, what a viewer without it draws (the fallback), versions.
   MV made layer types, layouts and sizings registry entries; build on that.
2. Custom navigators: a navigator's picture as a registry entry with the window public (`render/window.ts`,
   `NavigatorModel`, navigators.md).
3. Legends' library surface: placement, hiding, a compact form; editing stays the application's (themes.md).

**Reads**: extensibility.md, viewer-composition.md, navigators.md, themes.md (legends), layouts.md (track sizing's
`ext:` names and the port), ports.md; core `src/context/`, `src/render/` registry code as built.

**Output**: `.seq/execution/fourth-slice/drafts/extension-openings.md` only. No code, no commits. DR and SP draft
beside you in the same folder: touch only your file.

**Model**: strong (opus). **Where**: the main checkout, read only except that file.

**Report**: under 300 words per part: the decision, the public surface, the gaps with recommendations.
