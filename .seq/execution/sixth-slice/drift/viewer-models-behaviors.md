# Drift: the first behaviors list (viewer-models)

**What**: `viewer-models.md`, "Behaviors (the first list)", names `follow` (keeps followers in step) and `focus`
(focus and reveal) as behaviors with their params as config items. The code has neither: following is the
`ViewModel`'s own (`follow.followerMoves`, `follow.cyclicArc` are `ViewerConfig` items the model reads; viewer spec,
ViewModel, Follow), and focus and reveal are commands (`address.focus`, `address.reveal` of `addressing`; the
`navigate` behavior's double-click focus). The behaviors built (the providers' names in `behaviors/`) are
`addressing`, `cursor`, `editing`, `headers`, `hover`, `hoverCard`, `keysPanel`, `legend`, `linkedPairs`,
`links.hover`, `navigate`, `navigator`, `rows`, `select`, `window`.

**Where**: `packages/seq/viewer/src/behaviors/`, `models/view.ts`.

**Options**: (a) amend the design: follow is a model rule (a viewport derived through a map must hold with no behavior
listed, or a follower drifts), focus a command; (b) split a `follow` behavior out of `ViewModel` so a host can drop
it.

**Recommendation**: (a). A follow is part of what the view is, not a policy a host opts out of; the spec already says
so.
