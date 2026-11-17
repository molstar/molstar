# Drift: the protocol's details lane D settled

Where: protocol.md, data-store.md. Lane D's choices, each with its recommendation (the orchestrator agrees unless the
review says otherwise):

- A narrow asset reader (`readText(uri, signal)`, later `read(uri, { range, signal })`) instead of `AssetManager` in the
  open context; the default adapter wraps `AssetManager` behind the policy (ports.md, security.md). Edit protocol.md.
- `SourceEntry` minimal (`address, format, uri?, value?`) until lane E binds it to the language's source node.
- Data keys without shapes: one entry holds every shape the coverage lists for that system. Keep.
- A range crossing a tile edge appears in each tile; draw lists dedupe by id. Say so in layout-and-draw-lists.md.
- Tables carry a `form` field (rows or density; per position or per bin) and extra fields in a `fields` record. Add to
  protocol.md.
- The seq file is a tree value through the core JSON codec; a failure's `kind` is spelled `failure`.
- Malformed answers: recommended not retried (deterministic); lane D retried them. To fix.
- `seq.store.over_budget` (underscore, as repr's codes); edit protocol.md.
- Eviction may drop wanted margin or prefetch tiles, refetched on the next `want`; coarser levels are the planner's.

After lane D's review (the orchestrator's decisions, sent in the fix round): malformed answers are not retried (failed
with a cooldown, EC-174); the ranges form `rows` is renamed `list`; ranges sliced across tiles carry a stable row
identity; the store has per-consumer want sets (`wanter()`) so exports never touch the live view's wants (EC-125); when
view work waits on a full in-flight bound the oldest stale request is cancelled at once (EC-173). protocol.md and
data-store.md are edited at the merge.
