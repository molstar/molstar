# Drift: how a stored patch's refs are typed (lane K3's review)

**What**: patches.md types a stored `set` value through "the target's kind and key", with the kind read from the
step's dry run over the spec. So whether a value is a ref depends on the spec at that moment. Making law 3 hold under
that rule took lane K3 a retype of every step before each move, a loose pass for the inverse, and link integrity over
every step after any address change. The independent review (2026-11-18) found:
- **B1**: a forward refusal weakened: a move unrelated to a step switches off a dangling refusal the parent gave
  (`[remove /d, move /c → / rename d, remove /b]` with a step `set /d link "/b"`); fixing it breaks law 3 elsewhere
  unless the allowances also count a value typed in one root and not the other;
- **S1**: cost grows with moves × steps: 10 reorders on a p53-size spec with 30 steps 0.08 → 6.7 ms, with 120 steps
  25 ms; 10 renames with 120 steps 31 ms (the 16 ms edit budget);
- **S2**: nested stored patches (a step inside a step's inserted node) still break law 3 (4 in about 23,700 cases);
- **C1**: a host's unrelated edit can be refused because it changes a step's dry run (a step's own move now refused),
  while the common silent retarget (a step inserts `/x` and sets a key on it; the host inserts its own `/x`) passes.

**Where**: `packages/core/src/repr/patch/{apply,refs}.ts` on `seq/k3-stored-law3` (not merged).

**Options**:
1. **Type a stored `set` value without the spec** (recommended): a value is a ref when the target language types that
   key as a ref in any kind that has it (by kind and key across the target's kinds, as the inverse's loose pass
   already does), never by reading the dry run. Then typing never changes under a patch: no retype, no loose pass, no
   dependent chunks in the ref index, no per-move cost, nested steps by the same rule. The dry run stays for what a
   ref names (link integrity's `names`), and a step whose dry run changes is shown as not applying (option (a) below),
   not a refusal of the host's edit. Cost: a key that is a ref in one kind and not in another is taken as a ref
   (conservative: it follows moves and counts in link integrity).
2. **Keep typing by the dry run** and finish K3: B1 (judge by the new root's typing, count values typed in one root
   only), S1 (retype only dirty steps, skip reorders), S2 (nested steps by their own step, or declared opaque).
   More code, a cost per move that grows with steps, and C1 stays.

**And the step's meaning after a move** (the lane's question): (a) accept, a step that no longer applies is shown so
(recommended by the lane and the reviewer; matches `place.after` and changes.md); (b) type by the first op's state
(drops "insert, then set on it"); (c) refuse such moves (blocks host edits for story steps).

**Recommendation**: option 1 with (a). Lane K3 rebuilds on it (its hand cases and the fuzz kept); patches.md's
stored-patch rule and the repr `_spec` change; R+K's `dep` chunks go.

**Resolved** (the author, 2026-11-18): option 1 with (a), and stories redesigned around it (`stories.md`).
