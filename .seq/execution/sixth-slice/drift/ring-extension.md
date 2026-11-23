# Drift: two things `ring-extension` designs that `ext-rondo` does not do

**What**:

1. **The ribbons' "limit" issue.** The design says ribbons past `rondo.maxRibbons` are drawn in their stems' fans and
   said once (`rondo.ribbons.limit`, info: "612 links drawn in fans"). The code has no such issue: the fan is drawn and
   the card lists it, nothing is reported. The other `rondo.ribbons.*` issues exist (`layout`, `colors`, `filter`,
   `color`, `width`).
2. **The hue ring without bands.** The design says that with no `bands`, `rondoRing` writes a ranges annotation of whole
   sequences from the assembly's inline `sequences` and the hue ring over it. The code writes the bands and hue tracks
   only when `bands` is given; with none, a ring has no hue ring.

**Where**: `packages/seq/ext-rondo/src/ribbons.ts` (the issues its layout reports), `ring.ts` (`ringTracks`).

**Options**: (1) (a) add the info issue; (b) drop it from the design (the legend and the card already say what is in
fans). (2) (a) write the annotation (needs the assembly's sequences at call time: only inline ones, so a ring over a
source assembly could not); (b) amend the design: the hue ring is the bands' second layer, absent without bands, a host
that wants one passes a ranges node.

**Recommendation**: (1) (b): a fan is a normal drawing at the width rule, not a loss, and an info issue on every busy
ring is noise. (2) (b): inferring a node from inline lengths only works for one kind of assembly and hides a data
choice in a helper.

## Added in the fix round

3. **A hit is not `(rect :a :b)`.** The design says a ribbon is hit as its drawn run, `(rect :a ... :b ...)`, a stub as
   the stem's range by the partner ends' hull. The code's hit is `{ element: { owner, key }, system, position, bin, count }`
   (`shapes.ts`, `hitsOf`), keyed `ribbon:`, `stem:`, `stub:`; no rect is built. Options: (a) amend the design to the
   key hit; (b) build the rect (a bookmark of a ribbon needs one). Recommendation: (a) until a bookmark of a ribbon is
   asked for.
4. **"No lazy chunk inside the package."** The viewer entry lazy-loads `menu.js` and `card.js` (`viewer.ts`) once
   composed; the design says the package is itself the opt-in. Recommendation: amend the design (the chunk keeps the
   headless entry free of the card and the menu).
5. **`rondoRing` does not write the spec's `extensions`.** The design says it writes `extensions: { rondo: '0.1' }`; the
   code leaves it to the writer (`createSeqSpec({ extensions: [rondo] })`), as the spec says. Recommendation: amend the
   design; the writer owns what the spec requires.
