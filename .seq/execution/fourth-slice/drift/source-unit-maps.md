# Drift: resolve cannot hold a map onto a source's unit

**What**: structure-syncing.md reads a source's map onto one unit (SIFTS per chain) as the route a structure viewer's
hover takes. Lane ST1 reads an overlay unit's own maps (`Overlay.edgesFrom`), but resolve refuses a map whose `to` is
a unit a source's description gave (`seq.resolve.ref.missing`), so ST1's test injects the overlay edge by hand.

**Options**:
1. Resolve takes a map onto a described unit as an overlay edge (DP's description gives the unit before the map
   resolves; the map waits as pending while the description is).
2. Such maps stay the host's: it maps onto the static molecule, and the unit is reached by the copy rule.

**Recommendation**: 1, in ST2 or HD (the fake host's adapter serves SIFTS per chain): the design's route, and the
pending state is the one EC-169 already gives.
