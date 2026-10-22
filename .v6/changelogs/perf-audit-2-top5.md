## Performance: the top wins of perf audit 2

- Model: a polymer residue without a backbone atom (a 5' nucleotide without P, as in 1tup) no longer reads coordinates
  at index −1 while its ranges are built. That read made V8's `Vec3` arrays generic for the rest of the session. 1tup
  ×32 (186k atoms): molecular surface 1.52–1.57 → 0.21–0.25 GB allocated per build, unit boundaries 50–56 → 16–33 MB
- Graphics: interactions impostors size their buffers once, by cylinders. The count was multiplied by six twice, so
  the first chunk was six times too large and was then copied. `repr:interactions` allocation: 1tup ×32 148–163 →
  76–89 MB, asm-192 275–292 → 72–81 MB; asm-192 build 429–505 → 268–290 ms
- Model: contacts are tested only for the feature type pairs a tester accepts (`ContactTester.acceptsTypes`, optional),
  and each feature is queried only as far as its type's partners need. Contacts are the same, in the same order.
  Interactions: 1tup ×32 0.71–0.73 → 0.41–0.54 s, asm-192 2.3–2.7 → 1.1–1.5 s
- Model: a moved structure (`Structure.transform`: superposition, Prism's structure move) and its components restrict
  the root's inter-unit bonds instead of computing their own at every move. If the root has none yet, it computes them
  once. asm-192: 320–730 → 0.6–1.8 ms for the moved structure, 310–540 → 17–38 ms for its polymer component
- State: a builder commit that removes several nodes (or one with dependents) is diffed from its edits, not from both
  whole trees. Re-applying a representation preset to 500 structures: 4.1–15.1 → 0.5–0.7 s, 3.1 → 0.4 GB allocated
