## Model: faster inter-unit interactions

- The inter-unit contacts loop takes the features of a unit near its partner from the unit's feature grid instead of
  transforming every feature for every unit pair; contacts and their order are unchanged. 1tup 2×2×2 cell blocks
  (186k atoms, 288 units): loop 333 → 316 ms; as one model of 416 units: 330 → 318 ms
