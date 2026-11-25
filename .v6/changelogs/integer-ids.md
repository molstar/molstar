## Graphics: integer instance ids

- Instance ids and data indices as integers on WebGL2 and WebGPU (`integer-ids.md`): marking, colors and other data of
  instances and groups past 2^24 (16.7 M; the yeast cell's 31.5 M instances) read their own, not a neighbour's.
- Picking instances past 2^24: the instance id in 32 bits (the pick target's alpha its top byte).
