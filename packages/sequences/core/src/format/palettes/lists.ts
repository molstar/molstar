/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

/**
 * MVS's named color lists (`ColorListNameT`), keyed by MVS's names.
 *
 * Copied from Mol* (MIT, (c) 2018-2025 mol* contributors; Mol*'s src): the lists of
 * `mol-util/color/lists.ts` (authors Alexander Rose, Adam Midlik) that `extensions/mvs/helpers/colors.ts`
 * (`MvsNamedColorLists`, Adam Midlik) names, in its order. The ColorBrewer lists are (c) 2002 Cynthia Brewer, Mark
 * Harrower and The Pennsylvania State University, Apache License 2.0; others from Matplotlib and D3.js. Written as
 * hex digits and read once at load.
 */

import { colorList } from "./color.js";

const list = colorList;

export const ColorLists = {
  Reds: list("fff5f0 fee0d2 fcbba1 fc9272 fb6a4a ef3b2c cb181d a50f15 67000d"),
  Oranges: list("fff5eb fee6ce fdd0a2 fdae6b fd8d3c f16913 d94801 a63603 7f2704"),
  Greens: list("f7fcf5 e5f5e0 c7e9c0 a1d99b 74c476 41ab5d 238b45 006d2c 00441b"),
  Blues: list("f7fbff deebf7 c6dbef 9ecae1 6baed6 4292c6 2171b5 08519c 08306b"),
  Purples: list("fcfbfd efedf5 dadaeb bcbddc 9e9ac8 807dba 6a51a3 54278f 3f007d"),
  Greys: list("ffffff f0f0f0 d9d9d9 bdbdbd 969696 737373 525252 252525 000000"),
  OrRd: list("fff7ec fee8c8 fdd49e fdbb84 fc8d59 ef6548 d7301f b30000 7f0000"),
  BuGn: list("f7fcfd e5f5f9 ccece6 99d8c9 66c2a4 41ae76 238b45 006d2c 00441b"),
  PuBuGn: list("fff7fb ece2f0 d0d1e6 a6bddb 67a9cf 3690c0 02818a 016c59 014636"),
  GnBu: list("f7fcf0 e0f3db ccebc5 a8ddb5 7bccc4 4eb3d3 2b8cbe 0868ac 084081"),
  PuBu: list("fff7fb ece7f2 d0d1e6 a6bddb 74a9cf 3690c0 0570b0 045a8d 023858"),
  BuPu: list("f7fcfd e0ecf4 bfd3e6 9ebcda 8c96c6 8c6bb1 88419d 810f7c 4d004b"),
  RdPu: list("fff7f3 fde0dd fcc5c0 fa9fb5 f768a1 dd3497 ae017e 7a0177 49006a"),
  PuRd: list("f7f4f9 e7e1ef d4b9da c994c7 df65b0 e7298a ce1256 980043 67001f"),
  YlOrRd: list("ffffcc ffeda0 fed976 feb24c fd8d3c fc4e2a e31a1c bd0026 800026"),
  YlOrBr: list("ffffe5 fff7bc fee391 fec44f fe9929 ec7014 cc4c02 993404 662506"),
  YlGn: list("ffffe5 f7fcb9 d9f0a3 addd8e 78c679 41ab5d 238443 006837 004529"),
  YlGnBu: list("ffffd9 edf8b1 c7e9b4 7fcdbb 41b6c4 1d91c0 225ea8 253494 081d58"),
  Magma: list(`
    000004 02020b 050416 090720 0e0b2b 140e36 1a1042 21114d 29115a 331067 3b0f70 440f76 4c117a 54137d 5b157f
    641a80 6b1d81 752181 7c2382 842681 8c2981 942c80 9c2e7f a5317e ad347c b73779 bf3a77 c73d73 ce3f70 d6456b
    de4968 e44f64 ea5661 ef5d5e f4675c f7705c f9795d fb835f fc8c63 fd9668 fe9f6d fea973 feb47b febd82 fec68a
    fecf92 fed89a fde2a3 fdebac fcf4b6 fcfdbf
  `),
  Inferno: list(`
    000004 02020c 050417 0a0722 10092d 160b39 1e0c45 250c51 2f0a5b 390963 420a68 4a0c6b 520e6d 5a116e 62136e
    6a176e 721a6e 7c1d6d 84206b 8c2369 932667 9b2964 a32c61 ab2f5e b3325a bc3754 c33b4f ca404a cf4445 d74b3e
    dd513a e25734 e75e2e eb6628 f06f20 f37819 f68013 f8890c fa9207 fb9b06 fca50a fcae12 fbba1f fac42a f8cd37
    f6d746 f4e156 f2ea69 f2f27d f5f992 fcffa4
  `),
  Plasma: list(`
    1b068c 250591 2f0495 380499 40039c 49029f 5101a2 5901a4 6100a6 6800a7 7000a8 7801a8 7f03a7 8607a6 8e0ca4
    9511a1 9b179e a21c9a a82296 ae2791 b42d8d b93388 be3883 c33e7f c8447a cd4975 d14f71 d6556d da5a68 de6064
    e26660 e56c5b e97257 ec7853 ef7e4e f2854a f58b46 f79241 f9993d faa039 fca735 fdaf31 fdb62d fdbe29 fdc626
    fcce25 fad624 f8df24 f5e726 f2f026 eff821
  `),
  Viridis: list(`
    45085b 470f62 471669 481d6f 482374 472a79 46307d 453681 433c84 414286 3e4888 3c4d8a 3a538b 37588c 355d8c
    32628d 30678d 2e6c8e 2c718e 2a768e 287a8e 267f8e 24848d 23898d 218d8c 1f928c 1e978a 1e9b89 1ea087 20a585
    23a982 28ae7f 2eb27c 35b778 3dbb74 45bf6f 4fc369 59c764 64cb5d 70ce56 7cd24f 88d547 95d73f a2da37 afdc2e
    bdde26 cae01e d7e219 e4e318 f1e51c fde724
  `),
  Cividis: list(`
    002051 002459 00275f 012b64 052e68 0a326a 10366c 15396d 1d3d6e 24406e 2b446e 32486e 394b6e 404f6e 45526d
    4d566d 535a6e 595e6e 5e616e 64656f 696970 6e6c71 727072 777473 7b7874 7f7c75 837f76 888377 8c8677 908b78
    948f78 999378 9d9778 a29b78 a7a077 ada476 b2a874 b8ac72 beb170 c4b56d caba6a d1be66 d7c363 dec85f e4cc5a
    ead156 efd652 f4db4e f8e04b fbe548 fdea45
  `),
  Turbo: list(`
    23171b 362356 422f85 493caa 4b4ac7 4a58dd 4666eb 4174f4 3b82f8 3590f8 2f9df5 2aaaef 27b6e6 25c2dc 25ccd0
    27d7c4 2be0b7 31e8aa 38ee9d 42f490 4df884 59fb78 67fd6d 76fe63 85fd59 95fb51 a4f749 b4f242 c2ec3c d1e436
    dedd32 ead32e f4c92b fcbd28 ffb125 ffa423 ff9721 ff891e ff7b1c fd6d1a f65f18 ec5115 e14412 d5370f c72c0c
    ba2208 ad1905 a21202 980e00 920c00 900c00
  `),
  Warm: list(`
    6e40aa 753fad 7d3faf 863eb1 8e3eb2 963db3 9e3db3 a73cb3 af3cb2 b73cb1 bf3caf c73dac cf3da9 d63ea6 dc3fa2
    e4419d ea4299 f04494 f5468e fa4988 fe4b83 ff4e7c ff5276 ff5670 ff5a6a ff5e63 ff635d ff6757 ff6c51 ff724b
    ff7847 ff7d42 ff833d ff8a39 ff9036 fb9633 f69d31 f2a42f edaa2e e7b12e e2b72f dcbe30 d7c432 d1ca34 ccd038
    c6d63c c1dc41 bce146 b7e64c b3eb53 aff05b
  `),
  Cool: list(`
    6e40aa 6c43b1 6947b7 674bbd 6450c3 6054c8 5d59cd 595ed1 5563d5 5169d9 4c6edb 4874de 437ae0 3f80e1 3a85e1
    368ce1 3292e1 2e98df 2a9fde 27a5db 23abd8 21b1d5 1eb7d1 1cbccc 1bc2c7 1ac7c2 19cdbc 19d1b6 1ad5b0 1bdbaa
    1ddfa3 20e29c 23e696 27e98f 2bec89 30ef82 36f17c 3cf276 43f471 4af56c 52f667 5af663 63f75f 6cf65c 75f65a
    7ff658 88f557 92f457 9cf357 a6f159 aff05b
  `),
  CubehelixDefault: list(`
    000000 080308 0f0612 140a1c 180f26 1a1530 1b1c39 1a2241 192b47 17344b 163d4e 15464e 154f4d 17574a 195e46
    1f6642 266c3c 2f7137 397533 467830 54792f 637a2f 737b32 827a37 927a3f a07949 ad7955 b97963 c27a71 ca7c82
    d07e93 d382a3 d487b3 d48dc1 d294ce cf9cda cca5e3 c9aeea c5b8ef c3c1f2 c1caf3 c2d3f3 c3dbf2 c7e2f1 cce8f0
    d2eeef daf2ef e3f6f0 ecf9f3 f6fcf8 ffffff
  `),
  Rainbow: list(`
    6e40aa 7d3faf 8e3eb2 9e3db3 af3cb2 bf3caf cf3da9 dc3fa2 ea4299 f5468e fe4b83 ff5276 ff5a6a ff635d ff6c51
    ff7847 ff833d ff9036 f69d31 edaa2e e2b72f d7c432 ccd038 c1dc41 b7e64c aff05b 9cf357 88f557 75f659 62f75f
    52f667 43f471 36f17c 2bec89 23e696 1ddfa3 1ad6b0 19cdbc 1bc2c7 1eb7d1 23abd8 2a9fde 3292e1 3a86e1 437ae0
    4c6edb 5563d5 5d59cd 6450c3 6947b7 6e40aa
  `),
  Sinebow: list(`
    ff4040 fe4e32 fb5d26 f66d1c ef7d12 e78d0b dc9d05 d1ab02 c4bb00 b6c800 a7d503 97e007 88ea0d 77f215 68f81e
    58fc2a 49ff37 3bff44 2efd53 23fa62 18f472 10ed82 09e392 04d9a2 01cdb1 00bfbf 01b1cd 04a2d9 0892e2 1081ed
    1872f4 2362fa 2e53fd 3b44ff 4937ff 582afc 681ff8 7715f2 880dea 9707e0 a703d5 b600c8 c400bb d102ac dc059d
    e70b8d ef127d f61c6d fb265d fe324e ff4040
  `),
  RdBu: list("67001f b2182b d6604d f4a582 fddbc7 f7f7f7 d1e5f0 92c5de 4393c3 2166ac 053061"),
  RdGy: list("67001f b2182b d6604d f4a582 fddbc7 ffffff e0e0e0 bababa 878787 4d4d4d 1a1a1a"),
  PiYG: list("8e0152 c51b7d de77ae f1b6da fde0ef f7f7f7 e6f5d0 b8e186 7fbc41 4d9221 276419"),
  BrBG: list("543005 8c510a bf812d dfc27d f6e8c3 f5f5f5 c7eae5 80cdc1 35978f 01665e 003c30"),
  PRGn: list("40004b 762a83 9970ab c2a5cf e7d4e8 f7f7f7 d9f0d3 a6dba0 5aae61 1b7837 00441b"),
  PuOr: list("2d004b 542788 8073ac b2abd2 d8daeb f7f7f7 fee0b6 fdb863 e08214 b35806 7f3b08"),
  RdYlGn: list("a50026 d73027 f46d43 fdae61 fee08b ffffbf d9ef8b a6d96a 66bd63 1a9850 006837"),
  RdYlBu: list("a50026 d73027 f46d43 fdae61 fee090 ffffbf e0f3f8 abd9e9 74add1 4575b4 313695"),
  Spectral: list("9e0142 d53e4f f46d43 fdae61 fee08b ffffbf e6f598 abdda4 66c2a5 3288bd 5e4fa2"),
  Category10: list("1f77b4 ff7f0e 2ca02c d62728 9467bd 8c564b e377c2 7f7f7f bcbd22 17becf"),
  Observable10: list("4269d0 efb118 ff725c 6cc5b0 3ca951 ff8ab7 a463f2 97bbf5 9c6b4e 9498a0"),
  Tableau10: list("4e79a7 f28e2c e15759 76b7b2 59a14f edc949 af7aa1 ff9da7 9c755f bab0ab"),
  Set1: list("e41a1c 377eb8 4daf4a 984ea3 ff7f00 ffff33 a65628 f781bf 999999"),
  Set2: list("66c2a5 fc8d62 8da0cb e78ac3 a6d854 ffd92f e5c494 b3b3b3"),
  Set3: list("8dd3c7 ffffb3 bebada fb8072 80b1d3 fdb462 b3de69 fccde5 d9d9d9 bc80bd ccebc5 ffed6f"),
  Pastel1: list("fbb4ae b3cde3 ccebc5 decbe4 fed9a6 ffffcc e5d8bd fddaec f2f2f2"),
  Pastel2: list("b3e2cd fdcdac cbd5e8 f4cae4 e6f5c9 fff2ae f1e2cc cccccc"),
  Dark2: list("1b9e77 d95f02 7570b3 e7298a 66a61e e6ab02 a6761d 666666"),
  Paired: list("a6cee3 1f78b4 b2df8a 33a02c fb9a99 e31a1c fdbf6f ff7f00 cab2d6 6a3d9a ffff99 b15928"),
  Accent: list("7fc97f beaed4 fdc086 ffff99 386cb0 f0027f bf5b17 666666"),
  Chainbow: list(`
    4a41b5 4a58dd 426ff2 3987f9 2f9df5 28b2e9 25c6d8 27d7c4 2ee5ae 3bf098 4df884 62fd70 7bfe5f 95fb51 aff444
    c8ea3a dedd32 f0cc2c feb927 ffa423 ff8e1f ff761c f65f18 e54813 d0330e ba2208 a51403 960d00
  `),
};

export type ColorListName = keyof typeof ColorLists;
export const ColorListNames = Object.keys(ColorLists) as readonly ColorListName[];
