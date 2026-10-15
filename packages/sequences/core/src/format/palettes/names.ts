/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

/**
 * Mol*'s color names (X11 names, http://www.w3.org/TR/css3-color/#svg-color, with Mol*'s additions), the names a
 * document `Color` may use.
 *
 * Copied from Mol* (MIT, (c) 2018-2019 mol* contributors; Mol*'s src/mol-util/color/names.ts, author
 * Alexander Rose). Written as groups (a color, then its names) and read once at load.
 */

const src = `
    #f0f8ff aliceblue #faebd7 antiquewhite #00ffff aqua cyan #7fffd4 aquamarine #f0ffff azure #f5f5dc beige
    #ffe4c4 bisque #000000 black #ffebcd blanchedalmond #0000ff blue #8a2be2 blueviolet #a52a2a brown #deb887
    burlywood #5f9ea0 cadetblue #7fff00 chartreuse #d2691e chocolate #ff7f50 coral #6495ed cornflower
    cornflowerblue #fff8dc cornsilk #dc143c crimson #00008b darkblue #008b8b darkcyan #b8860b darkgoldenrod
    #a9a9a9 darkgray darkgrey #006400 darkgreen #bdb76b darkkhaki #8b008b darkmagenta #556b2f darkolivegreen
    #ff8c00 darkorange #9932cc darkorchid #8b0000 darkred #e9967a darksalmon #8fbc8f darkseagreen #483d8b
    darkslateblue #2f4f4f darkslategray darkslategrey #00ced1 darkturquoise #9400d3 darkviolet #ff1493 deeppink
    #00bfff deepskyblue #696969 dimgray dimgrey #1e90ff dodgerblue #b22222 firebrick #fffaf0 floralwhite #228b22
    forestgreen #ff00ff fuchsia magenta #dcdcdc gainsboro #f8f8ff ghostwhite #ffd700 gold #daa520 goldenrod
    #808080 gray grey #008000 green #adff2f greenyellow #f0fff0 honeydew #ff69b4 hotpink #cd5c5c indianred
    #4b0082 indigo #fffff0 ivory #f0e68c khaki #ffff54 laserlemon #e6e6fa lavender #fff0f5 lavenderblush #7cfc00
    lawngreen #fffacd lemonchiffon #add8e6 lightblue #f08080 lightcoral #e0ffff lightcyan #fafad2 lightgoldenrod
    lightgoldenrodyellow #d3d3d3 lightgray lightgrey #90ee90 lightgreen #ffb6c1 lightpink #ffa07a lightsalmon
    #20b2aa lightseagreen #87cefa lightskyblue #778899 lightslategray lightslategrey #b0c4de lightsteelblue
    #ffffe0 lightyellow #00ff00 lime #32cd32 limegreen #faf0e6 linen #800000 maroon #7f0000 maroon2 #b03060
    maroon3 #66cdaa mediumaquamarine #0000cd mediumblue #ba55d3 mediumorchid #9370db mediumpurple #3cb371
    mediumseagreen #7b68ee mediumslateblue #00fa9a mediumspringgreen #48d1cc mediumturquoise #c71585
    mediumvioletred #191970 midnightblue #f5fffa mintcream #ffe4e1 mistyrose #ffe4b5 moccasin #ffdead
    navajowhite #000080 navy #fdf5e6 oldlace #808000 olive #6b8e23 olivedrab #ffa500 orange #ff4500 orangered
    #da70d6 orchid #eee8aa palegoldenrod #98fb98 palegreen #afeeee paleturquoise #db7093 palevioletred #ffefd5
    papayawhip #ffdab9 peachpuff #cd853f peru #ffc0cb pink #dda0dd plum #b0e0e6 powderblue #800080 purple
    #7f007f purple2 #a020f0 purple3 #663399 rebeccapurple #ff0000 red #bc8f8f rosybrown #4169e1 royalblue
    #8b4513 saddlebrown #fa8072 salmon #f4a460 sandybrown #2e8b57 seagreen #fff5ee seashell #a0522d sienna
    #c0c0c0 silver #87ceeb skyblue #6a5acd slateblue #708090 slategray slategrey #fffafa snow #00ff7f
    springgreen #4682b4 steelblue #d2b48c tan #008080 teal #d8bfd8 thistle #ff6347 tomato #40e0d0 turquoise
    #ee82ee violet #f5deb3 wheat #ffffff white #f5f5f5 whitesmoke #ffff00 yellow #9acd32 yellowgreen
  `;

/** Lower-case name → 0xRRGGBB. */
export const ColorNames: ReadonlyMap<string, number> = (() => {
  const out = new Map<string, number>();
  let color = 0;
  for (const token of src.trim().split(/\s+/)) {
    if (token.startsWith('#')) color = Number.parseInt(token.slice(1), 16);
    else out.set(token, color);
  }
  return out;
})();
