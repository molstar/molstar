/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

/**
 * The spheres' second color (`SiteColorSchema`) blended into `vColor` by the sphere's radius on screen in CSS pixels
 * (as the levels of detail measure it): from none at `uSiteColorRange.x` to whole at `.y`, a table entry's alpha 0
 * none (a carbon its particle's). Expects `vColor` (`assign_color_varying`), `vSphere.w` the radius and `mvPosition`,
 * and for `siteType` the site (`site`; merged, its table's first entry `aSegmentItem.w`); not for palette colors. The
 * radius without the instance's scale (a coarse level's particles drawn up to 63 times larger to fill a far copy's
 * space, their atoms not near by that)
 */
export const assign_site_color = `
#if defined(dSiteColor_siteType) || defined(dSiteColor_group)
#if defined(dColorType_texture) && !defined(dUsePalette)
#if defined(dRenderVariant_color) || defined(dRenderVariant_tracing) || defined(dRenderVariant_emissive)
    {
        float siteColorDepth = uIsOrtho == 1.0 ? 1.0 : max(-mvPosition.z, uNear);
        float siteColorRadius = vSphere.w / max(length(aTransform[0].xyz), 1e-6);
        float siteColorPixels = siteColorRadius * uProjection[1][1] * 0.5 * uCssHeight / siteColorDepth;
        float siteColorT = smoothstep(uSiteColorRange.x, max(uSiteColorRange.y, uSiteColorRange.x + 0.001), siteColorPixels);
        if (siteColorT > 0.0) {
            #if defined(dSiteColor_siteType)
                #if defined(dSegmented)
                    vec4 siteColor = read_tSiteColor(int(aSegmentItem.w) + int(site.w));
                #else
                    vec4 siteColor = read_tSiteColor(int(site.w));
                #endif
            #else
                vec4 siteColor = read_tSiteColor(groupIndex(group));
            #endif
            vColor.rgb = mix(vColor.rgb, siteColor.rgb, siteColorT * siteColor.a);
        }
    }
#endif
#endif
#endif
`;
