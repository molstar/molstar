/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

/**
 * The spheres' second color (`SiteColorSchema`) blended into `out.vColor` by the sphere's radius on screen in CSS pixels
 * (as the levels of detail measure it): from none at `uSiteColorRange.x` to whole at `.y`, a table entry's alpha 0
 * none (a carbon its particle's). Expects `out.vColor` (`assign_varyings_wgsl`, its `groupIndex`), `radius` and
 * `mvPosition`, and for `siteType` the site (`site`) and, merged, `tableBase` its table's first entry; not
 * for palette colors (indices, not colors)
 */
export const assign_site_color_wgsl = (tableBase: string) => /* wgsl */ `
#if defined(dSiteColor_siteType) || defined(dSiteColor_group)
#if defined(dColorType_texture) && !defined(dUsePalette)
#if defined(dRenderVariant_color) || defined(dRenderVariant_tracing) || defined(dRenderVariant_emissive)
    {
        let siteColorDepth = select(max(-mvPosition.z, uGlobal.uNear), 1.0, uGlobal.uIsOrtho == 1.0);
        let siteColorPixels = radius * uGlobal.uProjection[1][1] * 0.5 * uGlobal.uCssHeight / siteColorDepth;
        let siteColorRange = uMaterial.uSiteColorRange;
        let siteColorT = smoothstep(siteColorRange.x, max(siteColorRange.y, siteColorRange.x + 0.001), siteColorPixels);
        if (siteColorT > 0.0) {
            #if defined(dSiteColor_siteType) && defined(dSegmented)
                let siteColor = read_tSiteColor(${tableBase} + i32(site.w));
            #elif defined(dSiteColor_siteType)
                let siteColor = read_tSiteColor(i32(site.w));
            #else
                let siteColor = read_tSiteColor(groupIndex);
            #endif
            out.vColor = vec4<f32>(mix(out.vColor.rgb, siteColor.rgb, siteColorT * siteColor.a), out.vColor.w);
        }
    }
#endif
#endif
#endif
`;
