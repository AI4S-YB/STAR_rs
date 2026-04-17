//! 1:1 port of `Parameters::samAttributes` (Parameters_samAttributes.cpp).
//!
//! Scope for M3 basic SE: resolves the user `outSAMattributes` list into
//! a vector of attribute IDs, filling the subset of `OutSAMAttrPresent`
//! flags that the SAM output path checks. BAM-only attributes are still
//! flagged, but the "requires BAM" enforcement uses the same logic as C++.

use star_core::types::{
    ATTR_AS, ATTR_CB, ATTR_CH, ATTR_CN, ATTR_CR, ATTR_CY, ATTR_GN, ATTR_GN_LOWER, ATTR_GX,
    ATTR_GX_LOWER, ATTR_HA, ATTR_HI, ATTR_JI, ATTR_JM, ATTR_MC, ATTR_MD, ATTR_NH, ATTR_NM,
    ATTR_NM_LOWER, ATTR_RB, ATTR_RG, ATTR_SF, ATTR_SM, ATTR_SQ, ATTR_SS, ATTR_UB, ATTR_UR,
    ATTR_UY, ATTR_VA, ATTR_VG, ATTR_VW, ATTR_XS,
};

use crate::parameters::Parameters;

#[derive(Debug, Default, Clone)]
pub struct OutSamAttrPresent {
    pub nh: bool,
    pub hi: bool,
    pub as_: bool,
    pub nm: bool,
    pub md: bool,
    pub n_m: bool,
    pub j_m: bool,
    pub j_i: bool,
    pub rg: bool,
    pub mc: bool,
    pub xs: bool,
    pub ch: bool,
    pub v_a: bool,
    pub v_g: bool,
    pub v_w: bool,
    pub r_b: bool,
    pub ha: bool,
    pub cr: bool,
    pub cy: bool,
    pub ur: bool,
    pub uy: bool,
    pub cb: bool,
    pub ub: bool,
    pub gx: bool,
    pub gn: bool,
    pub gx_lower: bool,
    pub gn_lower: bool,
    pub s_m: bool,
    pub s_s: bool,
    pub s_q: bool,
    pub s_f: bool,
    pub c_n: bool,
}

#[derive(Debug, Default, Clone)]
pub struct SamAttributes {
    pub order: Vec<u16>,
    pub order_quant: Vec<u16>,
    pub present: OutSamAttrPresent,
}

impl SamAttributes {
    /// Port of `Parameters::samAttributes` (Parameters_samAttributes.cpp).
    pub fn resolve(p: &mut Parameters) -> anyhow::Result<Self> {
        let mut out = SamAttributes::default();

        out.order_quant.push(ATTR_NH);
        out.order_quant.push(ATTR_HI);

        let v_attr1: Vec<&str> = match p.out_sam_attributes.first().map(|s| s.as_str()) {
            Some("None") => vec![],
            Some("All") => vec!["NH", "HI", "AS", "nM", "NM", "MD", "jM", "jI", "MC", "ch"],
            Some("Standard") => vec!["NH", "HI", "AS", "nM"],
            _ => p.out_sam_attributes.iter().map(|s| s.as_str()).collect(),
        };

        for attr in &v_attr1 {
            match *attr {
                "NH" => {
                    out.order.push(ATTR_NH);
                    out.present.nh = true;
                }
                "HI" => {
                    out.order.push(ATTR_HI);
                    out.present.hi = true;
                }
                "AS" => {
                    out.order.push(ATTR_AS);
                    out.present.as_ = true;
                }
                "NM" => {
                    out.order.push(ATTR_NM);
                    out.present.nm = true;
                }
                "MD" => {
                    out.order.push(ATTR_MD);
                    out.present.md = true;
                }
                "nM" => {
                    out.order.push(ATTR_NM_LOWER);
                    out.present.n_m = true;
                }
                "jM" => {
                    out.order.push(ATTR_JM);
                    out.present.j_m = true;
                }
                "jI" => {
                    out.order.push(ATTR_JI);
                    out.present.j_i = true;
                }
                "vA" => {
                    out.order.push(ATTR_VA);
                    out.present.v_a = true;
                }
                "vG" => {
                    out.order.push(ATTR_VG);
                    out.present.v_g = true;
                }
                "vW" => {
                    out.order.push(ATTR_VW);
                    out.present.v_w = true;
                }
                "ha" => {
                    out.order.push(ATTR_HA);
                    out.present.ha = true;
                }
                "RG" => {
                    out.order.push(ATTR_RG);
                    out.order_quant.push(ATTR_RG);
                    out.present.rg = true;
                }
                "rB" => {
                    out.order.push(ATTR_RB);
                    out.order_quant.push(ATTR_RB);
                    out.present.r_b = true;
                }
                "ch" => {
                    out.order.push(ATTR_CH);
                    out.order_quant.push(ATTR_CH);
                    out.present.ch = true;
                }
                "MC" => {
                    out.order.push(ATTR_MC);
                    out.order_quant.push(ATTR_MC);
                    out.present.mc = true;
                }
                "CR" => {
                    out.order.push(ATTR_CR);
                    out.order_quant.push(ATTR_CR);
                    out.present.cr = true;
                }
                "CY" => {
                    out.order.push(ATTR_CY);
                    out.order_quant.push(ATTR_CY);
                    out.present.cy = true;
                }
                "UR" => {
                    out.order.push(ATTR_UR);
                    out.order_quant.push(ATTR_UR);
                    out.present.ur = true;
                }
                "UY" => {
                    out.order.push(ATTR_UY);
                    out.order_quant.push(ATTR_UY);
                    out.present.uy = true;
                }
                "CB" => {
                    out.order.push(ATTR_CB);
                    out.order_quant.push(ATTR_CB);
                    out.present.cb = true;
                }
                "UB" => {
                    out.order.push(ATTR_UB);
                    out.order_quant.push(ATTR_UB);
                    out.present.ub = true;
                }
                "GX" => {
                    out.order.push(ATTR_GX);
                    out.order_quant.push(ATTR_GX);
                    out.present.gx = true;
                }
                "GN" => {
                    out.order.push(ATTR_GN);
                    out.order_quant.push(ATTR_GN);
                    out.present.gn = true;
                }
                "gx" => {
                    out.order.push(ATTR_GX_LOWER);
                    out.order_quant.push(ATTR_GX_LOWER);
                    out.present.gx_lower = true;
                }
                "gn" => {
                    out.order.push(ATTR_GN_LOWER);
                    out.order_quant.push(ATTR_GN_LOWER);
                    out.present.gn_lower = true;
                }
                "sM" => {
                    out.order.push(ATTR_SM);
                    out.order_quant.push(ATTR_SM);
                    out.present.s_m = true;
                }
                "sS" => {
                    out.order.push(ATTR_SS);
                    out.order_quant.push(ATTR_SS);
                    out.present.s_s = true;
                }
                "sF" => {
                    out.order.push(ATTR_SF);
                    out.order_quant.push(ATTR_SF);
                    out.present.s_f = true;
                }
                "sQ" => {
                    out.order.push(ATTR_SQ);
                    out.order_quant.push(ATTR_SQ);
                    out.present.s_q = true;
                }
                "cN" => {
                    out.order.push(ATTR_CN);
                    out.order_quant.push(ATTR_CN);
                    out.present.c_n = true;
                }
                "XS" => {
                    out.order.push(ATTR_XS);
                    out.present.xs = true;
                    if p.out_sam_strand_field_type != 1 {
                        p.out_sam_strand_field_type = 1;
                    }
                }
                other => {
                    anyhow::bail!(
                        "EXITING because of FATAL INPUT ERROR: unknown/unimplemented SAM attribute (tag): {}",
                        other
                    );
                }
            }
        }

        if p.out_sam_strand_field_type == 1 && !out.present.xs {
            out.order.push(ATTR_XS);
        }

        Ok(out)
    }
}
