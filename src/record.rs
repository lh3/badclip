//! The shared `extract` record format: helpers used by every subcommand that
//! reads `extract` output (`getfa`, `fltbnd`, `merge`), so they agree on it.

/// The breakend type of an `extract` record (its `ctype=` INFO tag).
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Ctype {
    Trim,
    Split,
    Del,
    Ins,
}

impl Ctype {
    /// The `ctype=` value of an INFO column. It is required: a record without
    /// one (`extract` output from before 0.2) or with an unknown value is an
    /// error.
    pub fn from_info(info: &str) -> Result<Ctype, String> {
        match info.split(';').find_map(|kv| kv.strip_prefix("ctype=")) {
            Some("trim") => Ok(Ctype::Trim),
            Some("split") => Ok(Ctype::Split),
            Some("del") => Ok(Ctype::Del),
            Some("ins") => Ok(Ctype::Ins),
            Some(v) => Err(format!("unknown ctype `{v}`")),
            None => Err("no `ctype=` tag (extract output from before badclip 0.2?)".to_string()),
        }
    }

    pub fn as_str(self) -> &'static str {
        match self {
            Ctype::Trim => "trim",
            Ctype::Split => "split",
            Ctype::Del => "del",
            Ctype::Ins => "ins",
        }
    }
}
