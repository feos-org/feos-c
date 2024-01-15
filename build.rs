use std::path::PathBuf;

fn main() {
    let crate_dir = std::env::var("CARGO_MANIFEST_DIR").unwrap();

    let mut config: cbindgen::Config = Default::default();
    config.language = cbindgen::Language::C;
    config.cpp_compat = true;
    config.include_version = false;
    config.documentation = true;
    config.documentation_style = cbindgen::DocumentationStyle::Doxy;

    let result = cbindgen::Builder::new()
        .with_crate(crate_dir)
        .with_config(config)
        .generate()
        .map(|data| {
            let mut path = PathBuf::from("include");
            path.push("feos.h");
            data.write_to_file(&path);
        });

    if result.is_ok() {
        println!("cargo:rerun-if-changed=src");
    } else {
    }
}
