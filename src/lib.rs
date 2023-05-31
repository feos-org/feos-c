use feos::pcsaft::*;
use feos::EosVariant;
use feos_core::joback::JobackRecord;
use feos_core::parameter::BinaryRecord;
use feos_core::parameter::Identifier;
use feos_core::parameter::IdentifierOption;
use feos_core::parameter::Parameter;
use feos_core::parameter::PureRecord;
use feos_core::Contributions;
use feos_core::EquationOfState;
use feos_core::StateBuilder;
use libc::c_char;
use libc::size_t;
use ndarray::ArrayView1;
use quantity::si::*;
use serde_json::Value;
use std::ffi::CStr;
use std::slice;
use std::sync::Arc;

/// Return a pointer to equation of state given a json string
#[no_mangle]
pub extern "C" fn eos_from_json(json: *const c_char) -> *mut Arc<EosVariant> {
    let c_str = unsafe {
        assert!(!json.is_null());
        CStr::from_ptr(json)
    };
    let r_str = c_str.to_str().unwrap();

    let v: Value = serde_json::from_str(r_str).unwrap();
    let model_name = match &v["model"] {
        Value::Null => return std::ptr::null_mut(),
        Value::String(s) => s,
        _ => return std::ptr::null_mut(),
    };

    match model_name.to_lowercase().as_str() {
        "pc-saft" | "pcsaft" => {
            let pure_records: Vec<PureRecord<PcSaftRecord, JobackRecord>> =
                serde_json::from_value(v["substance_parameters"].clone()).unwrap();
            let binary_records: Vec<BinaryRecord<Identifier, PcSaftBinaryRecord>> =
                serde_json::from_value(v["binary_parameters"].clone()).unwrap();
            let binary_matrix = <PcSaftParameters as Parameter>::binary_matrix_from_records(
                &pure_records,
                &binary_records,
                IdentifierOption::Name,
            );
            let parameters = Arc::new(PcSaftParameters::from_records(pure_records, binary_matrix));
            return Box::into_raw(Box::new(Arc::new(EosVariant::PcSaft(PcSaft::new(
                parameters,
            )))));
        }
        _ => return std::ptr::null_mut(),
    }
}

#[no_mangle]
pub extern "C" fn free_eos_ptr(eos_ptr: *mut Arc<EosVariant>) {
    if eos_ptr.is_null() {
        return;
    }
    unsafe { drop(Box::from_raw(eos_ptr)) }
}

#[no_mangle]
pub extern "C" fn pressure_bar(
    eos_ptr: *mut Arc<EosVariant>,
    temperature: f64,
    density: f64,
    molefracs: *const f64,
    len: size_t,
) -> f64 {
    let x = ArrayView1::from(unsafe {
        assert!(!molefracs.is_null());
        slice::from_raw_parts(molefracs, len as usize)
    });
    let eos = unsafe {
        assert!(!eos_ptr.is_null());
        &*eos_ptr
    };
    assert!(eos.components() == len as usize);
    let pressure = StateBuilder::new(&eos)
        .temperature(temperature * KELVIN)
        .density(density * MOL / METER.powi(3))
        .molefracs(&x.to_owned())
        .build()
        .unwrap()
        .pressure(Contributions::Total);
    return pressure.to_reduced(BAR).unwrap();
}

#[no_mangle]
pub extern "C" fn da_dv(
    eos_ptr: *mut Arc<EosVariant>,
    temperature: f64,
    density: f64,
    molefracs: *const f64,
    len: size_t,
) -> f64 {
    let x = ArrayView1::from(unsafe {
        assert!(!molefracs.is_null());
        slice::from_raw_parts(molefracs, len as usize)
    });
    let eos = unsafe {
        assert!(!eos_ptr.is_null());
        &*eos_ptr
    };
    assert!(eos.components() == len as usize);
    let state = StateBuilder::new(&eos)
        .temperature(temperature * KELVIN)
        .density(density * MOL / METER.powi(3))
        .molefracs(&x.to_owned())
        .build()
        .unwrap();
    let state_hd = state.derive1(feos_core::Derivative::DV);
    let da_dv = eos.evaluate_residual(&state_hd);
    return da_dv.eps[0];
}
