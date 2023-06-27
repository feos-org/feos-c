use feos::pcsaft::*;
use feos::ResidualModel;
use feos_core::parameter::BinaryRecord;
use feos_core::parameter::Identifier;
use feos_core::parameter::IdentifierOption;
use feos_core::parameter::Parameter;
use feos_core::parameter::PureRecord;
use feos_core::Components;
use feos_core::Contributions;
use feos_core::Derivative;
use feos_core::EosUnit;
use feos_core::Residual;
use feos_core::State;
use feos_core::StateBuilder;
use libc::c_char;
use libc::size_t;
use ndarray::ArrayView1;
use quantity::si::*;
use serde_json::Value;
use std::ffi::CStr;
use std::ffi::CString;
use std::slice;
use std::sync::Arc;

/// Return a pointer to equation of state given a json string
#[no_mangle]
pub extern "C" fn eos_from_json(json: *const c_char) -> *mut Arc<ResidualModel> {
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
            let pure_records: Vec<PureRecord<PcSaftRecord>> =
                serde_json::from_value(v["substance_parameters"].clone()).unwrap();
            let binary_records: Vec<BinaryRecord<Identifier, PcSaftBinaryRecord>> =
                serde_json::from_value(v["binary_parameters"].clone()).unwrap();
            let binary_matrix = <PcSaftParameters as Parameter>::binary_matrix_from_records(
                &pure_records,
                &binary_records,
                IdentifierOption::Name,
            );
            let parameters = Arc::new(PcSaftParameters::from_records(pure_records, binary_matrix));
            return Box::into_raw(Box::new(Arc::new(ResidualModel::PcSaft(PcSaft::new(
                parameters,
            )))));
        }
        _ => return std::ptr::null_mut(),
    }
}

#[no_mangle]
pub extern "C" fn free_eos_ptr(eos_ptr: *mut Arc<ResidualModel>) {
    if eos_ptr.is_null() {
        return;
    }
    unsafe { drop(Box::from_raw(eos_ptr)) }
}

fn calculate_derivatives(state: &State<ResidualModel>, nt: i32, nd: i32) -> Result<f64, String> {
    let rho = state
        .density
        .to_reduced(SIUnit::reference_density())
        .unwrap();
    let t = state
        .temperature
        .to_reduced(SIUnit::reference_temperature())
        .unwrap();
    match (nt, nd) {
        (0, 0) => Ok(state.eos.evaluate_residual(&state.derive0())),
        (1, 0) => Ok(state
            .eos
            .evaluate_residual(&state.derive1(Derivative::DT)).eps
            * -t),
        (0, 1) => Ok(state
            .eos
            .evaluate_residual(&state.derive1(Derivative::DV))
            .eps / -rho),
        (1, 1) => Ok(state
            .eos
            .evaluate_residual(&state.derive2_mixed(Derivative::DT, Derivative::DV))
            .eps1eps2),
        (2, 0) => Ok(state
            .eos
            .evaluate_residual(&state.derive2(Derivative::DT))
            .v2),
        (0, 2) => Ok(state
            .eos
            .evaluate_residual(&state.derive2(Derivative::DV))
            .v2),
        _ => Err(format!(
            "Derivative of degree dt: {}, drho: {} is not implemented.\nHighest derivatives:\ndt:   2\ndrho: 2",
            nt, nd
        )),
    }
}

#[no_mangle]
pub extern "C" fn get_Arxy(
    eos_ptr: *mut Arc<ResidualModel>,
    nt: i32,
    nd: i32,
    t: f64,
    rho: f64,
    molefrac: *const f64,
    ncomp: i32,
    val: &mut f64,
    errmsg: *mut u8,
    errmsg_length: i32,
) -> i32 {
    let eos = unsafe {
        assert!(!eos_ptr.is_null());
        &*eos_ptr
    };
    let x = ArrayView1::from(unsafe {
        assert!(!molefrac.is_null());
        slice::from_raw_parts(molefrac, ncomp as usize)
    });

    let state = StateBuilder::new(&eos)
        .temperature(t * KELVIN)
        .density(rho * MOL / METER.powi(3))
        .molefracs(&x.to_owned())
        .build();

    match state {
        Ok(s) => match calculate_derivatives(&s, nt, nd) {
            Ok(v) => *val = v,
            Err(e) => unsafe {
                if errmsg.is_null() {
                    return -1;
                }
                let err_msg = e.to_string();
                let err_msg_c = CString::new(err_msg).unwrap();
                let bytes = err_msg_c.as_bytes_with_nul();
                let msg_bytes = slice::from_raw_parts_mut(errmsg, errmsg_length as usize);
                msg_bytes[..bytes.len()].copy_from_slice(bytes);
                return -1;
            },
        },
        Err(e) => {
            unsafe {
                if errmsg.is_null() {
                    return -1;
                }
                // https://stackoverflow.com/questions/51320714/what-is-the-correct-way-to-fill-a-c-string-pointer-from-rust
                let err_msg = e.to_string();
                let err_msg_c = CString::new(err_msg).unwrap();
                let bytes = err_msg_c.as_bytes_with_nul();
                let msg_bytes = slice::from_raw_parts_mut(errmsg, errmsg_length as usize);
                msg_bytes[..bytes.len()].copy_from_slice(bytes);
                return -1;
            };
        }
    }
    0
}

#[no_mangle]
pub extern "C" fn pressure_bar(
    eos_ptr: *const Arc<ResidualModel>,
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
    eos_ptr: *const Arc<ResidualModel>,
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
    return da_dv.eps;
}

#[no_mangle]
pub extern "C" fn reduced_entropy(
    eos_ptr: *const Arc<ResidualModel>,
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
    let residual = StateBuilder::new(&eos)
        .temperature(temperature * KELVIN)
        .density(density * MOL / METER.powi(3))
        .molefracs(&x.to_owned())
        .build()
        .unwrap()
        .residual_entropy();
    return residual.to_reduced(KB).unwrap();
}
