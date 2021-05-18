
function errs = get_default_atm_errs(d)

err_k = d.out_k.haf_err;
err_di = d.out_di.haf_err;
err_ecf = d.out_ecf.haf_err;
err_ecf_p = d.out_ecf_p.haf_err;
if isfield(d, 'out_hyb')
    err_ecfh = d.out_hyb.haf_err;
    err_ecfh_p = d.out_hyb_p.haf_err;
else
    err_ecfh = d.out_ecfh.haf_err;
    err_ecfh_p = d.out_ecfh_p.haf_err;
end

errs = [err_k, err_di, err_ecf, err_ecfh, err_ecf_p, err_ecfh_p];

end