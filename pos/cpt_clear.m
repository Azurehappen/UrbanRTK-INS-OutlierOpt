function cpt = cpt_clear(cpt)
% Clear the data where don't have diff correction
i0 = find(cpt.diff_corr==0);
if ~isempty(i0)
    ind_prn = find(cpt.prn_record~=0);
    cpt.prn_record(ind_prn(i0))=0;
    cpt.svprn_mark(ind_prn(i0))=0;
    cpt.diff_corr(i0)=[];
    cpt.corr_range(i0)=[];
    cpt.doppler(i0)=[];
    cpt.s_pos_ecef(:,i0)=[];
    cpt.elev(i0) = [];
    cpt.az(i0) = [];
end

end