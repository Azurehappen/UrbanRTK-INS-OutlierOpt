function Es = ephpos_err(cpt,p)

Es = NaN(length(cpt.corr_range),1);
for i = 1:length(cpt.corr_range)
    Es(i) = norm(p.P_base-cpt.s_pos_prc(:,i))+sagnac(p,cpt.s_pos_prc(:,i),p.P_base)...
        - (norm(p.P_base-cpt.s_pos_ecef(:,i))+sagnac(p,cpt.s_pos_ecef(:,i),p.P_base));
end

end