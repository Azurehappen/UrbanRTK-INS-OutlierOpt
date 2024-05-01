los = [-0.666281913791409	0.742465936998521	-0.0693592498067464
-0.288898123603725	0.921498863679181	-0.259570064661909
-0.607752643413810	0.0548124218039874	-0.792231461184353
0.187829277541013	0.723907067652485	0.663836544652036
-0.506624141949879	0.729335321713035	-0.459782946689689
0.694018274618246	-0.106036190130782	-0.712107325730084
0.631690513333651	0.746484128580154	0.209117073926706
-0.768058536079231	0.475688272632016	-0.428724316962579
0.681758965691222	0.491305490521296	0.542056512004035
0.550919025261321	0.355476569936566	-0.755066865670995
0.138203286218740	0.717046894805025	-0.683186439442803];

% double difference los
dd_los = los(1:end-1,:);
dd_los = dd_los - los(end,:);

%% Double difference approach
% It estimates the double-diff ambiguity N_i - N_p.
% We don't need information about N_p.
H_amb1 = 0.1903 * eye(size(dd_los,1));

H1 = [dd_los,zeros(size(H_amb1));  
    dd_los,H_amb1];

cov1 = (H1'*H1)^(-1);

%% Estimating clock approach
% It estimates each ambiguity.
H_amb2 = 0.1903 * eye(size(los,1));
H2 = [los,ones(size(los,1),1),zeros(size(los,1),1),zeros(size(H_amb2));
    los,ones(size(los,1),1),ones(size(los,1),1),H_amb2];
H2 = [los,ones(size(los,1),1),zeros(size(H_amb2));
    los,ones(size(los,1),1),H_amb2];
cov2 = (H2'*H2)^(-1);

%% DD for phase only
H_amb3 = 0.1903 * eye(size(dd_los,1));

H3 = [los,ones(size(los,1),1),zeros(size(los,1),size(H_amb3,2));  
    dd_los,zeros(size(dd_los,1),1),H_amb3];

cov3 = (H3'*H3)^(-1);