function obs = simplyobs(obs,count)
% Delete empty data in obs strcut

if count < length(obs.tr_sow)
    obs.tr_prime(:,count+1:end)=[]; 
    obs.tr_sow(:,count+1:end)=[];
    obs.tr_week(:,count+1:end)=[];
    obs.gps(1).data.P(:,count+1:end)=[];  obs.gps(2).data.P(:,count+1:end)=[];
    obs.gps(1).data.C(:,count+1:end)=[];  obs.gps(2).data.C(:,count+1:end)=[];
    obs.gps(1).data.D(:,count+1:end)=[];  obs.gps(2).data.D(:,count+1:end)=[];
    obs.gps(1).data.S(:,count+1:end)=[];  obs.gps(2).data.S(:,count+1:end)=[];
    obs.gps(3).data.P(:,count+1:end)=[];  obs.gps(4).data.P(:,count+1:end)=[];
    obs.gps(3).data.C(:,count+1:end)=[];  obs.gps(4).data.C(:,count+1:end)=[];
    obs.gps(3).data.D(:,count+1:end)=[];  obs.gps(4).data.D(:,count+1:end)=[];
    obs.gps(3).data.S(:,count+1:end)=[];  obs.gps(4).data.S(:,count+1:end)=[];
    obs.glo(1).data.P(:,count+1:end)=[];  obs.glo(2).data.P(:,count+1:end)=[];
    obs.glo(1).data.C(:,count+1:end)=[];  obs.glo(2).data.C(:,count+1:end)=[];
    obs.glo(1).data.D(:,count+1:end)=[];  obs.glo(2).data.D(:,count+1:end)=[];
    obs.glo(1).data.S(:,count+1:end)=[];  obs.glo(2).data.S(:,count+1:end)=[];
    obs.glo(3).data.P(:,count+1:end)=[];  obs.glo(4).data.P(:,count+1:end)=[];
    obs.glo(3).data.C(:,count+1:end)=[];  obs.glo(4).data.C(:,count+1:end)=[];
    obs.glo(3).data.D(:,count+1:end)=[];  obs.glo(4).data.D(:,count+1:end)=[];
    obs.glo(3).data.S(:,count+1:end)=[];  obs.glo(4).data.S(:,count+1:end)=[];
    obs.gal(1).data.P(:,count+1:end)=[];  obs.gal(2).data.P(:,count+1:end)=[];
    obs.gal(1).data.C(:,count+1:end)=[];  obs.gal(2).data.C(:,count+1:end)=[];
    obs.gal(1).data.D(:,count+1:end)=[];  obs.gal(2).data.D(:,count+1:end)=[];
    obs.gal(1).data.S(:,count+1:end)=[];  obs.gal(2).data.S(:,count+1:end)=[];
	obs.bds(1).data.P(:,count+1:end)=[];  obs.bds(2).data.P(:,count+1:end)=[];
    obs.bds(1).data.C(:,count+1:end)=[];  obs.bds(2).data.C(:,count+1:end)=[];
    obs.bds(1).data.D(:,count+1:end)=[];  obs.bds(2).data.D(:,count+1:end)=[];
    obs.bds(1).data.S(:,count+1:end)=[];  obs.bds(2).data.S(:,count+1:end)=[];
end

end